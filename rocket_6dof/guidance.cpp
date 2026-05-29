//  guidance.cpp  --  Time-programmed pitch guidance + LTG ascent guidance
//
//  Cross-reference to Zipfel guidance.cpp (Section 10.3.x):
//
//  LTG (Linear Tangent Guidance) is a closed-form optimal trajectory
//  solver for multi-stage powered ascent.  Given the desired end-state
//  (radius, speed, flight-path angle) and the per-stage propulsion
//  characteristics, it iteratively computes:
//    * tgo  -- time remaining until end-state is reached
//    * VGO  -- velocity still to be gained
//    * ULAM -- unit thrust vector to drive VGO to zero efficiently
//    * LAMD -- turning rate of the thrust vector
//    * UTIC -- final unit thrust vector command in inertial coords
//    * UTBC -- UTIC rotated into body coords (for display / future use)
//
//  Reference: Jaggers, R., "Multistage Linear Tangent Ascent Guidance
//  as Baselined for the Space Shuttle Vehicle," NASA MSC, Internal
//  Note MSC-EG-72-39, June 1972.
//
//  LTG runs as a parallel diagnostic computation.  Its outputs (UTBC
//  in particular) are not currently wired into Control's command path.

#include "guidance.h"
#include "control.h"
#include "newton.h"
#include "ins.h"
#include "propulsion.h"
#include "cad_utils.h"
#include <cmath>
#include <cstdio>

namespace rocket6dof {

namespace {
constexpr double RAD = osk::PI / 180.0;
constexpr double DEG = 180.0 / osk::PI;

// Helper: unit vector of v (returns zero vector if v has zero magnitude)
osk::Vec univec3(const osk::Vec& v) {
    double m = v.mag();
    if (m < 1e-30) return osk::Vec(0, 0, 0);
    return v * (1.0 / m);
}
} // anon

Guidance::Guidance()
    : control(nullptr), newton(nullptr), ins(nullptr), prop(nullptr),
      mguide(0),
      ltg_drives_control(0),
      phase2_yaw_sign(1.0),
      t_pitch_start(5.0), t_pitch_end(8.0),
      ancomx_program(0.5), alcomx_program(0.0),
      t_att_start(5.0), t_att_end(15.0),
      theta_com_start(89.0), theta_com_end(89.0),
      psi_com_start(90.0), psi_com_end(90.0),
      ancomx_cmd(0.0), alcomx_cmd(0.0), phase(0),
      ltg_step(0.5),
      dbi_desired(0.0), dvbi_desired(0.0), thtvdx_desired(0.0),
      num_stages(1),
      delay_ignition(0.0),
      amin(2.0),
      lamd_limit(0.1),
      time_ltg(0.0),
      init_flag(1), inisw_flag(1), skip_flag(1),
      ipas_flag(1), ipas2_flag(1), print_flag(1),
      beco_flag(0), ltg_count(0), nst(1),
      VGO(0,0,0), RBIAS(0,0,0), RGRAV(0,0,0), RGO(0,0,0),
      SDII(0,0,0), UD(0,0,0), UY(0,0,0), UZ(0,0,0),
      tgo(0.0), tgop(0.0),
      UTIC(0,0,0), UTBC(0,0,0), ULAM(0,0,0), LAMD(0,0,0),
      vgom(0.0), lamd(0.0),
      ddb(0.0), dvdb(0.0), thtvddbx(0.0),
      nstmax(0),
      SPII_save(0,0,0), VPII_save(0,0,0)
{
    // Zero per-stage arrays (compile-time-sized MAX_LTG_STAGES)
    for (int k = 0; k < MAX_LTG_STAGES; ++k) {
        char_time[k]      = 0.0;
        exhaust_vel[k]    = 0.0;
        burnout_epoch[k]  = 0.0;
        TAUN_save[k]      = 0.0;
    }
}

void Guidance::init() {
    if (initCount == 0) {
        ancomx_cmd = 0.0;
        alcomx_cmd = 0.0;
        phase      = 0;

        time_ltg   = 0.0;
        init_flag  = 1;
        inisw_flag = 1;
        skip_flag  = 1;
        ipas_flag  = 1;
        ipas2_flag = 1;
        print_flag = 1;
        beco_flag  = 0;
        ltg_count  = 0;
        nst        = 1;

        for (int k = 0; k < MAX_LTG_STAGES; ++k) {
            TAUN_save[k] = char_time[k];
        }
    }
}

void Guidance::update() {
    double t = osk::State::t;

    // ---- mguide=0: off ----
    if (mguide == 0) {
        ancomx_cmd = 0.0;
        alcomx_cmd = 0.0;
        phase      = 0;
        UTBC       = osk::Vec(0, 0, 0);
        UTIC       = osk::Vec(0, 0, 0);
        return;
    }

    // ---- mguide=1: time-programmed pitch profile ----
    if (mguide == 1) {
        if (t < t_pitch_start) {
            phase = 1;
            ancomx_cmd = 0.0;
            alcomx_cmd = 0.0;
        }
        else if (t < t_pitch_end) {
            phase = 2;
            ancomx_cmd = ancomx_program;
            alcomx_cmd = alcomx_program;
        }
        else {
            phase = 3;
            ancomx_cmd = 0.0;
            alcomx_cmd = 0.0;
        }
        if (control) {
            control->ancomx = ancomx_cmd;
            control->alcomx = alcomx_cmd;
        }
        return;
    }

    // ---- mguide=2: time-programmed attitude profile ----
    // Linearly ramps the commanded body attitude between two waypoints.
    // Outputs are written directly into Control->theta_com_inertial and
    // Control->psi_com_inertial, which the cascaded autopilot (maut=50
    // or 60) reads as its attitude target.
    if (mguide == 2) {
        double theta_cmd, psi_cmd;
        if (t < t_att_start) {
            phase = 1;
            theta_cmd = theta_com_start;
            psi_cmd   = psi_com_start;
        }
        else if (t < t_att_end) {
            phase = 2;
            double dt = t_att_end - t_att_start;
            double f = (dt > 1.0e-6) ? (t - t_att_start) / dt : 1.0;
            theta_cmd = theta_com_start + f * (theta_com_end - theta_com_start);
            psi_cmd   = psi_com_start   + f * (psi_com_end   - psi_com_start);
        }
        else {
            phase = 3;
            theta_cmd = theta_com_end;
            psi_cmd   = psi_com_end;
        }
        if (control) {
            control->theta_com_inertial = theta_cmd;
            control->psi_com_inertial   = psi_cmd;
        }
        return;
    }

    // ---- mguide=5: LTG ----
    if (mguide == 5) {
        // Start the LTG clock on first call
        if (init_flag) {
            init_flag = 0;
            time_ltg  = 0.0;
        } else {
            time_ltg += osk::State::dt;
        }

        // Run LTG only every ltg_step seconds (not every kernel sub-stage)
        if (time_ltg > ltg_step * ltg_count) {
            ltg_count++;
            UTIC = ltg();

            // Rotate UTIC into body frame for display
            if (ins) {
                UTBC = ins->TBIC * UTIC;
            } else if (newton) {
                // Fallback to truth TBI through Newton's stored kin pointer.
                // Newton doesn't expose TBI directly; use INS if available.
                UTBC = UTIC;  // identity fallback
            }
        }

        // ---- Phase-2 closed-loop: convert UTBC to ancomx/alcomx ----
        // Only when:
        //   * ltg_drives_control flag is set
        //   * LTG settling period (skip_flag) is over -- UTBC is non-zero
        //   * Propulsion is wired so we know thrust and mass
        // Conversion uses the steady-state relation between commanded
        // normal acceleration and body attitude angle:
        //   For body to pitch up by angle alpha, the autopilot must
        //   command ancomx such that ancomx * G0 = (T/m) * sin(alpha).
        //   If UTBC = (1, ey, ez), the desired pitch-up angle is -ez
        //   (since UTBC.z > 0 means desired direction below body-x,
        //   i.e., nose-down command).  So
        //     ancomx = -(T/m/G0) * UTBC.z
        //     alcomx =  yaw_sign * (T/m/G0) * UTBC.y
        // gnmax/gymax limiters in Control will clip any over-command.
        if (ltg_drives_control && skip_flag == 0 && prop && control) {
            const double G0 = 9.80675445;
            double m = prop->vmass;
            double T = prop->thrust;
            if (m > 1.0e-3 && T > 1.0) {
                double scale = T / (m * G0);
                ancomx_cmd = -scale * UTBC.z;
                alcomx_cmd =  phase2_yaw_sign * scale * UTBC.y;
                control->ancomx = ancomx_cmd;
                control->alcomx = alcomx_cmd;
            }
        }
    }
}

void Guidance::rpt() {
    if (osk::State::sample(1.0)) {
        if (mguide == 5) {
            std::printf("Guid t=%7.3f  LTG tgo=%7.3f vgom=%8.2f  "
                        "UTIC=(%+.3f, %+.3f, %+.3f)\n",
                        osk::State::t, tgo, vgom,
                        UTIC.x, UTIC.y, UTIC.z);
        } else if (mguide == 1) {
            const char* phase_name = "pre";
            switch (phase) {
                case 1: phase_name = "vertical"; break;
                case 2: phase_name = "pitch-over"; break;
                case 3: phase_name = "gravity-turn"; break;
            }
            std::printf("Guid t=%7.3f  phase=%-12s  ancomx=%+.3f g\n",
                        osk::State::t, phase_name, ancomx_cmd);
        } else if (mguide == 2) {
            const char* phase_name = "pre";
            switch (phase) {
                case 1: phase_name = "hold-start"; break;
                case 2: phase_name = "ramping";    break;
                case 3: phase_name = "hold-end";   break;
            }
            double tc = control ? control->theta_com_inertial : 0.0;
            double pc = control ? control->psi_com_inertial   : 0.0;
            std::printf("Guid t=%7.3f  phase=%-12s  theta_com=%+.2f  psi_com=%+.2f deg\n",
                        osk::State::t, phase_name, tc, pc);
        }
    }
}

//  ltg() -- main LTG step; mirrors Zipfel's guidance_ltg().
osk::Vec Guidance::ltg() {
    osk::Vec UTIC_local(0, 0, 0);

    // Once BECO has fired, freeze outputs and skip the algorithm.  The
    // arrays and state aren't designed for "past end of burn" operation
    // and continuing to call the sub-functions risks nst overrun and
    // numerical breakdown.
    if (beco_flag) {
        return UTIC;     // hold last good direction
    }

    // ---- Read INS-computed state (or truth if INS not available) ----
    osk::Vec SBIIC, VBIIC, FSPCB;
    osk::Mat TBIC;
    double dbi_in = 0.0, dvbi_in = 0.0, thtvdx_in = 0.0;
    if (ins) {
        SBIIC = ins->SBIIC;
        VBIIC = ins->VBIIC;
        FSPCB = ins->FSPCB;
        TBIC  = ins->TBIC;
    } else if (newton) {
        SBIIC = newton->SBII;
        VBIIC = newton->VBII;
        FSPCB = newton->FSPB;
        TBIC  = osk::Mat(1,0,0, 0,1,0, 0,0,1);
    }
    if (newton) {
        dbi_in    = newton->dbi;
        dvbi_in   = newton->dvbi;
        thtvdx_in = newton->thtvdx;
    }

    // Current vehicle acceleration in inertial frame: ABII = ~TBIC * FSPCB
    osk::Vec ABII = TBIC.transpose() * FSPCB;
    double   amag1 = ABII.mag();

    // ---- Initialization on first cycle ----
    osk::Vec VMISS(0, 0, 0);
    osk::Vec SPII = SPII_save;
    osk::Vec VPII = VPII_save;

    if (inisw_flag) {
        inisw_flag = 0;
        SPII = SBIIC;
        VPII = VBIIC;
        ltg_crct(SDII, UD, UY, UZ, VMISS, VGO,
                 dbi_desired, dvbi_desired, thtvdx_desired,
                 SPII, VPII, SBIIC, VBIIC);
    } else {
        // Update velocity-to-go: VGO -= ABII * ltg_step
        VGO = VGO - ABII * ltg_step;
    }

    vgom = VGO.mag();

    // ---- Build per-stage data vectors (sized MAX_LTG_STAGES) ----
    double TAUN[MAX_LTG_STAGES];
    double VEXN[MAX_LTG_STAGES];
    double BOTN[MAX_LTG_STAGES + 1];  // BOTN[0] is unused per Zipfel indexing
    BOTN[0] = 0.0;
    for (int k = 0; k < MAX_LTG_STAGES; ++k) {
        TAUN[k]     = TAUN_save[k];
        VEXN[k]     = exhaust_vel[k];
        BOTN[k + 1] = burnout_epoch[k];
    }

    // ---- _tgo(): compute time-to-go, burn times, l_igrl, nstmax ----
    double tgop_loc = tgo;
    double BURNTN[MAX_LTG_STAGES] = {0};
    double L_IGRLN[MAX_LTG_STAGES] = {0};
    double TGON[MAX_LTG_STAGES] = {0};
    double l_igrl = 0.0;
    int    nstmax_loc = 0;
    ltg_tgo(tgop_loc, BURNTN, L_IGRLN, TGON, l_igrl, nstmax_loc,
            tgo, nst, TAUN, VEXN, BOTN,
            delay_ignition, vgom, amag1, amin, time_ltg, num_stages);
    tgop   = tgop_loc;
    nstmax = nstmax_loc;
    // Save mutated TAUN back
    for (int k = 0; k < MAX_LTG_STAGES; ++k) {
        TAUN_save[k] = TAUN[k];
    }

    // ---- _igrl(): five thrust integrals ----
    double s_igrl = 0, j_igrl = 0, q_igrl = 0, h_igrl = 0, p_igrl = 0;
    double j_over_l = 0, tlam = 0, qprime = 0;
    ltg_igrl(s_igrl, j_igrl, q_igrl, h_igrl, p_igrl, j_over_l, tlam, qprime,
             nst, nstmax, BURNTN, L_IGRLN, TGON, TAUN, VEXN, l_igrl, time_ltg);

    // ---- _trate(): turning rate ULAM, LAMD, RGO ----
    ltg_trate(ULAM, LAMD, RGO, ipas2_flag,
              VGO, s_igrl, q_igrl, j_over_l, lamd_limit, vgom, time_ltg,
              tgo, tgop, SDII, SBIIC, VBIIC, RBIAS, UD, UY, UZ, RGRAV);
    lamd = LAMD.mag();

    // ---- Compute the thrust command vector ----
    osk::Vec TC = ULAM + LAMD * (time_ltg - tlam);

    // Settling: discard first 10 cycles
    if (skip_flag) {
        skip_flag++;
        if (skip_flag == 10) {
            skip_flag = 0;
        }
    } else {
        UTIC_local = univec3(TC);
    }

    // ---- _pdct(): predict end-state ----
    ltg_pdct(SPII, VPII, RGRAV, RBIAS,
             LAMD, ULAM, l_igrl, s_igrl, j_igrl, q_igrl, h_igrl, p_igrl,
             j_over_l, qprime, SBIIC, VBIIC, RGO, tgo);
    SPII_save = SPII;
    VPII_save = VPII;

    // ---- _crct(): correct VGO based on predicted vs desired ----
    ltg_crct(SDII, UD, UY, UZ, VMISS, VGO,
             dbi_desired, dvbi_desired, thtvdx_desired,
             SPII, VPII, SBIIC, VBIIC);

    // ---- BECO logic ----
    if (tgo < 10.0 * osk::State::dt) {
        if (!beco_flag) {
            beco_flag = 1;
            if (print_flag) {
                print_flag = 0;
                ddb      = dbi_desired - dbi_in;
                dvdb     = dvbi_desired - dvbi_in;
                thtvddbx = thtvdx_desired - thtvdx_in;
                std::printf("\n *** LTG Boost engine cut-off: t=%.3f sec ***\n",
                            osk::State::t);
                std::printf("     Orbital position dbi = %.1f m    "
                            "Inertial speed dvbi = %.2f m/s    "
                            "FPA thtvdx = %.3f deg\n",
                            dbi_in, dvbi_in, thtvdx_in);
                std::printf("     Position error    ddb = %+.2f m    "
                            "Speed error    dvdb = %+.3f m/s    "
                            "FPA error  thtvddbx = %+.4f deg\n\n",
                            ddb, dvdb, thtvddbx);
            }
        }
    }

    return UTIC_local;
}

//  ltg_tgo  --  Time-to-go to desired end-state
void Guidance::ltg_tgo(double& tgop_out,
                       double BURNTN[MAX_LTG_STAGES], double L_IGRLN[MAX_LTG_STAGES],
                       double TGON[MAX_LTG_STAGES], double& l_igrl, int& nstmax_out,
                       double& tgo_io, int& nst_io,
                       double TAUN[MAX_LTG_STAGES],
                       const double VEXN[MAX_LTG_STAGES], const double BOTN[MAX_LTG_STAGES+1],
                       double delay_ign, double vgom_in, double amag1,
                       double amin_in, double time_ltg_in, int num_stages_in) {
    if (ipas_flag) {
        nst_io = 1;
    }
    tgop_out = tgo_io;

    tgo_io = 0.0;
    l_igrl = 0.0;
    nstmax_out = num_stages_in;

    // Guard: don't advance past last stage (BOTN array is size num_stages+1)
    if (nst_io < num_stages_in
        && time_ltg_in >= BOTN[nst_io]) {
        nst_io++;
    }
    // Hard cap to keep nst_io in [1, num_stages]
    if (nst_io > num_stages_in) {
        nst_io = num_stages_in;
    }

    int i = nst_io - 1;
    for (i = nst_io - 1; i < nstmax_out; i++) {
        if (i == nst_io - 1) {
            TAUN[nst_io - 1] = TAUN[nst_io - 1] - (time_ltg_in - BOTN[nst_io - 1]);
        }
        if ((amag1 >= amin_in)
            && (time_ltg_in > (BOTN[nst_io - 1] + delay_ign))) {
            TAUN[nst_io - 1] = VEXN[nst_io - 1] / amag1;
        }
        if (i == nst_io - 1) {
            BURNTN[i] = BOTN[i + 1] - time_ltg_in;
        } else {
            BURNTN[i] = BOTN[i + 1] - BOTN[i];
        }

        double dum1 = TAUN[i];
        double dum2 = BURNTN[i];
        L_IGRLN[i] = -VEXN[i] * std::log(1.0 - dum2 / dum1);

        l_igrl += L_IGRLN[i];
        if (l_igrl < vgom_in) {
            tgo_io += BURNTN[i];
            TGON[i] = tgo_io;
        } else {
            i++;
            break;
        }
    }
    nstmax_out = i;

    // Adjust last stage to exactly match vgom
    if (i > 0) {
        l_igrl -= L_IGRLN[i - 1];
        double almx = vgom_in - l_igrl;
        L_IGRLN[i - 1] = almx;
        double dum3 = VEXN[i - 1];
        BURNTN[i - 1] = TAUN[i - 1] * (1.0 - std::exp(-almx / dum3));
        // recompute tgo properly: subtract old burnt and add new
        // Zipfel just accumulates: tgo += BURNTN[i-1].  But we've also
        // already accumulated the OLD BURNTN[i-1] into tgo before
        // entering the else-break.  Zipfel doesn't undo it because for
        // i==nstmax (only the LAST stage), the accumulation never
        // happened (the break occurred BEFORE the accumulation).
        // We replicate the same logic: tgo_io was last updated for
        // i < (nstmax-1) stages; the (nstmax-1)-th stage burn time
        // is now being computed and added.
        tgo_io += BURNTN[i - 1];
        TGON[i - 1] = tgo_io;
    }
    l_igrl = vgom_in;

    if (ipas_flag) {
        tgop_out = tgo_io;
        ipas_flag = 0;
    }
}

//  ltg_igrl  --  Five thrust integrals over remaining boost phase
void Guidance::ltg_igrl(double& s_igrl, double& j_igrl, double& q_igrl,
                        double& h_igrl, double& p_igrl, double& j_over_l,
                        double& tlam, double& qprime,
                        int nst_in, int nstmax_in,
                        const double BURNTN[MAX_LTG_STAGES], const double L_IGRLN[MAX_LTG_STAGES],
                        const double TGON[MAX_LTG_STAGES], const double TAUN[MAX_LTG_STAGES],
                        const double VEXN[MAX_LTG_STAGES], double l_igrl, double time_ltg_in) {
    s_igrl = j_igrl = q_igrl = h_igrl = p_igrl = 0.0;
    double ls_igrl = 0.0;

    for (int i = nst_in - 1; i < nstmax_in; i++) {
        double tb  = BURNTN[i];
        double tga = TGON[i];
        double dummy = TAUN[i];
        double x = tb / dummy;

        double a1 = 0.0;
        if (std::fabs(x - 2.0) < 1e-12) {
            a1 = 1.0 / (1.0 - 0.5 * x * 1.001);
        } else {
            a1 = 1.0 / (1.0 - 0.5 * x);
        }
        double a2 = 0.0;
        if (std::fabs(x - 1.0) < 1e-12) {
            std::fprintf(stderr, " *** LTG: end-state cannot be reached *** \n");
            return;
        } else {
            a2 = 1.0 / (1.0 - x);
        }
        double aa = VEXN[i] / dummy;   // longitudinal accel
        double ll_igrl = L_IGRLN[i];
        double a1x   = 4.0 * a1 - a2 - 3.0;
        double a2xsq = 2.0 * a2 - 4.0 * a1 + 2.0;
        double sa = (aa * tb * tb / 2.0)         * (1.0 + a1x / 3.0 + a2xsq / 6.0);
        double ha = (aa * tb * tb * tb / 3.0)    * (1.0 + a1x * 0.75 + a2xsq * 0.6);
        double ja = (aa * tb * tb / 2.0)         * (1.0 + a1x * (2.0 / 3.0) + a2xsq / 2.0);
        double qa = (aa * tb * tb * tb / 6.0)    * (1.0 + a1x / 2.0 + a2xsq * 0.3);
        double pa = (aa * tb * tb * tb * tb / 12.0) * (1.0 + a1x * 0.6 + a2xsq * 0.4);

        if (i != nst_in - 1) {
            double t1 = TGON[i - 1];
            ha = ha + 2.0 * t1 * ja + t1 * t1 * ll_igrl;
            ja = ja + t1 * ll_igrl;
            pa = pa + 2.0 * t1 * qa + t1 * t1 * sa;
            qa = qa + t1 * sa;
        }
        ha = ja * tga - qa;
        sa = sa + ls_igrl * tb;
        qa = qa + j_igrl * tb;
        pa = pa + h_igrl * tb;

        s_igrl += sa;
        q_igrl += qa;
        p_igrl += pa;
        h_igrl += ha;
        ls_igrl += ll_igrl;
        j_igrl += ja;
    }
    if (std::fabs(l_igrl) > 1e-30) {
        j_over_l = j_igrl / l_igrl;
    } else {
        j_over_l = 0.0;
    }
    tlam   = time_ltg_in + j_over_l;
    qprime = q_igrl - s_igrl * j_over_l;
}

//  ltg_trate  --  Turning rate ULAM, LAMD and range-to-go RGO
void Guidance::ltg_trate(osk::Vec& ULAM_out, osk::Vec& LAMD_out,
                         osk::Vec& RGO_out, int& ipas2_io,
                         const osk::Vec& VGO_in, double s_igrl,
                         double q_igrl, double j_over_l, double lamd_lim,
                         double vgom_in, double /*time_ltg_in*/,
                         double tgo_in, double tgop_in,
                         const osk::Vec& SDII_in, const osk::Vec& SBIIC,
                         const osk::Vec& VBIIC, const osk::Vec& RBIAS_in,
                         const osk::Vec& UD_in, const osk::Vec& UY_in,
                         const osk::Vec& UZ_in, osk::Vec& RGRAV_io) {
    if (vgom_in == 0.0) return;

    ULAM_out = univec3(VGO_in);
    LAMD_out = osk::Vec(0, 0, 0);

    if (ipas2_io) {
        ipas2_io = 0;
        RGO_out = ULAM_out * s_igrl;
    }
    ltg_trate_rtgo(RGO_out, RGRAV_io, tgo_in, tgop_in, SDII_in, SBIIC, VBIIC,
                   RBIAS_in, ULAM_out, UD_in, UY_in, UZ_in, s_igrl);

    double denom = q_igrl - s_igrl * j_over_l;
    if (std::fabs(denom) > 1e-30) {
        LAMD_out = (RGO_out - ULAM_out * s_igrl) * (1.0 / denom);
    }

    double lamd_mag = LAMD_out.mag();
    if (lamd_mag >= lamd_lim && lamd_mag > 1e-30) {
        LAMD_out = LAMD_out * (lamd_lim / lamd_mag);
    }
}

void Guidance::ltg_trate_rtgo(osk::Vec& RGO_out, osk::Vec& RGRAV_io,
                              double tgo_in, double tgop_in,
                              const osk::Vec& SDII_in, const osk::Vec& SBIIC,
                              const osk::Vec& VBIIC, const osk::Vec& RBIAS_in,
                              const osk::Vec& ULAM_in,
                              const osk::Vec& UD_in, const osk::Vec& UY_in,
                              const osk::Vec& UZ_in, double s_igrl) {
    // Scale gravity contribution
    if (tgop_in > 1e-30) {
        double r = tgo_in / tgop_in;
        RGRAV_io = RGRAV_io * (r * r);
    }
    osk::Vec RGO_LOCAL = SDII_in - (SBIIC + VBIIC * tgo_in + RGRAV_io) - RBIAS_in;

    double rgox = RGO_LOCAL.dot(UD_in);
    double rgoy = RGO_LOCAL.dot(UY_in);
    osk::Vec RGOXY = UD_in * rgox + UY_in * rgoy;

    double num   = RGOXY.dot(ULAM_in);
    double denom = ULAM_in.dot(UZ_in);
    if (std::fabs(denom) < 1e-30) {
        return;   // keep previous RGO_out
    }
    double rgoz = (s_igrl - num) / denom;
    RGO_out = RGOXY + UZ_in * rgoz;
}

//  ltg_pdct  --  Predict end-state via Keplerian + thrust corrections
void Guidance::ltg_pdct(osk::Vec& SPII_out, osk::Vec& VPII_out,
                        osk::Vec& RGRAV_io, osk::Vec& RBIAS_out,
                        const osk::Vec& LAMD_in, const osk::Vec& ULAM_in,
                        double l_igrl, double s_igrl, double j_igrl,
                        double q_igrl, double h_igrl, double p_igrl,
                        double j_over_l, double qprime,
                        const osk::Vec& SBIIC, const osk::Vec& VBIIC,
                        const osk::Vec& RGO_in, double tgo_in) {
    double lmdsq = LAMD_in.dot(LAMD_in);

    osk::Vec VTHRUST = ULAM_in * (l_igrl - 0.5 * lmdsq * (h_igrl - j_igrl * j_over_l));
    osk::Vec RTHRUST = ULAM_in * (s_igrl - 0.5 * lmdsq * (p_igrl - j_over_l * (q_igrl + qprime)))
                     + LAMD_in * qprime;

    RBIAS_out = RGO_in - RTHRUST;

    osk::Vec SBIIC1 = SBIIC - RTHRUST * 0.1 - VTHRUST * (tgo_in / 30.0);
    osk::Vec VBIIC1 = VBIIC + RTHRUST * (1.2 / std::max(tgo_in, 1e-30)) - VTHRUST * 0.1;

    osk::Vec SBIIC2(0,0,0), VBIIC2(0,0,0);
    int flag = cad_kepler(SBIIC2, VBIIC2, SBIIC1, VBIIC1, tgo_in);
    if (flag) {
        // Kepler failed -- keep best-effort estimate
        std::fprintf(stderr, " *** Warning: cad_kepler did not converge in LTG _pdct ***\n");
        SBIIC2 = SBIIC1 + VBIIC1 * tgo_in;
        VBIIC2 = VBIIC1;
    }
    osk::Vec VGRAV = VBIIC2 - VBIIC1;
    RGRAV_io = SBIIC2 - SBIIC1 - VBIIC1 * tgo_in;

    SPII_out = SBIIC + VBIIC * tgo_in + RGRAV_io + RTHRUST;
    VPII_out = VBIIC + VGRAV + VTHRUST;
}

//  ltg_crct  --  End-state corrector: update VGO to drive miss to zero
void Guidance::ltg_crct(osk::Vec& SDII_out, osk::Vec& UD_out,
                        osk::Vec& UY_out, osk::Vec& UZ_out,
                        osk::Vec& VMISS_out, osk::Vec& VGO_io,
                        double dbi_des, double dvbi_des, double thtvdx_des,
                        const osk::Vec& SPII_in, const osk::Vec& VPII_in,
                        const osk::Vec& SBIIC, const osk::Vec& VBIIC) {
    // Desired position lies along the predicted-position direction
    UD_out = univec3(SPII_in);
    SDII_out = UD_out * dbi_des;

    // Trajectory-plane base vectors
    UY_out = univec3(VBIIC.cross(SBIIC));
    UZ_out = univec3(UD_out.cross(UY_out));

    // Velocity-to-be-gained
    double th = thtvdx_des * RAD;
    osk::Vec VDII = (UD_out * std::sin(th) + UZ_out * std::cos(th)) * dvbi_des;
    VMISS_out = VPII_in - VDII;
    VGO_io = VGO_io - VMISS_out;
}

} // namespace rocket6dof
