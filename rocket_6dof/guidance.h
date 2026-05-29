//  guidance.h  --  Sounding-rocket / orbital guidance
//
//  Public interface matches Zipfel's guidance.cpp (Section 10.3.x,
//  hyper[400-499]).  Modes:
//
//    mguide = 0  no guidance (Control commands stay at whatever value
//                they're set to externally)
//    mguide = 1  time-programmed pitch profile (sounding-rocket flight
//                program)
//    mguide = 5  Linear Tangent Guidance (LTG) for orbital ascent.
//                Computes the optimal unit thrust vector UTIC (inertial
//                coords) and UTBC (body coords) to reach a desired
//                end-state (radius dbi_desired, speed dvbi_desired,
//                flight-path angle thtvdx_desired).
//                Reference: Jaggers, R., "Multistage Linear Tangent
//                Ascent Guidance as Baselined for the Space Shuttle
//                Vehicle," NASA MSC, Internal Note MSC-EG-72-39, 1972.
//
//  LTG in v1 (Phase 1):
//    * Computes UTBC, UTIC, tgo, vgom each LTG cycle (default 1 Hz)
//    * Does NOT yet drive Control's ancomx/alcomx -- the outputs are
//      diagnostic only, like in Zipfel's SSR6 example.  Closed-loop
//      LTG-to-Control wiring is Phase 2.
//    * Triggers a print at "boost engine cut-off" (BECO) when tgo
//      drops below ~10 integration steps.

#ifndef ROCKET6DOF_GUIDANCE_H
#define ROCKET6DOF_GUIDANCE_H

#include "../osk/osk.h"

namespace rocket6dof {

class Control;
class Newton;
class INS;
class Propulsion;

class Guidance : public osk::Block {
public:
    // ---- Inputs ----
    Control*    control;
    Newton*     newton;     // truth state (for diagnostics: dbi, dvbi, thtvdx)
    INS*        ins;        // INS estimates (LTG runs on these in real ops)
    Propulsion* prop;       // for thrust, mass (used by Phase-2 UTBC->ancomx)
    void getsFrom(Control* c) {
        control = c; newton = nullptr; ins = nullptr; prop = nullptr;
    }
    void getsFrom(Control* c, Newton* n, INS* i) {
        control = c; newton = n; ins = i; prop = nullptr;
    }
    void getsFrom(Control* c, Newton* n, INS* i, Propulsion* p) {
        control = c; newton = n; ins = i; prop = p;
    }

    // ---- Mode selector ----
    int mguide;             // 0 = off, 1 = time-programmed accel,
                            // 2 = time-programmed attitude (writes
                            //     theta_com_inertial/psi_com_inertial
                            //     into Control), 5 = LTG

    // ---- Phase-2 LTG closed-loop flag ----
    // When 1 and mguide==5, the LTG output UTBC is converted into
    // ancomx/alcomx and written into Control.  When 0 (default), LTG
    // remains diagnostic only (Phase 1 behavior).
    int ltg_drives_control;

    // ---- Phase-2 converter gain ----
    // The UTBC->accel conversion uses the steady-state relation
    //   ancomx = -(T/m/g) * UTBC.z
    //   alcomx = +(T/m/g) * UTBC.y * yaw_sign
    // Where yaw_sign is +1 by default (set via this parameter for testing).
    double phase2_yaw_sign;

    // ---- Time-program (mguide=1) parameters ----
    double t_pitch_start;
    double t_pitch_end;
    double ancomx_program;
    double alcomx_program;

    // ---- Attitude-program (mguide=2) parameters ----
    // The attitude commands ramp linearly between two waypoints during
    // the interval [t_att_start, t_att_end].  Outside that interval the
    // commands are held at theta_com_start/end or psi_com_start/end as
    // appropriate.
    double t_att_start;     // [s]   start of attitude ramp
    double t_att_end;       // [s]   end of attitude ramp
    double theta_com_start; // [deg] commanded body pitch at t<=t_att_start
    double theta_com_end;   // [deg] commanded body pitch at t>=t_att_end
    double psi_com_start;   // [deg] commanded body yaw at t<=t_att_start
    double psi_com_end;     // [deg] commanded body yaw at t>=t_att_end

    // ---- Time-program outputs ----
    double ancomx_cmd;
    double alcomx_cmd;
    int    phase;

    // ---- LTG (mguide=5) parameters ----
    double ltg_step;        // [s]   LTG cycle time, default 0.5 s
    double dbi_desired;     // [m]   desired radius from Earth center
    double dvbi_desired;    // [m/s] desired inertial speed
    double thtvdx_desired;  // [deg] desired flight-path angle
    int    num_stages;      // number of remaining boost stages (1..MAX_LTG_STAGES)
    double delay_ignition;  // [s]   ignition delay after staging
    double amin;            // [m/s^2] min accel below which tau is not updated

    // Per-stage propulsion characteristics (Zipfel's char_time = m0/m_dot).
    // Use the same compile-time MAX_STAGES as Propulsion (default 4) so
    // LTG can handle any vehicle Propulsion can model.  Backward-compat
    // accessors char_time1/2/3 etc. are still loaded by mission_config.
#ifdef ROCKET6DOF_MAX_STAGES
    static constexpr int MAX_LTG_STAGES = ROCKET6DOF_MAX_STAGES;
#else
    static constexpr int MAX_LTG_STAGES = 4;
#endif
    double char_time[MAX_LTG_STAGES];      // [s] m0/m_dot per stage
    double exhaust_vel[MAX_LTG_STAGES];    // [m/s] effective exhaust velocity per stage
    double burnout_epoch[MAX_LTG_STAGES];  // [s] time-from-LTG-start of stage burnout
    double lamd_limit;      // [1/s]   limit on |LAMD|

    // ---- LTG state (saved across calls) ----
    double time_ltg;        // [s]   time since LTG clock started
    int    init_flag;       // 1 = need to start clock
    int    inisw_flag;      // 1 = need to initialize SPII/VPII/etc.
    int    skip_flag;       // 1..10 = output suppression for transient
    int    ipas_flag;       // 1 = first call to _tgo()
    int    ipas2_flag;      // 1 = first call to _trate()
    int    print_flag;      // 1 = BECO printout still pending
    int    beco_flag;       // 1 = boost engine cut-off has fired
    int    ltg_count;       // counter of LTG cycles
    int    nst;             // current stage (1, 2, or 3)

    osk::Vec VGO;           // velocity to be gained [m/s]
    osk::Vec RBIAS;         // position bias [m]
    osk::Vec RGRAV;         // gravity position loss [m]
    osk::Vec RGO;           // range-to-go [m]
    osk::Vec SDII;          // desired inertial position [m]
    osk::Vec UD, UY, UZ;    // unit base vectors of the trajectory frame
    double   tgo;           // [s] time-to-go
    double   tgop;          // [s] previous time-to-go

    // ---- LTG diagnostics outputs ----
    osk::Vec UTIC;          // unit thrust vector, inertial coords
    osk::Vec UTBC;          // unit thrust vector, body coords
    osk::Vec ULAM;          // unit thrust along VGO
    osk::Vec LAMD;          // turning rate vector [1/s]
    double   vgom;          // |VGO|
    double   lamd;          // |LAMD|
    double   ddb;           // position error at BECO [m]
    double   dvdb;          // speed error at BECO [m/s]
    double   thtvddbx;      // FPA error at BECO [deg]
    int      nstmax;        // # of stages needed to meet end-state

    Guidance();
    void init()   override;
    void update() override;
    void rpt()    override;

    ACCESS_FN(double,   ancomx_cmd)
    ACCESS_FN(double,   alcomx_cmd)
    ACCESS_FN(osk::Vec, UTBC)
    ACCESS_FN(osk::Vec, UTIC)
    ACCESS_FN(double,   tgo)
    ACCESS_FN(double,   vgom)

private:
    // ---- LTG sub-functions (ported from Zipfel) ----
    osk::Vec ltg();   // main LTG step; returns UTIC

    // Returns nothing; updates l_igrl, nstmax, BURNTN[MAX_LTG_STAGES], L_IGRLN[MAX_LTG_STAGES],
    // TGON[MAX_LTG_STAGES], TAUN[MAX_LTG_STAGES] (in/out), tgo (in/out), nst (in/out).
    void ltg_tgo(double& tgop_out,
                 double BURNTN[MAX_LTG_STAGES], double L_IGRLN[MAX_LTG_STAGES],
                 double TGON[MAX_LTG_STAGES], double& l_igrl, int& nstmax,
                 double& tgo_io, int& nst_io,
                 double TAUN[MAX_LTG_STAGES],
                 const double VEXN[MAX_LTG_STAGES], const double BOTN[MAX_LTG_STAGES+1],
                 double delay_ign, double vgom_in, double amag1,
                 double amin_in, double time_ltg_in, int num_stages_in);

    void ltg_igrl(double& s_igrl, double& j_igrl, double& q_igrl,
                  double& h_igrl, double& p_igrl, double& j_over_l,
                  double& tlam, double& qprime,
                  int nst_in, int nstmax_in,
                  const double BURNTN[MAX_LTG_STAGES], const double L_IGRLN[MAX_LTG_STAGES],
                  const double TGON[MAX_LTG_STAGES], const double TAUN[MAX_LTG_STAGES],
                  const double VEXN[MAX_LTG_STAGES], double l_igrl, double time_ltg_in);

    void ltg_trate(osk::Vec& ULAM_out, osk::Vec& LAMD_out, osk::Vec& RGO_out,
                   int& ipas2_io,
                   const osk::Vec& VGO_in, double s_igrl, double q_igrl,
                   double j_over_l, double lamd_lim, double vgom_in,
                   double time_ltg_in, double tgo_in, double tgop_in,
                   const osk::Vec& SDII_in, const osk::Vec& SBIIC,
                   const osk::Vec& VBIIC, const osk::Vec& RBIAS_in,
                   const osk::Vec& UD_in, const osk::Vec& UY_in,
                   const osk::Vec& UZ_in, osk::Vec& RGRAV_io);

    void ltg_trate_rtgo(osk::Vec& RGO_out, osk::Vec& RGRAV_io,
                        double tgo_in, double tgop_in,
                        const osk::Vec& SDII_in, const osk::Vec& SBIIC,
                        const osk::Vec& VBIIC, const osk::Vec& RBIAS_in,
                        const osk::Vec& ULAM_in,
                        const osk::Vec& UD_in, const osk::Vec& UY_in,
                        const osk::Vec& UZ_in, double s_igrl);

    void ltg_pdct(osk::Vec& SPII_out, osk::Vec& VPII_out,
                  osk::Vec& RGRAV_io, osk::Vec& RBIAS_out,
                  const osk::Vec& LAMD_in, const osk::Vec& ULAM_in,
                  double l_igrl, double s_igrl, double j_igrl,
                  double q_igrl, double h_igrl, double p_igrl,
                  double j_over_l, double qprime,
                  const osk::Vec& SBIIC, const osk::Vec& VBIIC,
                  const osk::Vec& RGO_in, double tgo_in);

    void ltg_crct(osk::Vec& SDII_out, osk::Vec& UD_out, osk::Vec& UY_out,
                  osk::Vec& UZ_out, osk::Vec& VMISS_out,
                  osk::Vec& VGO_io,
                  double dbi_des, double dvbi_des, double thtvdx_des,
                  const osk::Vec& SPII_in, const osk::Vec& VPII_in,
                  const osk::Vec& SBIIC, const osk::Vec& VBIIC);

    // Saved state from prior LTG cycles (SPII/VPII from _pdct)
    osk::Vec SPII_save, VPII_save;
    // Saved characteristic times TAUN (mutated by _tgo)
    double TAUN_save[MAX_LTG_STAGES];
};

} // namespace rocket6dof

#endif
