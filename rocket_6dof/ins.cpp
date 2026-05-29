//  ins.cpp  --  Inertial Navigation System block
//
//  Cross-reference: Zipfel ins.cpp (Section 10.3.x, hyper[300-399]).
//
//  Modes:
//    0 (ideal): outputs equal truth.  This is the default and is what
//               existing tests assume.
//    1 (simple-error): adds constant biases to accel and gyro measurements;
//               propagates INS's own SBIIC, VBIIC, TBIC by integrating
//               the BIASED measurements.  Drifts over time, as a real
//               INS does between GPS / StarTrack updates.
//
//  When mins=1, the INS holds its own integrator states for position,
//  velocity, and the 9 components of TBIC.  Initial values come from
//  truth at t=0 (assumes "aligned" INS at startup), then the integration
//  diverges from truth at a rate proportional to the bias magnitudes.

#include "ins.h"
#include "newton.h"
#include "euler.h"
#include "kinematics.h"
#include "gps.h"
#include "startrack.h"
#include <cmath>
#include <cstdio>

namespace rocket6dof {

namespace {
constexpr double DEG = 180.0 / osk::PI;

osk::Mat skew(const osk::Vec& w) {
    return osk::Mat(
         0.0, -w.z,  w.y,
         w.z,  0.0, -w.x,
        -w.y,  w.x,  0.0
    );
}

osk::Mat get_mat(double a, double b, double c,
                 double d, double e, double f,
                 double g, double h, double i) {
    return osk::Mat(a,b,c, d,e,f, g,h,i);
}
} // anon

INS::INS()
    : newton(nullptr), euler(nullptr), kin(nullptr), gps(nullptr), startrack(nullptr),
      mins(0),
      gps_update_count(0),
      startrack_update_count(0),
      bias_accel(0, 0, 0),
      bias_gyro(0, 0, 0),
      SBIIC(0, 0, 0), VBIIC(0, 0, 0),
      TBIC(1,0,0, 0,1,0, 0,0,1),
      sxc(0), syc(0), szc(0),
      sxcd(0), sycd(0), szcd(0),
      vxc(0), vyc(0), vzc(0),
      vxcd(0), vycd(0), vzcd(0),
      t00c(1), t01c(0), t02c(0),
      t10c(0), t11c(1), t12c(0),
      t20c(0), t21c(0), t22c(1),
      t00cd(0), t01cd(0), t02cd(0),
      t10cd(0), t11cd(0), t12cd(0),
      t20cd(0), t21cd(0), t22cd(0),
      FSPCB(0, 0, 0), WBICB(0, 0, 0), WBICI(0, 0, 0),
      ppcx(0), qqcx(0), rrcx(0),
      dbic(0), dvbec(0),
      ins_pos_err(0), ins_vel_err(0), ins_att_err(0)
{
    // Register integrator states (only relevant when mins=1).  Even in
    // mins=0 mode they are registered but stay at zero derivative.
    addIntegrator(sxc, sxcd);
    addIntegrator(syc, sycd);
    addIntegrator(szc, szcd);
    addIntegrator(vxc, vxcd);
    addIntegrator(vyc, vycd);
    addIntegrator(vzc, vzcd);
    addIntegrator(t00c, t00cd);
    addIntegrator(t01c, t01cd);
    addIntegrator(t02c, t02cd);
    addIntegrator(t10c, t10cd);
    addIntegrator(t11c, t11cd);
    addIntegrator(t12c, t12cd);
    addIntegrator(t20c, t20cd);
    addIntegrator(t21c, t21cd);
    addIntegrator(t22c, t22cd);
}

void INS::init() {
    if (initCount == 0) {
        gps_update_count = 0;
        startrack_update_count = 0;
        // At startup, the INS is "aligned" with truth: its initial
        // position, velocity, and attitude match the physics state.
        if (newton) {
            sxc = newton->SBII.x;
            syc = newton->SBII.y;
            szc = newton->SBII.z;
            vxc = newton->VBII.x;
            vyc = newton->VBII.y;
            vzc = newton->VBII.z;
            SBIIC = newton->SBII;
            VBIIC = newton->VBII;
        }
        if (kin) {
            osk::Mat T = kin->TBI;
            t00c = T[0][0]; t01c = T[0][1]; t02c = T[0][2];
            t10c = T[1][0]; t11c = T[1][1]; t12c = T[1][2];
            t20c = T[2][0]; t21c = T[2][1]; t22c = T[2][2];
            TBIC = T;
        }
    }
}

void INS::update() {
    if (mins == 0) {
        // ---- Ideal INS: outputs equal truth ----
        if (newton) {
            SBIIC = newton->SBII;
            VBIIC = newton->VBII;
            FSPCB = newton->FSPB;
            dbic  = SBIIC.mag();
        }
        if (euler) {
            WBICB = euler->WBIB;
            WBICI = euler->WBII;
        }
        if (kin) {
            TBIC = kin->TBI;
        }
        // Body rates from WBEB so they're Earth-relative (what an
        // accelerometer-coupled gyro reports for autopilot use).
        if (euler) {
            ppcx = euler->ppx;
            qqcx = euler->qqx;
            rrcx = euler->rrx;
        }
        // Hold integrator states pinned to truth (zero derivative).
        // This way the registered states don't drift when mode=0.
        if (newton) {
            sxc = newton->SBII.x;  syc = newton->SBII.y;  szc = newton->SBII.z;
            vxc = newton->VBII.x;  vyc = newton->VBII.y;  vzc = newton->VBII.z;
        }
        if (kin) {
            osk::Mat T = kin->TBI;
            t00c = T[0][0]; t01c = T[0][1]; t02c = T[0][2];
            t10c = T[1][0]; t11c = T[1][1]; t12c = T[1][2];
            t20c = T[2][0]; t21c = T[2][1]; t22c = T[2][2];
        }
        sxcd = sycd = szcd = 0.0;
        vxcd = vycd = vzcd = 0.0;
        t00cd = t01cd = t02cd = 0.0;
        t10cd = t11cd = t12cd = 0.0;
        t20cd = t21cd = t22cd = 0.0;

        ins_pos_err = 0.0;
        ins_vel_err = 0.0;
        ins_att_err = 0.0;
        dvbec = newton ? newton->dvbe : 0.0;
        return;
    }

    // ---- mins=1: simple-bias INS, self-integrating ----
    // Read truth measurements from physics blocks
    osk::Vec FSPB_true(0,0,0), WBIB_true(0,0,0), GRAVG(0,0,0);
    osk::Vec SBII_true(0,0,0), VBII_true(0,0,0);
    if (newton) {
        FSPB_true = newton->FSPB;
        SBII_true = newton->SBII;
        VBII_true = newton->VBII;
        // Pull gravity from Environment via Newton's stored copy.  Since
        // Newton applies gravity inertially with ABII = (1/m)*~TBI*FAPB
        // + GRAV_inertial, we need GRAV_inertial here.  Approximate as
        // inverse-square gravity along -SBII direction.
        double dbi = SBII_true.mag();
        if (dbi > 1.0e-3) {
            const double MU_E = 3.986004418e14;
            double g_mag = MU_E / (dbi * dbi);
            GRAVG = SBII_true * (-g_mag / dbi);
        }
    }
    if (euler) {
        WBIB_true = euler->WBIB;
    }

    // Apply biases to produce "measured" values
    osk::Vec FSPCB_meas = FSPB_true  + bias_accel;
    osk::Vec WBICB_meas = WBIB_true  + bias_gyro;

    // Pull integrator scalars into outputs
    TBIC  = get_mat(t00c, t01c, t02c, t10c, t11c, t12c, t20c, t21c, t22c);
    SBIIC = osk::Vec(sxc, syc, szc);
    VBIIC = osk::Vec(vxc, vyc, vzc);

    // ---- GPS update correction ----
    // If a GPS measurement is available this step, replace the INS
    // position/velocity state with the GPS measurement (loosely-coupled
    // full-state update).  This bounds the INS drift to GPS noise level
    // between updates.  We also reset the integrator scalars so RK4
    // continues from the corrected state, and clear the GPS flag.
    //
    // Important: GPS only produces measurements at stepstart=true (see
    // gps.cpp), which is the only place where RK4 will respect the
    // overwrite of integrator scalars (propagate(0) caches state_ to
    // save_ at this stage; subsequent stages use save_).
    if (gps && gps->gps_update_avail) {
        SBIIC = gps->SBII_meas;
        VBIIC = gps->VBII_meas;
        sxc = SBIIC.x;  syc = SBIIC.y;  szc = SBIIC.z;
        vxc = VBIIC.x;  vyc = VBIIC.y;  vzc = VBIIC.z;
        gps->gps_update_avail = 0;
        gps_update_count++;
    }

    // ---- Star-tracker attitude correction ----
    // URIC is a small-angle tilt-error vector.  Apply as
    //   TBIC := (I + [URIC]_x) * TBIC
    // followed by one Bar-Itzhack first-order orthonormalization step.
    // Then refresh the t**c integrator scalars from the corrected DCM.
    // Same RK4-stepstart-gating rationale as GPS: Startrack produces
    // URIC only at stepstart=true so propagate(0) sees the corrected
    // state.
    if (startrack && startrack->star_update_avail) {
        osk::Mat I3(1, 0, 0, 0, 1, 0, 0, 0, 1);
        osk::Mat R_corr = I3 + skew(startrack->URIC);
        osk::Mat TBIC_new = R_corr * TBIC;
        // Bar-Itzhack one-step orthonormalization:
        //   T <- T - 0.5 * (T*T^T - I) * T
        osk::Mat TT = TBIC_new * TBIC_new.transpose();
        osk::Mat E  = TT - I3;
        TBIC = TBIC_new - (E * TBIC_new) * 0.5;
        // Refresh integrator scalars
        t00c = TBIC[0][0]; t01c = TBIC[0][1]; t02c = TBIC[0][2];
        t10c = TBIC[1][0]; t11c = TBIC[1][1]; t12c = TBIC[1][2];
        t20c = TBIC[2][0]; t21c = TBIC[2][1]; t22c = TBIC[2][2];
        startrack->star_update_avail = 0;
        startrack_update_count++;
    }

    // ---- Compute derivatives from MEASURED values ----
    // Inertial position: dS/dt = V
    sxcd = vxc;  sycd = vyc;  szcd = vzc;

    // Inertial velocity: dV/dt = ~TBIC * FSPCB + GRAVG(based on SBIIC)
    osk::Vec accel_I = TBIC.transpose() * FSPCB_meas;
    // Use measured gravity at SBIIC (not at truth position) -- INS
    // doesn't know truth.
    osk::Vec GRAVG_c(0,0,0);
    double dbic_loc = SBIIC.mag();
    if (dbic_loc > 1.0e-3) {
        const double MU_E = 3.986004418e14;
        double g_mag = MU_E / (dbic_loc * dbic_loc);
        GRAVG_c = SBIIC * (-g_mag / dbic_loc);
    }
    accel_I = accel_I + GRAVG_c;
    vxcd = accel_I.x;  vycd = accel_I.y;  vzcd = accel_I.z;

    // Attitude: dTBIC/dt = -[WBICB_meas]_x * TBIC
    osk::Mat dT = skew(WBICB_meas) * TBIC * (-1.0);
    t00cd = dT[0][0]; t01cd = dT[0][1]; t02cd = dT[0][2];
    t10cd = dT[1][0]; t11cd = dT[1][1]; t12cd = dT[1][2];
    t20cd = dT[2][0]; t21cd = dT[2][1]; t22cd = dT[2][2];

    // ---- Compute outputs from INS state and measurements ----
    FSPCB = FSPCB_meas;
    WBICB = WBICB_meas;
    WBICI = TBIC.transpose() * WBICB;
    dbic  = SBIIC.mag();

    // Body rates in deg/s.  In the non-rotating Earth model used here,
    // WBEB = WBIB (Earth-relative equals inertial body rates).
    ppcx = WBICB.x * DEG;
    qqcx = WBICB.y * DEG;
    rrcx = WBICB.z * DEG;

    dvbec = newton ? newton->dvbe : 0.0;

    // Diagnostics: error magnitudes vs truth
    ins_pos_err = (SBIIC - SBII_true).mag();
    ins_vel_err = (VBIIC - VBII_true).mag();
    // Attitude error: small-angle magnitude of (TBI_true * TBIC^T - I)
    if (kin) {
        osk::Mat I3(1, 0, 0, 0, 1, 0, 0, 0, 1);
        osk::Mat E_att = kin->TBI * TBIC.transpose() - I3;
        // The skew-symmetric part encodes the small-angle vector.
        osk::Vec u(E_att[2][1], E_att[0][2], E_att[1][0]);
        ins_att_err = u.mag();
    }
}

void INS::rpt() {
    if (osk::State::sample(1.0)) {
        std::printf("INS  t=%7.3f  mins=%d  pos_err=%8.3f m  vel_err=%8.4f m/s  "
                    "qqcx=%+.3f deg/s\n",
                    osk::State::t, mins, ins_pos_err, ins_vel_err, qqcx);
    }
}

} // namespace rocket6dof
