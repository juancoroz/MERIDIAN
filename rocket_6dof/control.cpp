//  control.cpp  --  Two-channel acceleration autopilot
//
//  Cross-reference to Zipfel control.cpp (Section 10.3.x):
//
//      Zipfel name      OSK member        Meaning
//      -----------      ----------        -------
//      maut             maut              mode = |mauty|mautp|
//      waclp,zaclp,     waclp,zaclp,      pitch closed-loop pole params
//        paclp            paclp
//      wacly,zacly,     wacly,zacly,      yaw closed-loop pole params
//        pacly            pacly
//      delimx,drlimx    delimx,drlimx     deflection limiters
//      gnmax,gymax      gnmax,gymax       accel command limiters
//      ancomx,alcomx    ancomx,alcomx     accel commands [g's]
//      zz, zzd          zz, zzd           pitch integrator state
//      yy, yyd          yy, yyd           yaw integrator state
//      delecx, delrcx   delecx, delrcx    pitch / yaw deflection commands
//      GAINFP, GAINFY   gainfp1..3,       diagnostic feedback gains
//                       gainfy1..3
//
//  Equations (Zipfel control.cpp lines 188-199 for pitch; 269-280 for yaw):
//
//    Pole-placement gains:
//      gainfp3 = waclp^2 * paclp / (dla * dmde)
//      gainfp2 = (2*zaclp*waclp + paclp + dmq - dla/V) / dmde
//      gainfp1 = (waclp^2 + 2*zaclp*waclp*paclp + dma + dmq*dla/V
//                 - gainfp2*dmde*dla/V) / (dla * dmde)
//
//    Pitch loop:
//      zzd     = g0 * ancomx + fspb_z         (accel-error integrator)
//      dqc_rad = -gainfp1 * (-fspb_z)
//                 - gainfp2 * qqcx_rad
//                 + gainfp3 * zz
//      delecx  = dqc_rad * (180/pi)
//
//  Note: Zipfel's gain expressions implicitly assume specific signs.
//  Pitch: positive ancomx -> rocket nose up -> positive delecx command
//  produces a NEGATIVE pitching moment (dmde is negative) which tilts
//  thrust to push the nose up.

#include "control.h"
#include "environment.h"
#include "newton.h"
#include "euler.h"
#include "kinematics.h"
#include "aerodynamics.h"
#include "propulsion.h"
#include "ins.h"
#include <cmath>
#include <cstdio>

namespace rocket6dof {

namespace {
constexpr double G0   = 9.80675445;
constexpr double DEG  = 180.0 / osk::PI;
constexpr double RAD  = osk::PI / 180.0;

inline double sign(double x) { return (x > 0.0) - (x < 0.0); }

// Wrap an angle difference into (-180, +180] degrees.
inline double wrap180(double d) {
    while (d >  180.0) d -= 360.0;
    while (d <= -180.0) d += 360.0;
    return d;
}
} // anon

Control::Control()
    : env(nullptr), newton(nullptr), euler(nullptr),
      kin(nullptr),
      aero(nullptr), prop(nullptr), ins(nullptr),
      maut(0),
      vac_rate_damp(0.5),
      vac_max_gain(1.0),
      tau_att(2.0),
      q_max(20.0),
      theta_com_inertial(0.0),
      psi_com_inertial(0.0),
      waclp(0.5), zaclp(0.7), paclp(0.7), factwaclp(0.0),
      wacly(0.5), zacly(0.7), pacly(0.7), factwacly(0.0),
      delimx(30.0), drlimx(30.0),
      gnmax(30.0), gymax(30.0),
      ancomx(0.0), alcomx(0.0),
      zz(0.0), zzd(0.0),
      yy(0.0), yyd(0.0),
      delecx(0.0), delrcx(0.0),
      ancomx_actual(0.0), alcomx_actual(0.0),
      gainfp1(0), gainfp2(0), gainfp3(0),
      gainfy1(0), gainfy2(0), gainfy3(0)
{
    // Register the two integrator states.
    addIntegrator(zz, zzd);
    addIntegrator(yy, yyd);
}

void Control::init() {
    if (initCount == 0) {
        zz = 0.0;  zzd = 0.0;
        yy = 0.0;  yyd = 0.0;
        delecx = 0.0;
        delrcx = 0.0;
        // Washout filter states for accel-feedback load relief.
        // Initialize to 0; first ~5*tau worth of samples are the
        // filter's startup transient, during which it adapts to
        // the actual steady-state anx.  At dt=0.01s, tau=1s gives
        // about 500 sim steps of startup -- well before max-Q.
        anx_lp_state = 0.0;
        ayx_lp_state = 0.0;
    }
}

void Control::update() {
    // ===== VACUUM-MODE AUTOPILOT (maut = 99) ============================
    // Attitude-hold rate damper for the exo-atmospheric phase.  The
    // standard maut=53 pole-placement autopilot fails in vacuum
    // (dla*dmde -> 0 when q -> 0), so this mode replaces it with a
    // simple gimbal-based body-rate damper that does not depend on aero.
    //
    // Closed-loop dynamics for a single axis (small-angle, no aero):
    //
    //   I * q_dot = -thrust * arm * delta_e_rad                   (moment EOM)
    //   delta_e_rad = (I / (thrust * arm * tau)) * q_rad          (control law)
    //
    //   => q_dot = -q / tau    (exponential decay with time constant tau)
    //
    //  vac_rate_damp parameterizes tau (in seconds).  Typical value 0.5 s.
    //  Sign convention: positive delecx -> FPB.z negative -> FMPB.y
    //  negative -> qqd negative.  Since positive qqx is nose-DOWN
    //  rotation (Zipfel body frame), positive delecx drives qqx down
    //  (nose-up).  For damping a positive qqx, we want positive delecx
    //  -> the control law sign is naturally correct:
    //
    //      delecx = +(I/(thrust*arm*tau)) * q_rad
    //
    //  Yaw channel: positive delrcx -> FPB.y positive -> FMPB.z negative
    //  -> rrd negative.  Since positive rrx is nose-LEFT rotation
    //  (around +z body axis, RHR), positive delrcx drives rrx down
    //  (nose-right).  To damp positive rrx, we want positive delrcx.
    //
    //  ancomx/alcomx are ignored in this mode -- this is purely an
    //  attitude-hold rate damper for the coast/post-burnout phase.
    if (maut == 99) {
        delecx = 0.0;
        delrcx = 0.0;
        zzd = 0.0;
        yyd = 0.0;

        // Read body rates [rad/s] from INS when available, else Euler truth
        double q_rad = 0.0, r_rad = 0.0;
        if (ins) {
            q_rad = ins->qqcx * RAD;
            r_rad = ins->rrcx * RAD;
        } else if (euler) {
            q_rad = euler->qqx * RAD;
            r_rad = euler->rrx * RAD;
        }

        // Pull TVC/Propulsion data needed for the gain
        double thr = (prop ? prop->thrust : 0.0);
        double Iy  = (prop ? prop->IBBB[1][1] : 1.0);
        double Iz  = (prop ? prop->IBBB[2][2] : 1.0);

        // The gimbal-arm distance.  Control doesn't carry a TVC pointer,
        // so it uses a representative default (gimbal aft of CG by ~2 m),
        // or zeros out when there's no thrust.
        double arm = 2.0;
        if (prop) {
            // Compute gimbal moment arm from the thrust attachment
            // location (thrust_loc) minus the current CG location.
            // thrust_loc defaults to 5.0 m (small-vehicle assumption);
            // override via control.thrust_loc in config for larger
            // vehicles where thrust is attached further aft.  Should
            // match tvc.parm when TVC is used.
            arm = thrust_loc - prop->xcg;
            if (arm < 0.5) arm = 0.5;
        }

        // Time-constant parameterization: vac_rate_damp is tau (seconds)
        double tau = (vac_rate_damp > 1.0e-3) ? vac_rate_damp : 0.5;

        if (thr > 1.0) {
            double K_inv_p = Iy / (thr * arm * tau);   // [rad of defl per rad/s]
            double K_inv_y = Iz / (thr * arm * tau);
            // Cap gain to prevent runaway when thrust is small
            double max_gain = (vac_max_gain > 0.0) ? vac_max_gain : 1.0;
            if (K_inv_p > max_gain) K_inv_p = max_gain;
            if (K_inv_y > max_gain) K_inv_y = max_gain;
            delecx = K_inv_p * q_rad * DEG;
            delrcx = K_inv_y * r_rad * DEG;
        }

        // Saturate
        if (std::fabs(delecx) > delimx) delecx = delimx * sign(delecx);
        if (std::fabs(delrcx) > drlimx) delrcx = drlimx * sign(delrcx);
        return;
    }

    // ===== CASCADED ATTITUDE-TRACKING AUTOPILOT (maut = 50 or 60) =======
    // Two-loop cascaded control:
    //
    //   OUTER LOOP (attitude tracking):
    //     theta_err = theta_com - theta_act        [deg]
    //     q_target  = clamp(theta_err / tau_att, +/- q_max)   [deg/s]
    //
    //   INNER LOOP (rate tracking, same physics as maut=99):
    //     q_err      = q_target - q_act            [deg/s]
    //     delecx_rad = -(I_y / (T*arm*tau_rate)) * q_err_rad
    //
    // Sign discipline (verified by simulation):
    //   theta positive = nose ABOVE horizontal (geodetic body pitch)
    //   q_act = qqx = WBIB.y * DEG, positive = nose-UP rotation
    //     (rotation around +body-y, RHR takes +z down -> +x forward)
    //   positive delecx -> FPB.z = -seta*thr negative -> FMPB.y =
    //     arm*FPB.z = +*neg = NEGATIVE moment -> WBIBD.y negative ->
    //     WBIB.y becomes negative = nose-DOWN.
    //
    //   So:
    //   * To increase theta (raise nose), we want positive q_target
    //   * For positive q_target with q_act=0, q_err > 0
    //   * positive q_err -> negative delecx (to produce nose-up moment)
    //   * Hence: delecx = -K_inv * q_err
    //
    //   Yaw analogue: positive delrcx -> FPB.y positive -> FMPB.z =
    //     -arm*FPB.y negative -> WBIB.z becomes negative = nose-LEFT.
    //     But psibdx positive = nose-RIGHT (Zipfel/aerospace).  So
    //     positive psi_err (need nose-right) needs r_target positive,
    //     but +WBIB.z = nose-LEFT.  So r_target = -psi_err/tau_att.
    //     For positive r_err (need MORE positive r = MORE nose-left),
    //     we need positive delrcx (more negative WBIB.z = nose-left).
    //
    //   maut = 60: theta_com_inertial, psi_com_inertial set externally
    //   maut = 50: converted from ancomx/alcomx using aero
    if (maut == 50 || maut == 60) {
        delecx = 0.0;
        delrcx = 0.0;
        zzd = 0.0;
        yyd = 0.0;

        // Read actual body attitude from Kinematics (truth).
        double theta_act = (kin ? kin->thtbdx : 0.0);     // [deg]
        double psi_act   = (kin ? kin->psibdx : 0.0);     // [deg]

        // Read body rates [rad/s] from INS when available
        double q_rad = 0.0, r_rad = 0.0;
        if (ins) {
            q_rad = ins->qqcx * RAD;
            r_rad = ins->rrcx * RAD;
        } else if (euler) {
            q_rad = euler->qqx * RAD;
            r_rad = euler->rrx * RAD;
        }

        // Determine attitude targets
        double theta_com_eff = theta_com_inertial;
        double psi_com_eff   = psi_com_inertial;

        if (maut == 50) {
            // Acceleration-driven mode: ancomx -> alpha_command -> theta_com
            // Limit ancomx/alcomx to max-g
            double anc = ancomx;
            if (anc >  gnmax) anc =  gnmax;
            if (anc < -gnmax) anc = -gnmax;
            ancomx_actual = anc;
            double alc = alcomx;
            if (alc >  gymax) alc =  gymax;
            if (alc < -gymax) alc = -gymax;
            alcomx_actual = alc;

            // Aero gains
            double dla_eff = (aero ? aero->dla : 0.0);
            double dyb_eff = (aero ? aero->dyb : 0.0);
            // Velocity-vector attitude (FPA and heading)
            double thtvdx = (newton ? newton->thtvdx : 0.0);  // [deg]
            double psivdx = (newton ? newton->psivdx : 0.0);  // [deg]

            // Compute alpha command (positive alpha = body nose above velocity)
            // Only invert when dla is large enough (i.e. q is non-trivial)
            double alpha_cmd_deg = 0.0;
            double beta_cmd_deg  = 0.0;
            if (dla_eff > 0.1) {  // [m/s^2 per rad of alpha]
                alpha_cmd_deg = (anc * G0 / dla_eff) * DEG;
            }
            if (std::fabs(dyb_eff) > 0.1) {
                beta_cmd_deg  = (alc * G0 / dyb_eff) * DEG;
            }
            // Limit to reasonable alpha range (avoid stall)
            const double ALPHA_MAX = 10.0;  // [deg]
            if (alpha_cmd_deg >  ALPHA_MAX) alpha_cmd_deg =  ALPHA_MAX;
            if (alpha_cmd_deg < -ALPHA_MAX) alpha_cmd_deg = -ALPHA_MAX;
            if (beta_cmd_deg  >  ALPHA_MAX) beta_cmd_deg  =  ALPHA_MAX;
            if (beta_cmd_deg  < -ALPHA_MAX) beta_cmd_deg  = -ALPHA_MAX;

            theta_com_eff = thtvdx + alpha_cmd_deg;
            psi_com_eff   = psivdx + beta_cmd_deg;
        }

        // ---- Load relief (maut=60 only, opt-in via q_threshold > 0) ----
        // Blend theta_com_eff toward the velocity-vector pitch angle
        // (FPA) when dynamic pressure is high.  At full blend the
        // autopilot is commanding alpha=0; at zero blend it follows
        // the open-loop schedule unchanged.  Smooth sigmoid in
        // q_dyn keeps the transition continuous so the autopilot
        // doesn't see a step change in command at any specific Mach.
        //
        // This is a simplified version of real launcher load relief
        // (which typically uses accelerometer feedback to detect
        // transverse loads).  The velocity-frame-tracking variant
        // here is sufficient to demonstrate the principle and
        // reduce peak q*alpha in MC studies.
        if (maut == 60 && load_relief_q_threshold > 0.0 && env) {
            double q = std::fabs(env->pdynmc);
            // Smooth sigmoid: w(q) = 1 / (1 + exp(-(q - q0) / (q_width/2)))
            // q_width is the half-width, so the sigmoid scale uses
            // q_width directly (gives 0.5 at q=q0, ~0.88 at q0+q_width,
            // ~0.12 at q0-q_width).
            double scale = (load_relief_q_width > 1.0) ? load_relief_q_width : 1.0;
            double x = (q - load_relief_q_threshold) / scale;
            // Clip exp argument to avoid overflow on extreme inputs
            if      (x >  30.0) x =  30.0;
            else if (x < -30.0) x = -30.0;
            double w = 1.0 / (1.0 + std::exp(-x));
            // Blend toward FPA (= thtvdx in Newton's geographic frame)
            double thtvdx = (newton ? newton->thtvdx : theta_com_eff);
            double psivdx = (newton ? newton->psivdx : psi_com_eff);
            theta_com_eff = (1.0 - w) * theta_com_eff + w * thtvdx;
            psi_com_eff   = (1.0 - w) * psi_com_eff   + w * psivdx;
        }

        // ---- OUTER LOOP: attitude -> rate command ----
        double tau_a = (tau_att > 1.0e-3) ? tau_att : 2.0;

        // ---- Load-relief tau detuning ----
        // When accel feedback is engaged (gain > 0 and q above threshold),
        // SLOW the outer attitude loop by load_relief_tau_factor so the
        // feedback term has room to act without the autopilot fighting
        // it.  Uses the SAME q_dyn sigmoid as the feedback term, so the
        // two effects engage/disengage together.  Default factor=1.0
        // means no detuning (backward compatible).
        if (maut == 60 && accel_relief_gain > 0.0
            && load_relief_tau_factor > 1.0
            && env && env->pdynmc > 0.0)
        {
            double q_lr = env->pdynmc;
            double scale = (accel_relief_q_width > 1.0) ? accel_relief_q_width : 1.0;
            double x = (q_lr - accel_relief_q_threshold) / scale;
            if      (x >  30.0) x =  30.0;
            else if (x < -30.0) x = -30.0;
            double w_lr = 1.0 / (1.0 + std::exp(-x));
            // Blend tau between nominal and slowed value
            tau_a *= (1.0 + w_lr * (load_relief_tau_factor - 1.0));
        }

        double theta_err = wrap180(theta_com_eff - theta_act);  // [deg]
        double psi_err   = wrap180(psi_com_eff   - psi_act);    // [deg]

        // q_target sign: positive theta_err means we need more nose-up.
        // Positive WBIB.y = nose-up, so q_target should be POSITIVE for
        // positive theta_err.
        double q_target_dps = +theta_err / tau_a;
        // Yaw: positive psi_err means we need more nose-right (psi
        // increasing).  By RHR, +z rotation takes +x to +y (forward
        // to right) = nose-right.  So positive WBIB.z = nose-right,
        // and r_target = +psi_err/tau.
        double r_target_dps = +psi_err / tau_a;
        // Saturate to q_max
        double qm = (q_max > 0.0) ? q_max : 20.0;
        if (q_target_dps >  qm) q_target_dps =  qm;
        if (q_target_dps < -qm) q_target_dps = -qm;
        if (r_target_dps >  qm) r_target_dps =  qm;
        if (r_target_dps < -qm) r_target_dps = -qm;

        // ---- INNER LOOP: rate error -> deflection ----
        // Plant: Iy*qqd = M_tvc + M_aero
        //   M_tvc  = -arm*thrust*delecx_rad (for small delecx)
        //   M_aero = q*S*L*clm  (clm = Cm from Aero tables, with damping)
        //
        // Desired behavior: qqd = -(qqx - q_target)/tau_r, i.e., closed-
        // loop time constant tau_r.  Solving for delecx_rad:
        //
        //   delecx_rad = [M_aero + Iy*q_err/tau_r] / (arm*thrust)
        //
        // This includes feed-forward aero-moment compensation.  Without
        // the M_aero term, the inner loop loses authority at high q.
        double thr = (prop ? prop->thrust : 0.0);
        double Iy  = (prop ? prop->IBBB[1][1] : 1.0);
        double Iz  = (prop ? prop->IBBB[2][2] : 1.0);
        double arm = 2.0;
        if (prop) {
            arm = thrust_loc - prop->xcg;
            if (arm < 0.5) arm = 0.5;
        }
        double tau_r = (vac_rate_damp > 1.0e-3) ? vac_rate_damp : 0.5;

        if (thr > 1.0) {
            // Rate error in rad/s.  Sign chosen so that positive q_err
            // means q_act is too positive (need to reduce), wanting
            // positive delecx (which creates nose-DOWN moment).
            double q_err_rad = q_rad - (q_target_dps * RAD);
            double r_err_rad = r_rad - (r_target_dps * RAD);

            // Feed-forward aero moment compensation.  Read aero's
            // current clm/cln coefficients to anticipate the aero
            // moment we need to cancel.  Only when aero is wired and
            // in atmospheric flight (q > 1 kPa).
            double M_aero_y = 0.0;
            double M_aero_z = 0.0;
            if (aero && env && env->pdynmc > 1000.0) {
                double qSL = env->pdynmc * aero->refa * aero->refd;
                M_aero_y = qSL * aero->clm;
                M_aero_z = qSL * aero->cln;
            }

            double delecx_rad = (M_aero_y + Iy*q_err_rad/tau_r) / (arm*thr);
            double delrcx_rad = (M_aero_z + Iz*r_err_rad/tau_r) / (arm*thr);

            // Saturate at deflection limit
            double max_rad = delimx * RAD;
            if (delecx_rad >  max_rad) delecx_rad =  max_rad;
            if (delecx_rad < -max_rad) delecx_rad = -max_rad;
            max_rad = drlimx * RAD;
            if (delrcx_rad >  max_rad) delrcx_rad =  max_rad;
            if (delrcx_rad < -max_rad) delrcx_rad = -max_rad;

            delecx = delecx_rad * DEG;
            delrcx = delrcx_rad * DEG;
        }

        // ---- Accelerometer-feedback load relief (maut=60, opt-in) ----
        // Production-style load relief with washout filter:
        // subtract a gimbal-deflection term proportional to the
        // HIGH-PASS-FILTERED body lateral acceleration.  The
        // high-pass filter (washout) is the key element that
        // distinguishes this from the broken naive version:
        // it isolates transient lateral loads (gusts, q*alpha
        // excursions, wind shear) from the steady-state lateral
        // acceleration the autopilot needs to shape the gravity
        // turn.  Without it, the feedback fights the autopilot.
        //
        // First-order discrete IIR low-pass: anx_lp tracks the
        // slow DC component of anx with time constant tau.
        // The high-pass output anx_hp = anx - anx_lp captures
        // only the rapid deviations.
        //
        // Sign convention: positive anx = body accelerating
        // "up" (nose pushed up); to counter, decrease delecx
        // (less nose-up gimbal command).  Feedback subtracts.
        //
        // Gated by q_dyn through a sigmoid: feedback only active
        // during high-q phase (no point applying when q ~ 0).
        if (maut == 60 && accel_relief_gain > 0.0 && env && newton) {
            // Sample time for the IIR update.  osk::Sim drives at
            // dt=0.01s nominally; reading State::dt keeps the filter
            // physically correct if dt changes.
            double dt = osk::State::dt;
            double tau = (accel_relief_tau > 1.0e-6) ? accel_relief_tau : 1.0;
            double alpha = dt / (tau + dt);  // IIR coefficient
            // Update low-pass states
            double anx_now = newton->anx;
            double ayx_now = newton->ayx;
            anx_lp_state += alpha * (anx_now - anx_lp_state);
            ayx_lp_state += alpha * (ayx_now - ayx_lp_state);
            // High-pass = raw - low-pass
            double anx_hp = anx_now - anx_lp_state;
            double ayx_hp = ayx_now - ayx_lp_state;
            // q_dyn sigmoid gate
            double q = std::fabs(env->pdynmc);
            double scale = (accel_relief_q_width > 1.0) ? accel_relief_q_width : 1.0;
            double x = (q - accel_relief_q_threshold) / scale;
            if      (x >  30.0) x =  30.0;
            else if (x < -30.0) x = -30.0;
            double w = 1.0 / (1.0 + std::exp(-x));
            // Apply feedback using HIGH-PASS signal
            delecx -= w * accel_relief_gain * anx_hp;
            delrcx -= w * accel_relief_gain * ayx_hp;
        }

        // Saturate
        if (std::fabs(delecx) > delimx) delecx = delimx * sign(delecx);
        if (std::fabs(delrcx) > drlimx) delrcx = drlimx * sign(delrcx);
        return;
    }

    // Decode mode: maut = mauty*10 + mautp
    int mauty = maut / 10;
    int mautp = maut % 10;

    // ---- Read upstream sensors ----
    // When INS is wired, sensor signals come from INS (which models bias,
    // drift, etc.).  Otherwise fall back to direct truth from Newton/Euler.
    double pdynmc = (env ? env->pdynmc : 0.0);
    double dvbe   = (ins    ? ins->dvbec
                    : newton ? newton->dvbe
                    : 0.0);
    double V      = (dvbe > 1.0 ? dvbe : 1.0);

    // Body specific force (m/s^2)
    double fspb_y, fspb_z;
    if (ins) {
        fspb_y = ins->FSPCB.y;
        fspb_z = ins->FSPCB.z;
    } else if (newton) {
        fspb_y = newton->FSPB.y;
        fspb_z = newton->FSPB.z;
    } else {
        fspb_y = fspb_z = 0.0;
    }

    // Body rates (deg/s)
    double qqcx_deg, rrcx_deg;
    if (ins) {
        qqcx_deg = ins->qqcx;
        rrcx_deg = ins->rrcx;
    } else if (euler) {
        qqcx_deg = euler->qqx;
        rrcx_deg = euler->rrx;
    } else {
        qqcx_deg = rrcx_deg = 0.0;
    }

    // Aero dimensional derivatives
    double dla  = (aero ? aero->dla  : 0.0);
    double dma  = (aero ? aero->dma  : 0.0);
    double dmq  = (aero ? aero->dmq  : 0.0);
    double dmde = (aero ? aero->dmde : 0.0);
    double dyb  = (aero ? aero->dyb  : 0.0);
    double dnb  = (aero ? aero->dnb  : 0.0);
    double dnr  = (aero ? aero->dnr  : 0.0);
    double dndr = (aero ? aero->dndr : 0.0);

    // Thrust state -- control fires only while thrusting (TVC needs thrust)
    bool thrust_on = (prop && prop->thrust > 0.0);

    // ---- Gain-schedule the closed-loop poles with dynamic pressure ----
    // Zipfel's exact form (control.cpp lines 184-185 / 265-266).  Below
    // q=20 kPa the formula gives negative waclp; clamp to a small positive
    // value to keep gains finite.
    double waclp_eff = (0.1 + 0.5e-5 * (pdynmc - 20.0e3)) * (1.0 + factwaclp);
    double paclp_eff = (0.7 + 1.0e-5 * (pdynmc - 20.0e3)) * (1.0 + factwaclp);
    if (waclp_eff < 0.05) waclp_eff = 0.05;
    if (paclp_eff < 0.1)  paclp_eff = 0.1;
    waclp = waclp_eff;
    paclp = paclp_eff;
    double wacly_eff = (0.1 + 0.5e-5 * (pdynmc - 20.0e3)) * (1.0 + factwacly);
    double pacly_eff = (0.7 + 1.0e-5 * (pdynmc - 20.0e3)) * (1.0 + factwacly);
    if (wacly_eff < 0.05) wacly_eff = 0.05;
    if (pacly_eff < 0.1)  pacly_eff = 0.1;
    wacly = wacly_eff;
    pacly = pacly_eff;

    // ---- Pitch channel (mautp == 3) ----
    delecx = 0.0;
    if (mautp == 3) {
        // Limit accel command to max-g capability
        double anc = ancomx;
        if (anc >  gnmax) anc =  gnmax;
        if (anc < -gnmax) anc = -gnmax;
        ancomx_actual = anc;

        // Low-q guard: this autopilot's pole-placement formula was
        // derived assuming nontrivial aerodynamic feedback.  When q is
        // very low, the dimensional derivatives dla, dmde are tiny and
        // the gain formula amplifies noise.  Hold delecx at zero below
        // 5 kPa (well below the Zipfel gain-schedule break point of
        // 20 kPa where waclp/paclp clamp out of the formula).  Above
        // 5 kPa, run the pole-placement formula normally.
        if (thrust_on && pdynmc > 5000.0) {
            // Avoid divide-by-zero when aero derivatives are still zero
            // (early in flight, before q has built up)
            double denom = dla * dmde;
            if (std::fabs(denom) > 1.0e-10) {
                // Zipfel's pole-placement formula was derived for an aero
                // elevator (dmde > 0).  For TVC, dmde < 0 in the raw
                // Aero derivative because the moment arm is aft of CG.
                // Empirically with s = +1, the closed loop is stable
                // but the response is slow (gains are dominated by
                // gravity-induced fspb_z disturbance rather than
                // ancomx command).  The MUCH preferred path for our
                // vehicle is the cascaded attitude-tracking autopilot
                // (maut = 50 or 60), which has clean inner/outer loop
                // separation and works in both atmosphere and vacuum.
                // This maut=53 mode is retained for legacy / reference.
                double s = 1.0;

                gainfp3 = s * waclp*waclp*paclp / denom;
                gainfp2 = s * (2.0*zaclp*waclp + paclp + dmq - dla/V) / dmde;
                gainfp1 = s * (waclp*waclp + 2.0*zaclp*waclp*paclp + dma
                            + dmq*dla/V - gainfp2*s*dmde*dla/V) / denom;

                // Integrator derivative
                zzd = G0 * anc + fspb_z;

                // Pitch deflection command, in radians, then convert to deg
                double dqc = -gainfp1*(-fspb_z)
                             - gainfp2*qqcx_deg*RAD
                             + gainfp3*zz;
                delecx = dqc * DEG;
            } else {
                zzd = G0 * anc + fspb_z;
            }
        } else {
            zzd = 0.0;
        }
    } else {
        zzd = 0.0;
    }

    // ---- Yaw channel (mauty == 5) ----
    delrcx = 0.0;
    if (mauty == 5) {
        double alc = alcomx;
        if (alc >  gymax) alc =  gymax;
        if (alc < -gymax) alc = -gymax;
        alcomx_actual = alc;

        // Same low-q guard as pitch channel.
        if (thrust_on && pdynmc > 5000.0) {
            double denom = dyb * dndr;
            if (std::fabs(denom) > 1.0e-10) {
                // Same sign-convention reasoning as pitch channel: use
                // unflipped formula (s=+1) for stable feedback.
                double s = 1.0;

                gainfy3 = -s * wacly*wacly*pacly / denom;
                gainfy2 = s * (2.0*zacly*wacly + pacly + dnr + dyb/V) / dndr;
                gainfy1 = -s * (-wacly*wacly - 2.0*zacly*wacly*pacly + dnb
                            + dnr*dyb/V - gainfy2*s*dndr*dnb/V) / denom;

                yyd = G0 * alc - fspb_y;

                double drc = -gainfy1*fspb_y
                             - gainfy2*rrcx_deg*RAD
                             + gainfy3*yy;
                delrcx = drc * DEG;
            } else {
                yyd = G0 * alc - fspb_y;
            }
        } else {
            yyd = 0.0;
        }
    } else {
        yyd = 0.0;
    }

    // ---- Saturate deflections ----
    if (std::fabs(delecx) > delimx) delecx = delimx * sign(delecx);
    if (std::fabs(delrcx) > drlimx) delrcx = drlimx * sign(delrcx);
}

void Control::rpt() {
    if (osk::State::sample(1.0)) {
        std::printf("Ctrl t=%7.3f  ancomx=%7.3f  delecx=%+7.3f deg  "
                    "alcomx=%7.3f  delrcx=%+7.3f deg\n",
                    osk::State::t,
                    ancomx, delecx, alcomx, delrcx);
    }
}

} // namespace rocket6dof
