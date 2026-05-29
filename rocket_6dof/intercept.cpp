//  intercept.cpp  --  Trajectory termination logic
//
//  Cross-reference to Zipfel intercept.cpp (Section 10.3.x):
//
//      Zipfel name      OSK member       Meaning
//      -----------      ----------       -------
//      write (one-shot) fired_           prevents re-printing on
//                                        subsequent updates
//      (none -- alt<=0  trcond           condition code:
//       is the only        0 = running
//       condition)         1 = ground impact
//                          2 = apogee
//                          3 = max altitude
//                          4 = aero decay (M < trmach && q < trdynm)
//                          5 = max time
//
//  Sets Sim::stop = -1 on the first firing condition; the kernel then
//  terminates cleanly at the next step boundary.

#include "intercept.h"
#include "environment.h"
#include "newton.h"
#include "euler.h"
#include <cstdio>
#include <cmath>

namespace rocket6dof {

Intercept::Intercept()
    : env(nullptr), newton(nullptr), eul(nullptr),
      check_ground_impact(true),
      check_apogee(false),
      check_alt_max(false),
      check_aero_decay(false),
      check_time_max(false),
      check_tumble(false),
      alt_max(1.0e9),
      trmach(0.0),
      trdynm(0.0),
      t_max(1.0e9),
      rate_max(3.49),            // ~200 deg/s; vehicles tumbling this fast
                                 // are structurally lost
      t_no_check_until(0.5),     // small delay so launch IC doesn't trigger
      trcond(0),
      t_terminate(0.0),
      fired_(false),
      vz_prev_(0.0)
{
    // No integrator states.
}

void Intercept::init() {
    if (initCount == 0) {
        fired_      = false;
        trcond      = 0;
        t_terminate = 0.0;
        vz_prev_    = 0.0;
    }
    // On stage transitions (initCount > 0), keep fired_ as-is so a single
    // Intercept persists across stages.
}

void Intercept::update() {
    // ---- One-shot: don't keep firing after we've already triggered ----
    if (fired_) return;

    // ---- Guard: don't check until we've moved past launch ----
    if (osk::State::t < t_no_check_until) return;

    if (!newton || !env) return;

    // ---- Read state ----
    double alt    = newton->alt;
    double vmach  = env->vmach;
    double pdyn   = env->pdynmc;
    double dvbe   = newton->dvbe;
    double psivdx = newton->psivdx;
    double thtvdx = newton->thtvdx;
    double t      = osk::State::t;

    // Vertical velocity in NED frame: +up = -VBED.z (since z is Down).
    double vz_up  = -newton->VBED.z;

    int    code   = 0;
    const char* label = "";

    // ---- Check conditions in priority order ----
    if (check_ground_impact && alt <= 0.0) {
        code  = 1;  label = "Ground impact";
    }
    else if (check_apogee && osk::State::stepstart
             && vz_prev_ > 0.0 && vz_up <= 0.0) {
        code  = 2;  label = "Apogee reached";
    }
    else if (check_alt_max && alt > alt_max) {
        code  = 3;  label = "Maximum altitude exceeded";
    }
    else if (check_aero_decay && vmach < trmach && pdyn < trdynm) {
        code  = 4;  label = "Aero decay (low Mach + low q)";
    }
    else if (check_tumble && eul) {
        // |WBEB| = magnitude of body angular velocity vector (rad/s).
        // We use std::sqrt of x^2+y^2+z^2 directly rather than calling
        // a method, to be explicit and not depend on osk::Vec's API.
        double wx = eul->WBEB[0];
        double wy = eul->WBEB[1];
        double wz = eul->WBEB[2];
        double rate_mag = std::sqrt(wx*wx + wy*wy + wz*wz);
        if (rate_mag > rate_max) {
            code  = 6;  label = "Tumble detected (|body rate| > rate_max)";
        }
    }

    if (code == 0 && check_time_max && t > t_max) {
        code  = 5;  label = "Maximum flight time exceeded";
    }

    if (code != 0) {
        fired_      = true;
        trcond      = code;
        t_terminate = t;

        std::printf("\n");
        std::printf(" *** Intercept: %s   Time = %.3f sec ***\n", label, t);
        std::printf("       alt = %.1f m  speed = %.2f m/s\n", alt, dvbe);
        std::printf("       heading = %.2f deg  flight-path = %.2f deg\n",
                    psivdx, thtvdx);
        std::printf("       Mach = %.3f  q = %.2f Pa\n\n", vmach, pdyn);

        // Signal clean termination
        osk::Sim::stop = -1;
    }

    // ---- Update apogee-detection memory at step boundaries only ----
    // Updating every sub-stage would cause us to compare against the
    // sub-stage-3 value at the next sub-stage-0, which is unstable.
    if (osk::State::stepstart) {
        vz_prev_ = vz_up;
    }
}

} // namespace rocket6dof
