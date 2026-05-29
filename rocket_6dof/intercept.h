//  intercept.h  --  Trajectory termination conditions
//
//  Watches the simulation state for termination conditions and signals
//  Sim::stop = -1 when any condition fires.  Modeled on Zipfel's
//  intercept.cpp (Section 10.3.x, hyper[650-699]); we extend it with
//  apogee and aero-decay conditions since this is a sounding-rocket
//  sim, not a homing missile.
//
//  Supported conditions (each opt-in via its own flag):
//
//    1. Ground impact      (alt <= 0)                  -- on by default
//    2. Apogee             (vertical velocity crosses 0 downward)
//    3. Max altitude       (alt > alt_max)
//    4. Min Mach           (vmach < trmach)            -- aero decay
//    5. Min dynamic pres   (pdynmc < trdynm)           -- aero decay
//    6. Max flight time    (t > t_max)
//
//  When a condition fires, a one-line summary is printed to stdout and
//  Sim::stop = -1 is written.  The 'trcond' member records which
//  condition fired (1..6 above); 0 = still running.  The 'fired_' flag
//  is a one-shot to prevent re-printing on successive updates.

#ifndef ROCKET6DOF_INTERCEPT_H
#define ROCKET6DOF_INTERCEPT_H

#include "../osk/osk.h"

namespace rocket6dof {

class Environment;
class Newton;
class Euler;

class Intercept : public osk::Block {
public:
    // ---- Inputs ----
    Environment* env;
    Newton*      newton;
    Euler*       eul;        // optional: needed only if check_tumble is on
    void getsFrom(Environment* e, Newton* n) { env = e; newton = n; eul = nullptr; }
    void getsFrom(Environment* e, Newton* n, Euler* eu) { env = e; newton = n; eul = eu; }

    // ---- Condition selectors (per-condition opt-in flags) ----
    bool   check_ground_impact;  // alt <= 0                  default: true
    bool   check_apogee;         // vertical velocity sign change
    bool   check_alt_max;        // alt > alt_max
    bool   check_aero_decay;     // vmach < trmach && pdynmc < trdynm
    bool   check_time_max;       // t > t_max
    bool   check_tumble;         // |WBEB| > rate_max         default: false

    // ---- Threshold parameters ----
    double alt_max;       // [m]    max altitude before termination
    double trmach;        // [-]    min Mach before aero-decay termination
    double trdynm;        // [Pa]   min dyn pressure (paired with trmach)
    double t_max;         // [s]    max flight time (independent of Sim's tmax)
    double rate_max;      // [rad/s] tumble threshold; default 3.49 (~200 deg/s)

    // ---- Suppress-on-init guard ----
    // The apogee detector needs a previous vz value.  We don't want it
    // to fire spuriously when the rocket starts at ground level (alt=0)
    // with zero velocity (alt<=0 ground-impact would trigger immediately).
    // 't_no_check_until' delays all checks until simulation time >= this
    // value.  Default = 0.0 (check from the start).
    double t_no_check_until;

    // ---- Outputs ----
    int    trcond;            // 0=still running, 1=ground, 2=apogee,
                              // 3=max alt, 4=aero decay, 5=time,
                              // 6=tumble
    double t_terminate;       // [s]   time of termination

    Intercept();
    void init()   override;
    void update() override;
    void rpt()    override {}    // no periodic reporting

private:
    bool   fired_;             // one-shot guard
    double vz_prev_;           // previous step's vertical velocity (= -VBED.z)
};

} // namespace rocket6dof

#endif
