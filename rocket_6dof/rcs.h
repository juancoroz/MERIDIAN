//  rcs.h  --  Reaction Control System (vacuum-phase attitude actuator)
//
//  The RCS produces body-frame torques and side forces via small
//  thrusters firing.  In atmosphere, TVC and aero deflections handle
//  attitude control; out of atmosphere, those aren't effective and the
//  vehicle needs RCS for any attitude maneuver.
//
//  Two independent capabilities:
//
//    Attitude control (mrcs_moment = 10*rcs_type + rcs_mode):
//        rcs_type = 0  off
//                 = 1  proportional thrusters (continuous Newton output)
//                 = 2  on-off thrusters (Schmitt trigger discrete output)
//        rcs_mode = 0  no control
//                 = 1  geodetic Euler-angle control (phi, theta, psi)
//                 = 2  thrust-vector direction + roll control (uses UTBC)
//                 = 3  incidence + roll control (uses alpha, beta)
//
//    Side-force control (mrcs_force):
//        = 0 off
//        = 1 proportional side thrusters
//        = 2 on-off (Schmitt) side thrusters
//
//  Outputs:
//      FMRCS = 3-vector, RCS-induced body moment [Nm]
//      FARCS = 3-vector, RCS-induced body force  [N]
//  These add to Forces' aggregated FAPB/FMB.
//
//  Cross-reference: Zipfel rcs.cpp (hyper[50-93]).

#ifndef ROCKET6DOF_RCS_H
#define ROCKET6DOF_RCS_H

#include "../osk/osk.h"

namespace rocket6dof {

class Propulsion;
class Newton;
class Kinematics;
class INS;
class Guidance;

class RCS : public osk::Block {
public:
    // ---- Inputs ----
    Propulsion* prop;       // for IBBB (vehicle moment of inertia)
    Newton*     newton;     // for alpha, beta (incidence)
    Kinematics* kin;        // for body Euler angles (truth)
    INS*        ins;        // for body rates and specific force
    Guidance*   guidance;   // for UTBC (thrust direction command)

    void getsFrom(Propulsion* p, Newton* n, Kinematics* k, INS* i, Guidance* g) {
        prop = p; newton = n; kin = k; ins = i; guidance = g;
    }

    // ---- Mode selectors ----
    int mrcs_moment;       // 10*rcs_type + rcs_mode
    int mrcs_force;        // 0/1/2

    // ---- Proportional-mode parameters ----
    double rcs_freq;       // [rad/s] closed-loop natural frequency
    double rcs_zeta;       // damping ratio

    // ---- Schmitt-trigger parameters ----
    double dead_zone;      // [deg or N]  dead band of trigger
    double hysteresis;     // [deg or N]
    double rcs_tau;        // [s] switching-function slope

    // ---- Saturation limits ----
    double roll_mom_max;   // [Nm]
    double pitch_mom_max;  // [Nm]
    double yaw_mom_max;    // [Nm]
    double side_force_max; // [N]
    double acc_gain;       // [N/(m/s^2)] side-force gain

    // ---- Attitude commands (set externally) ----
    double phibdcomx;      // [deg] roll-angle command
    double thtbdcomx;      // [deg] pitch-angle command
    double psibdcomx;      // [deg] yaw-angle command
    double alphacomx;      // [deg] AoA command (mode 3)
    double betacomx;       // [deg] sideslip command (mode 3)
    double aycomx;         // [g]   right-force command
    double azcomx;         // [g]   down-force command

    // ---- Internal Schmitt-trigger state ----
    double roll_save, pitch_save, yaw_save;
    double right_save, down_save;
    int    o_roll, o_pitch, o_yaw;
    int    o_right, o_down;
    int    roll_count, pitch_count, yaw_count;
    int    right_count, down_count;

    // ---- Outputs ----
    osk::Vec FMRCS;        // [Nm]
    osk::Vec FARCS;        // [N]

    // ---- Diagnostics ----
    double e_roll, e_pitch, e_yaw;
    double e_right, e_down;

    RCS();
    void init()   override;
    void update() override;
    void rpt()    override;

    ACCESS_FN(osk::Vec, FMRCS)
    ACCESS_FN(osk::Vec, FARCS)

private:
    // Proportional thruster with saturation
    static double rcs_prop(double input, double limit);
    // Schmitt trigger: dead-zone + hysteresis switching
    static int    rcs_schmitt(double input_new, double input_prev,
                              double dz, double hyst);
};

} // namespace rocket6dof

#endif
