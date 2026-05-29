//  ins.h  --  Inertial Navigation System block
//
//  Public members match Zipfel's ins.cpp (Section 10.3.x,
//  hyper[300-399]).  The INS reads truth state from the physics blocks,
//  adds sensor errors, and produces "computed" state that downstream
//  Control / Guidance can use INSTEAD of truth.  This more closely
//  models a real vehicle, where the autopilot only ever sees sensor
//  estimates (never true state).
//
//  Modes (mins):
//    0  ideal INS -- outputs identically equal truth (no noise, no bias).
//                    This is the default and the verification-friendly
//                    mode.  Useful to wire Control through INS without
//                    perturbing existing test results.
//    1  simple-error INS -- constant accelerometer bias and gyro bias
//                    are added to measurements; INS integrates the
//                    perturbed measurements to produce position,
//                    velocity, and attitude that drift over time.
//
//  Future modes (deferred):
//    2  Full Zipfel error model -- tilt error states, position error
//                    states, random walk, scale factor, misalignment.
//                    Requires GPS / StarTrack hookups for correction.
//
//  Cross-reference to Zipfel (selected variables):
//      Zipfel name      OSK member       Meaning
//      -----------      ----------       -------
//      TBIC             TBIC             computed body-from-inertial DCM
//      WBICB            WBICB            computed body rate, body coords [rad/s]
//      WBICI            WBICI            computed body rate, inertial coords
//      FSPCB            FSPCB            computed specific force, body [m/s^2]
//      SBIIC            SBIIC            computed inertial position [m]
//      VBIIC            VBIIC            computed inertial velocity [m/s]
//      qqcx,rrcx,ppcx   qqcx,rrcx,ppcx   computed body rates [deg/s]
//      EBIASA           bias_accel       accelerometer bias [m/s^2]
//      EBIASG           bias_gyro        gyro bias [rad/s]

#ifndef ROCKET6DOF_INS_H
#define ROCKET6DOF_INS_H

#include "../osk/osk.h"

namespace rocket6dof {

class Newton;
class Euler;
class Kinematics;
class GPS;
class Startrack;

class INS : public osk::Block {
public:
    // ---- Inputs (truth from physics blocks) ----
    Newton*     newton;
    Euler*      euler;
    Kinematics* kin;
    GPS*        gps;          // optional; provides discrete position/velocity
                              // updates that correct INS drift
    Startrack*  startrack;    // optional; provides attitude tilt corrections
    void getsFrom(Newton* n, Euler* e, Kinematics* k) {
        newton = n; euler = e; kin = k; gps = nullptr; startrack = nullptr;
    }
    void getsFrom(Newton* n, Euler* e, Kinematics* k, GPS* g) {
        newton = n; euler = e; kin = k; gps = g; startrack = nullptr;
    }
    void getsFrom(Newton* n, Euler* e, Kinematics* k, GPS* g, Startrack* s) {
        newton = n; euler = e; kin = k; gps = g; startrack = s;
    }

    // ---- Mode selector ----
    int mins;             // 0 = ideal, 1 = constant bias

    // ---- GPS-update tracking ----
    int gps_update_count;       // # of GPS updates INS has consumed
    int startrack_update_count; // # of star-track updates INS has consumed

    // ---- Sensor error parameters (used when mins=1) ----
    osk::Vec bias_accel;  // [m/s^2]  accelerometer bias, body axes
    osk::Vec bias_gyro;   // [rad/s]  gyro bias, body axes

    // ---- States (integrator-managed; used when mins=1) ----
    // These accumulate over time, drifting INS away from truth.
    osk::Vec SBIIC;       // [m]      computed inertial position
    osk::Vec VBIIC;       // [m/s]    computed inertial velocity
    osk::Mat TBIC;        // [-]      computed body-from-inertial DCM

    // Integrator scalars (for OSK addIntegrator)
    double sxc, syc, szc;       // SBIIC components (state)
    double sxcd, sycd, szcd;    // dot
    double vxc, vyc, vzc;       // VBIIC components (state)
    double vxcd, vycd, vzcd;    // dot
    // TBIC stored as 9 scalars (like Kinematics)
    double t00c, t01c, t02c, t10c, t11c, t12c, t20c, t21c, t22c;
    double t00cd, t01cd, t02cd, t10cd, t11cd, t12cd, t20cd, t21cd, t22cd;

    // ---- Outputs (read by Control / Guidance) ----
    osk::Vec FSPCB;       // [m/s^2]  computed specific force, body coords
    osk::Vec WBICB;       // [rad/s]  computed body rate, body coords
    osk::Vec WBICI;       // [rad/s]  computed body rate, inertial coords
    double   ppcx;        // [deg/s]  computed roll rate (Earth-relative approx)
    double   qqcx;        // [deg/s]  computed pitch rate
    double   rrcx;        // [deg/s]  computed yaw rate

    // Position/velocity diagnostics
    double dbic;          // [m]      computed |SBIIC|
    double dvbec;         // [m/s]    computed geographic speed
    double ins_pos_err;   // [m]      |SBIIC - SBII_truth|
    double ins_vel_err;   // [m/s]    |VBIIC - VBII_truth|
    double ins_att_err;   // [rad]    small-angle magnitude of (TBIC vs TBI)

    INS();
    void init()   override;
    void update() override;
    void rpt()    override;

    // Getters for Control
    ACCESS_FN(osk::Vec, FSPCB)
    ACCESS_FN(osk::Vec, WBICB)
    ACCESS_FN(double,   qqcx)
    ACCESS_FN(double,   rrcx)
    ACCESS_FN(double,   ppcx)
    ACCESS_FN(double,   dvbec)
};

} // namespace rocket6dof

#endif
