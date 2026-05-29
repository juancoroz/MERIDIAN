//  euler.h  --  Rotational equations of motion for the 6DOF rocket
//
//  Public members match the "out"-tagged module variables of Zipfel,
//  "Modeling and Simulation of Aerospace Vehicle Dynamics" 3rd ed.,
//  Section 10.3.x (round6[150-199]).  See euler.cpp for the full
//  cross-reference table.
//
//  Integrates the body angular velocity WBIB (3 scalar states) under
//  the kernel's default RK4.  The derivative is:
//
//      WBIBD = I^-1 * (FMB - omega x (I*omega))
//
//  where I = IBBB (inertia tensor in body frame), FMB is the applied
//  moment in body frame, and omega = WBIB.  This is Euler's rotational
//  equation solved for angular acceleration.
//
//  Reads FMB from Forces, IBBB from MassProps, TBI from Kinematics.
//  Writes WBIB which Kinematics reads back to drive the TBI integration.
//  Kinematics MUST appear before Euler in the stage vector so TBI is
//  fresh when init/update of Euler runs.

#ifndef ROCKET6DOF_EULER_H
#define ROCKET6DOF_EULER_H

#include "../osk/osk.h"

namespace rocket6dof {

class Forces;
class MassProps;
class Kinematics;
class Propulsion;

class Euler : public osk::Block {
public:
    // ---- Earth model selector ----
    //   0 = non-rotating Earth (WEII = 0); WBIB = WBEB
    //   1 = WGS84 rotating Earth (WBIB = WBEB + TBI*WEII)
    // Should match Newton::mearth, Kinematics::mearth, and
    // Environment::mgrav for consistency across the simulation.
    int mearth;

    // ---- Inputs ----
    Forces*     forces;
    MassProps*  mass;
    Kinematics* kin;
    Propulsion* prop;       // optional; takes priority over mass for IBBB

    void getsFrom(Forces* f, MassProps* m, Kinematics* k) {
        forces = f; mass = m; kin = k; prop = nullptr;
    }
    void getsFrom(Forces* f, Propulsion* p, Kinematics* k) {
        forces = f; mass = nullptr; kin = k; prop = p;
    }

    // ---- Initial-condition parameters (Zipfel input names) ----
    double ppx0;   // [deg/s] initial body roll rate  wrt Earth
    double qqx0;   // [deg/s] initial body pitch rate wrt Earth
    double rrx0;   // [deg/s] initial body yaw rate   wrt Earth

    // ---- States (integrated by kernel RK4) ----
    //   WBIB = body angular velocity wrt INERTIAL frame, body-axis components.
    // Stored as 3 scalars so they can be addIntegrator'd.
    double wx, wy, wz;       // WBIB components
    double wxd, wyd, wzd;    // WBIBD = derivative of WBIB

    // ---- Outputs (Zipfel Section 10.3.x names) ----
    osk::Vec WBIB;      // [rad/s] body angular velocity wrt inertial, body coords
    osk::Vec WBIBD;     // [rad/s^2] angular acceleration, body coords
    osk::Vec WBEB;      // [rad/s] body angular velocity wrt Earth, body coords
    osk::Vec WBII;      // [rad/s] body angular velocity wrt inertial, INERTIAL coords
    double   ppx;       // [deg/s] WBEB.x (roll rate)
    double   qqx;       // [deg/s] WBEB.y (pitch rate)
    double   rrx;       // [deg/s] WBEB.z (yaw rate)

    Euler();
    void init()   override;
    void update() override;
    void rpt()    override;

    // ---- Getters for downstream blocks (Kinematics reads WBIB) ----
    ACCESS_FN(osk::Vec, WBIB)
    ACCESS_FN(osk::Vec, WBII)
    ACCESS_FN(osk::Vec, WBEB)
};

} // namespace rocket6dof

#endif
