//  kinematics.h  --  Rotational kinematics for the 6DOF rocket
//
//  Public members match the "out"-tagged module variables of Zipfel,
//  "Modeling and Simulation of Aerospace Vehicle Dynamics" 3rd ed.,
//  Section 10.3.x (round6[100-149]).  See kinematics.cpp for the full
//  cross-reference table.
//
//  Integrates the body-from-inertial DCM TBI as 9 scalar states, then
//  orthonormalizes once per step using Bar-Itzhack's first-order
//  correction.  Reads body angular velocity WBIB from the Euler block.
//
//  DCM derivative (Zipfel Eq. 3.46):
//      d/dt TBI = -[WBIB]_x * TBI
//  where [w]_x is the skew-symmetric matrix of w.
//
//  v1 simplification: non-rotating Earth, so TDI = identity and
//  TBD = TBI.  When we add Earth rotation later, TDI will become
//  time-dependent and TBD = TBI * TDI^T.

#ifndef ROCKET6DOF_KINEMATICS_H
#define ROCKET6DOF_KINEMATICS_H

#include "../osk/osk.h"

namespace rocket6dof {

class Environment;
class Newton;
class Euler;     // forward; provides WBIB

class Kinematics : public osk::Block {
public:
    // ---- Earth model selector ----
    //   0 = spherical, non-rotating (TDI = identity, TBD = TBI)
    //   1 = WGS84 + Earth rotation (TBD = TBI * ~TDI(t))
    // Should match Newton::mearth and Environment::mgrav for consistency.
    int mearth;

    // ---- Inputs ----
    Environment* env;     // for VAED (wind) used in incidence-angle calc
    Newton*      newton;  // for VBED, VBII, dvba-related quantities
    Euler*       euler;   // for WBIB

    void getsFrom(Environment* e, Newton* n, Euler* eu) {
        env = e; newton = n; euler = eu;
    }

    // ---- Initial-condition parameters ----
    double psibdx0;  // [deg] initial yaw (3-2-1 Euler psi)
    double thtbdx0;  // [deg] initial pitch (theta)
    double phibdx0;  // [deg] initial roll (phi)

    // ---- States: TBI integrated as 9 scalars, row-major ----
    double t00, t01, t02;
    double t10, t11, t12;
    double t20, t21, t22;
    double t00d, t01d, t02d;
    double t10d, t11d, t12d;
    double t20d, t21d, t22d;

    // ---- Outputs (Zipfel Section 10.3.x names) ----
    osk::Mat TBI;       // body from inertial DCM
    osk::Mat TBD;       // body from geodetic DCM
    double   ortho_error;   // orthogonality residual after correction
    double   psibd, thtbd, phibd;     // Euler angles, rad
    double   psibdx, thtbdx, phibdx;  // Euler angles, deg
    double   alphax;    // angle of attack [deg]
    double   betax;     // sideslip [deg]
    double   alppx;     // total angle of attack [deg]
    double   phipx;     // aerodynamic roll [deg]
    double   alphaix;   // inertial alpha [deg]
    double   betaix;    // inertial beta [deg]

    Kinematics();
    void init()   override;
    void update() override;
    void rpt()    override;

    // ---- Getters for downstream blocks ----
    ACCESS_FN(osk::Mat, TBI)
    ACCESS_FN(osk::Mat, TBD)
    ACCESS_FN(double,   alphax)
    ACCESS_FN(double,   betax)
    ACCESS_FN(double,   alppx)
    ACCESS_FN(double,   phipx)
    ACCESS_FN(double,   thtbdx)
    ACCESS_FN(double,   psibdx)
    ACCESS_FN(double,   phibdx)
};

} // namespace rocket6dof

#endif
