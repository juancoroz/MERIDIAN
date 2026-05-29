//  environment.h  --  Atmosphere + gravity block for the 6DOF rocket
//
//  Public members match the "out"-tagged module variables of Zipfel,
//  "Modeling and Simulation of Aerospace Vehicle Dynamics" 3rd ed.,
//  Section 10.3.2 (round6[50-99]).  See environment.cpp for naming
//  cross-reference.
//
//  This is a pure feed-forward block: no integrator states are
//  registered.  Wind smoothing and Dryden turbulence (Zipfel's VAEDS,
//  taux1, taux2 states) are deliberately deferred to a later pass; in
//  v1 the wind output VAED is identically zero and dvba == |VBED|.
//
//  Inputs are PULLED from Newton via getsFrom(); they may alternatively
//  be set directly on the public members for tests that don't want a
//  Newton instance in the pipeline.

#ifndef ROCKET6DOF_ENVIRONMENT_H
#define ROCKET6DOF_ENVIRONMENT_H

#include "../osk/osk.h"

namespace rocket6dof {

class Newton;   // optional source for alt, SBII, VBED

class Environment : public osk::Block {
public:
    // ---- Atmosphere model selector (Zipfel 'matmo' digit of 'mair') ----
    //   0 = US 1976 Standard Atmosphere, 0 to 86 km
    int matmo;

    // ---- Gravity model selector ----
    //   0 = simple inverse-square, g = mu / r^2 along -r_hat
    int mgrav;

    // ---- Optional Newton input (preferred over direct member assignment) ----
    Newton* newton;
    void getsFrom(Newton* n) { newton = n; }

    // ---- Inputs (read each update; refreshed from Newton if wired) ----
    double   alt;     // [m]     geometric altitude above mean sea level
    osk::Vec SBII;    // [m]     vehicle inertial position from Earth center
    osk::Vec VBED;    // [m/s]   vehicle velocity in geodetic frame

    // ---- Outputs (Zipfel Section 10.3.2 names) ----
    double press;     // [Pa]    atmospheric pressure
    double rho;       // [kg/m^3] density
    double tempk;     // [K]     temperature
    double tempc;     // [degC]  temperature (diagnostic)
    double vsound;    // [m/s]   speed of sound
    double vmach;     // [-]     vehicle Mach number = |dvba| / vsound
    double pdynmc;    // [Pa]    dynamic pressure = 0.5 * rho * dvba^2
    double grav;      // [m/s^2] gravity magnitude
    osk::Vec GRAVG;   // [m/s^2] gravity vector in geocentric inertial frame
    osk::Vec VAED;    // [m/s]   wind velocity in geodetic frame (zero in v1)
    double dvba;      // [m/s]   vehicle speed wrt the air mass

    Environment();
    void init()   override;
    void update() override;
    void rpt()    override;

    // ---- Getters for downstream blocks (Aerodynamics, Forces, Newton) ----
    ACCESS_FN(double,   rho)
    ACCESS_FN(double,   press)
    ACCESS_FN(double,   tempk)
    ACCESS_FN(double,   vsound)
    ACCESS_FN(double,   vmach)
    ACCESS_FN(double,   pdynmc)
    ACCESS_FN(double,   grav)
    ACCESS_FN(double,   dvba)
    ACCESS_FN(osk::Vec, GRAVG)
    ACCESS_FN(osk::Vec, VAED)
};

// ---- Standalone helpers (also usable outside the Block) ----

// US 1976 Standard Atmosphere, valid for 0 <= h <= 86 km geometric.
// Writes density [kg/m^3], pressure [Pa], temperature [K] for the
// requested geometric altitude [m].  Returns 0 on success, nonzero
// if altitude is outside the supported range (caller decides what
// to do; this implementation clamps outputs to the boundary value).
int atmosphere76(double alt_m, double& rho, double& press, double& tempk);

// Inverse-square gravity: returns the acceleration vector at inertial
// position SBII [m].  Magnitude is mu_earth / |SBII|^2, direction is
// -SBII / |SBII| (i.e., toward Earth center).
osk::Vec grav_inv_sq(const osk::Vec& SBII);

} // namespace rocket6dof

#endif
