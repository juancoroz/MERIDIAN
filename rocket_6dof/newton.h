//  newton.h  --  Translational equations of motion for the 6DOF rocket
//
//  Public members match the "out"-tagged module variables of Zipfel,
//  "Modeling and Simulation of Aerospace Vehicle Dynamics" 3rd ed.,
//  Section 10.3.x (round6[210-299]).  See newton.cpp for the full
//  cross-reference table.
//
//  Integrates inertial position SBII and inertial velocity VBII as
//  6 scalar states under the kernel's default RK4.  The acceleration
//  derivative is built from body-frame applied force (FAPB), the
//  attitude DCM (TBI) from Kinematics, and the gravity vector (GRAVG)
//  from Environment.
//
//  v1 simplifications relative to Zipfel:
//    * Spherical, non-rotating Earth.  WEII = 0, TDI = TGI = identity.
//      This means VBED = VBII (geodetic == inertial frame).  Geodetic
//      lat/lon/alt are computed by spherical conversion from SBII.
//    * No mass-flow integration (vmass is supplied externally as a
//      time-varying scalar and read as-is).
//    * 'mfreeze' autopilot-freeze logic from Zipfel is omitted.

#ifndef ROCKET6DOF_NEWTON_H
#define ROCKET6DOF_NEWTON_H

#include "../osk/osk.h"

namespace rocket6dof {

class Environment;
class Kinematics;
class Forces;        // forward; provides FAPB
class MassProps;     // forward; provides vmass (fallback)
class Propulsion;    // forward; provides vmass (preferred when set)

class Newton : public osk::Block {
public:
    // ---- Earth model selector ----
    //   0 = spherical, non-rotating (v1; closed-form verification friendly)
    //   1 = WGS84 ellipsoid + Earth rotation (matches Zipfel's CADAC)
    // When mearth=1, this block expects Environment->mgrav=1 as well, so
    // that GRAVG comes back in geocentric coords; Newton will rotate via
    // ~TGI to inertial.
    int mearth;

    // ---- Inputs (read each update from upstream blocks) ----
    Environment* env;
    Kinematics*  kin;
    Forces*      forces;
    MassProps*   mass;        // optional, constant-mass fallback
    Propulsion*  prop;        // optional, time-varying mass (preferred if non-null)

    void getsFrom(Environment* e, Kinematics* k,
                  Forces* f, MassProps* m) {
        env = e; kin = k; forces = f; mass = m; prop = nullptr;
    }
    void getsFrom(Environment* e, Kinematics* k,
                  Forces* f, Propulsion* p) {
        env = e; kin = k; forces = f; mass = nullptr; prop = p;
    }

    // ---- States (integrated by the kernel) ----
    //   SBII : inertial position vector       [m]
    //   VBII : inertial velocity vector       [m/s]
    // Stored as 3 scalars each so they can be addIntegrator'd.
    double sx, sy, sz;          // SBII components
    double vx, vy, vz;          // VBII components
    double sxd, syd, szd;       // SBII derivatives = VBII
    double vxd, vyd, vzd;       // VBII derivatives = ABII (inertial accel)

    // ---- Initial-condition parameters (Zipfel "data" tag) ----
    double lonx0;   // [deg] initial longitude
    double latx0;   // [deg] initial latitude
    double alt0;    // [m]   initial altitude above mean spherical Earth
    double dvbe0;   // [m/s] initial geographic speed magnitude
    double psivdx0; // [deg] initial heading (azimuth of velocity, 0=N, 90=E)
    double thtvdx0; // [deg] initial flight-path angle (0=horizontal, +90=up)

    // ---- Outputs (Zipfel Section 10.3.x names) ----
    osk::Vec SBII;     // [m]    inertial position
    osk::Vec VBII;     // [m/s]  inertial velocity
    osk::Vec ABII;     // [m/s^2] inertial acceleration (saved for reporting)
    osk::Vec VBED;     // [m/s]  geographic velocity, geodetic frame
    osk::Vec FSPB;     // [m/s^2] specific force in body frame = FAPB/vmass
    double   dbi;      // [m]    distance from Earth center = |SBII|
    double   alt;      // [m]    altitude above mean spherical Earth
    double   altx;     // [kft]  altitude in thousand feet (diagnostic)
    double   lonx;     // [deg]  longitude
    double   latx;     // [deg]  latitude
    double   dvbe;     // [m/s]  geographic speed magnitude
    double   dvbi;     // [m/s]  inertial speed magnitude
    double   psivdx;   // [deg]  velocity heading
    double   thtvdx;   // [deg]  velocity flight-path angle
    double   axx;      // [g]    body-axial acceleration (FSPB.x / g0)
    double   ayx;      // [g]    body-lateral acceleration
    double   anx;      // [g]    body-normal acceleration (-FSPB.z / g0)

    Newton();
    void init()   override;
    void update() override;
    void rpt()    override;

    // ---- Getters for downstream blocks ----
    ACCESS_FN(osk::Vec, SBII)
    ACCESS_FN(osk::Vec, VBII)
    ACCESS_FN(osk::Vec, VBED)
    ACCESS_FN(osk::Vec, FSPB)
    ACCESS_FN(double,   alt)
    ACCESS_FN(double,   dvbi)
    ACCESS_FN(double,   dvbe)
    ACCESS_FN(double,   dbi)
};

} // namespace rocket6dof

#endif
