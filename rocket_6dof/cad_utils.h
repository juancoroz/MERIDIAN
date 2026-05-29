//  cad_utils.h  --  CADAC utility functions, ported to osk types
//
//  Source: Zipfel, "Modeling and Simulation of Aerospace Vehicle
//  Dynamics" 3rd ed., utility_functions.cpp.  Only the subset of
//  functions actually called by Newton, Kinematics, and Environment
//  is ported here.  The naming and equations are verbatim from
//  Zipfel; only the Matrix container is changed (Matrix(3,1) -> Vec,
//  Matrix(3,3) -> Mat).
//
//  Constants are WGS84 standard values.  Sign of C20 follows Zipfel:
//  C20 = -J2 / sqrt(5).

#ifndef ROCKET6DOF_CAD_UTILS_H
#define ROCKET6DOF_CAD_UTILS_H

#include "../osk/osk.h"

namespace rocket6dof {

// ---- WGS84 / Earth-rotation constants ----
// WGS84 standard values.  The flattening uses the IUGG/WGS84 reference
// value 3.3528106647e-3.  This is a deliberate departure from Zipfel's
// global_constants.hpp, which prints 3.33528106e-3 (a duplicated '3'
// in the mantissa, ~0.5% below the true value).  Other constants here
// match Zipfel byte-for-byte.  Sign of C20 follows Zipfel:
// C20 = -J2 / sqrt(5).
constexpr double WGS84_SMAJOR_AXIS = 6378137.0;          // [m]   equatorial radius
constexpr double WGS84_FLATTENING  = 3.3528106647e-3;    // [-]   IUGG/WGS84 standard
constexpr double WGS84_GM          = 3.9860044e14;       // [m^3/s^2]  (Zipfel uses 7-digit value)
constexpr double WGS84_WEII3       = 7.292115e-5;        // [rad/s] Earth spin rate (Zipfel's 7-digit value)
constexpr double WGS84_C20         = -4.8416685e-4;      // [-] normalized zonal harmonic
constexpr double REARTH_SPHERICAL  = 6370987.308;        // [m]   mean radius (Zipfel's value)
constexpr double GW_CLONG          = 0.0;                // [rad] Greenwich celestial longitude
                                                          //       at t=0.
constexpr double AGRAV             = 9.80675445;         // [m/s^2] standard gravity (Zipfel)

// Convergence threshold for iterative geodetic-latitude solve.
// Matches Zipfel's SMALL = 1e-7 from global_constants.hpp.
constexpr double CAD_SMALL = 1.0e-7;

// ---- Function prototypes ----

// Inertial position from geodetic (lon, lat, alt) using WGS84 ellipsoid.
osk::Vec cad_in_geo84(double lon, double lat, double alt, double time);

// Geodetic (lon, lat, alt) from inertial position.  Iterative.
void cad_geo84_in(double& lon, double& lat, double& alt,
                  const osk::Vec& SBII, double time);

// Geocentric (lon, lat, alt) from inertial position.  Spherical Earth.
void cad_geoc_in(double& lonc, double& latc, double& altc,
                 const osk::Vec& SBII, double time);

// WGS84 gravity vector in geographic (geocentric) frame.
//   GRAVG.x = lateral J2 component
//   GRAVG.y = 0 (axisymmetric)
//   GRAVG.z = primary radial inverse-square + J2 correction
osk::Vec cad_grav84(const osk::Vec& SBII, double time);

// DCM of geodetic (D-frame, NED at vehicle's lon/lat) wrt inertial (I).
osk::Mat cad_tdi84(double lon, double lat, double alt, double time);

// DCM of geographic/geocentric (G-frame) wrt inertial.
osk::Mat cad_tgi84(double lon, double lat, double alt, double time);

// DCM of Earth-fixed (E-frame) wrt inertial.
osk::Mat cad_tei(double time);

// Geodetic velocity diagnostics: dvbe, heading, flight-path angle.
void cad_geo84vel_in(double& dvbe, double& psivdx, double& thtvdx,
                     const osk::Vec& SBII, const osk::Vec& VBII,
                     double time);

// Project an inertial state forward in time along a Keplerian (two-body)
// trajectory.  Iteratively solves the universal Kepler equation.
//   SPII (out) = projected inertial position after tgo seconds [m]
//   VPII (out) = projected inertial velocity after tgo seconds [m/s]
//   SBII (in)  = initial inertial position [m]
//   VBII (in)  = initial inertial velocity [m/s]
//   tgo  (in)  = time-of-flight [s]
//   returns 0 = converged, 1 = failed to converge (SPII, VPII unchanged)
int cad_kepler(osk::Vec& SPII, osk::Vec& VPII,
               const osk::Vec& SBII, const osk::Vec& VBII,
               double tgo);

// Angle between two 3-vectors in radians (clamped to [0, pi]).  Returns 0
// if either vector has zero magnitude.  Ported from Zipfel's utility_functions.cpp.
double angle(const osk::Vec& v1, const osk::Vec& v2);

} // namespace rocket6dof

#endif
