//  newton.cpp  --  Translational EOM for the 6DOF rocket
//
//  Cross-reference to Zipfel "Modeling and Simulation of Aerospace
//  Vehicle Dynamics" 3rd ed., Section 10.3.x (round6[210-299]):
//
//      Zipfel name      OSK member      Meaning
//      -----------      ----------      -------
//      SBII             SBII            inertial position        [m]
//      VBII             VBII            inertial velocity        [m/s]
//      ABII             ABII            inertial acceleration    [m/s^2]
//      VBED             VBED            geographic velocity      [m/s]
//      FSPB             FSPB            specific force, body     [m/s^2]
//      dbi              dbi             |SBII|                   [m]
//      alt              alt             altitude                 [m]
//      lonx, latx       lonx, latx      lon, lat                 [deg]
//      dvbe             dvbe            |VBED|                   [m/s]
//      dvbi             dvbi            |VBII|                   [m/s]
//      psivdx, thtvdx   psivdx, thtvdx  heading, flight-path     [deg]
//      axx, ayx, anx    axx, ayx, anx   body accels              [g]
//
//  EOM (Zipfel Eq. 3.21):
//      d/dt VBII = (1/m) * T_IB * FAPB + GRAVG
//      d/dt SBII = VBII
//  where T_IB = ~TBI (transpose of body-from-inertial DCM).

#include "newton.h"
#include "environment.h"
#include "kinematics.h"
#include "forces.h"
#include "mass_props.h"
#include "propulsion.h"
#include "cad_utils.h"
#include <cmath>
#include <cstdio>

namespace rocket6dof {

namespace {
constexpr double R_EARTH = 6378137.0;     // [m] mean spherical Earth radius
constexpr double G0      = 9.80665;       // [m/s^2] standard gravity (for g-unit reporting)
constexpr double FOOT    = 3.280839895;   // m -> ft
constexpr double DEG     = 180.0 / osk::PI;
constexpr double RAD     = osk::PI / 180.0;
} // anon

Newton::Newton()
    : env(nullptr), kin(nullptr), forces(nullptr), mass(nullptr), prop(nullptr)
{
    mearth = 0;   // spherical non-rotating Earth by default

    // Initial-condition defaults: stationary on the equator at sea level
    lonx0 = 0.0;  latx0 = 0.0;  alt0 = 0.0;
    dvbe0 = 0.0;  psivdx0 = 0.0;  thtvdx0 = 0.0;

    sx = sy = sz = 0.0;
    vx = vy = vz = 0.0;
    sxd = syd = szd = 0.0;
    vxd = vyd = vzd = 0.0;

    SBII = VBII = ABII = VBED = FSPB = osk::Vec(0,0,0);
    dbi = alt = altx = 0.0;
    lonx = latx = 0.0;
    dvbe = dvbi = psivdx = thtvdx = 0.0;
    axx = ayx = anx = 0.0;

    // Register the six scalar integrator states.  Kernel uses default RK4.
    addIntegrator(sx, sxd);
    addIntegrator(sy, syd);
    addIntegrator(sz, szd);
    addIntegrator(vx, vxd);
    addIntegrator(vy, vyd);
    addIntegrator(vz, vzd);
}

void Newton::init() {
    if (initCount == 0) {
        double lat = latx0 * RAD;
        double lon = lonx0 * RAD;
        double psi   = psivdx0 * RAD;
        double theta = thtvdx0 * RAD;

        // ---- Initial inertial position ----
        if (mearth == 1) {
            // WGS84 ellipsoid + Earth rotation
            osk::Vec S0 = cad_in_geo84(lon, lat, alt0, /*time=*/0.0);
            sx = S0.x;  sy = S0.y;  sz = S0.z;
        } else {
            // Spherical, non-rotating: x toward (lon=0,lat=0), z along spin axis
            double r = R_EARTH + alt0;
            sx = r * std::cos(lat) * std::cos(lon);
            sy = r * std::cos(lat) * std::sin(lon);
            sz = r * std::sin(lat);
        }

        // ---- Initial inertial velocity from heading + flight-path ----
        // VBED in local NED: x=N, y=E, z=Down.
        double vN =  dvbe0 * std::cos(theta) * std::cos(psi);
        double vE =  dvbe0 * std::cos(theta) * std::sin(psi);
        double vD = -dvbe0 * std::sin(theta);

        if (mearth == 1) {
            // Use TDI from cad_utils: VBII = ~TDI * VBED + WEII x SBII
            osk::Mat TDI = cad_tdi84(lon, lat, alt0, /*time=*/0.0);
            osk::Vec VBED0(vN, vE, vD);
            osk::Vec VBEI = TDI.transpose() * VBED0;     // Earth-relative in inertial
            osk::Vec WEII(0.0, 0.0, WGS84_WEII3);
            osk::Vec SBII0(sx, sy, sz);
            osk::Vec VBII0 = VBEI + WEII.cross(SBII0);   // add Earth rotation contribution
            vx = VBII0.x;  vy = VBII0.y;  vz = VBII0.z;
        } else {
            // Spherical non-rotating: NED basis at (lat, lon) into inertial
            double sl = std::sin(lat), cl = std::cos(lat);
            double so = std::sin(lon), co = std::cos(lon);
            vx = vN*(-sl*co) + vE*(-so) + vD*(-cl*co);
            vy = vN*(-sl*so) + vE*( co) + vD*(-cl*so);
            vz = vN*( cl)    + vE*( 0 ) + vD*(-sl);
        }

        SBII = osk::Vec(sx, sy, sz);
        VBII = osk::Vec(vx, vy, vz);
    }
    // On subsequent stage entries (initCount > 0), keep SBII/VBII as they
    // were at the end of the previous stage -- this is the OSK pattern
    // for multi-stage continuity (manual Section 5.3).
}

void Newton::update() {
    // ---- Refresh Vec outputs from the integrator scalars ----
    SBII = osk::Vec(sx, sy, sz);
    VBII = osk::Vec(vx, vy, vz);
    dbi  = SBII.mag();
    dvbi = VBII.mag();

    // ---- Geodetic conversion ----
    double lat, lon;
    if (mearth == 1) {
        // WGS84 ellipsoid + Earth rotation
        cad_geo84_in(lon, lat, alt, SBII, osk::State::t);
    } else {
        // Spherical, non-rotating
        lat = std::asin( sz / (dbi > 1.0 ? dbi : 1.0) );
        lon = std::atan2(sy, sx);
        alt = dbi - R_EARTH;
    }
    lonx = lon * DEG;
    latx = lat * DEG;
    altx = 0.001 * alt * FOOT;

    // ---- Build acceleration ABII ----
    //   ABII = (1/m) * (~TBI) * FAPB + (gravity rotated to inertial)
    //
    // In mearth=1 mode, GRAVG is in geocentric coords (cad_grav84 output)
    // and must be rotated to inertial via ~TGI.  In mearth=0 mode, GRAVG
    // is already in inertial coords (grav_inv_sq returns inertial vector).

    osk::Vec FAPB(0, 0, 0);
    double   vmass = 1.0;
    osk::Mat TBI(1,0,0, 0,1,0, 0,0,1);    // default identity; osk::Mat() is zero!
    osk::Vec GRAVG(0, 0, 0);

    if (forces) FAPB  = forces->FAPB_();
    if (prop)        vmass = prop->vmass_();
    else if (mass)   vmass = mass->vmass_();
    if (kin)    TBI   = kin->TBI_();
    if (env)    GRAVG = env->GRAVG_();

    if (vmass < 1.0e-9) vmass = 1.0e-9;

    FSPB = FAPB * (1.0 / vmass);
    osk::Mat TIB = TBI.transpose();
    osk::Vec FSPI = TIB * FSPB;

    osk::Vec GRAV_I;
    if (mearth == 1) {
        osk::Mat TGI = cad_tgi84(lon, lat, alt, osk::State::t);
        GRAV_I = TGI.transpose() * GRAVG;
    } else {
        GRAV_I = GRAVG;
    }
    ABII = FSPI + GRAV_I;

    // ---- Push derivatives ----
    sxd = vx;  syd = vy;  szd = vz;
    vxd = ABII.x;  vyd = ABII.y;  vzd = ABII.z;

    // ---- Geodetic velocity diagnostic ----
    if (mearth == 1) {
        // VBED = TDI * (VBII - WEII x SBII)
        osk::Mat TDI = cad_tdi84(lon, lat, alt, osk::State::t);
        osk::Vec WEII(0.0, 0.0, WGS84_WEII3);
        VBED = TDI * (VBII - WEII.cross(SBII));
    } else {
        // Non-rotating: VBED = VBII expressed in local NED basis
        double sl = std::sin(lat), cl = std::cos(lat);
        double so = std::sin(lon), co = std::cos(lon);
        double vN = vx*(-sl*co) + vy*(-sl*so) + vz*( cl);
        double vE = vx*(-so   ) + vy*( co   ) + vz*( 0 );
        double vD = vx*(-cl*co) + vy*(-cl*so) + vz*(-sl);
        VBED = osk::Vec(vN, vE, vD);
    }
    dvbe = VBED.mag();

    // Heading and flight-path angle from VBED
    double vN = VBED.x, vE = VBED.y, vD = VBED.z;
    double vh = std::sqrt(vN*vN + vE*vE);
    psivdx = std::atan2(vE, vN) * DEG;
    thtvdx = std::atan2(-vD, vh) * DEG;

    // Body-frame accel diagnostics in g's
    axx =  FSPB.x / G0;
    ayx =  FSPB.y / G0;
    anx = -FSPB.z / G0;
}

void Newton::rpt() {
    if (osk::State::sample(1.0)) {
        std::printf("Newt t=%7.3f  alt=%9.1f  dvbi=%8.2f  dvbe=%8.2f  "
                    "lon=%7.3f  lat=%7.3f  fpa=%6.2f\n",
                    osk::State::t, alt, dvbi, dvbe, lonx, latx, thtvdx);
    }
}

} // namespace rocket6dof
