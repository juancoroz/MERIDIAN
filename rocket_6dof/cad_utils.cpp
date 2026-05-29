//  cad_utils.cpp  --  CADAC utility functions, ported to osk types
//
//  Each function below is a line-by-line port of the corresponding
//  function in Zipfel's utility_functions.cpp.  Equation references in
//  comments cite Britting (1971) "Inertial Navigation Systems Analysis"
//  where used in the source.  The only substantive change is the
//  container type: Zipfel's Matrix(3,1) -> osk::Vec, Matrix(3,3) ->
//  osk::Mat.  Element access uses (.x, .y, .z) and operator[][] which
//  are equivalent to Zipfel's get_loc(i,j) / assign_loc(i,j,v).
//
//  Convention reminder (Zipfel naming):
//    I = inertial         D = geodetic (NED at vehicle lat/lon)
//    E = Earth-fixed      G = geographic/geocentric
//    TDI = DCM from I to D (Zipfel writes this as "T.M. of D wrt I")

#include "cad_utils.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace rocket6dof {

namespace {
constexpr double DEG = 180.0 / osk::PI;
constexpr double RAD = osk::PI / 180.0;

inline int sgn(double x) { return (x > 0) - (x < 0); }
} // anon

//  cad_in_geo84:  Inertial position from geodetic (lon, lat, alt).
//  Britting 1971, Eq. 4-15 and 4-21.  WGS84 ellipsoid.
osk::Vec cad_in_geo84(double lon, double lat, double alt, double time) {
    const double F = WGS84_FLATTENING;
    const double A = WGS84_SMAJOR_AXIS;

    // Britting Eq. 4-21: radius of ellipsoid surface at this latitude
    double r0 = A * (1.0 - F * (1.0 - std::cos(2.0*lat)) / 2.0
                         + 5.0 * F*F * (1.0 - std::cos(4.0*lat)) / 16.0);

    // Britting Eq. 4-15: deflection of normal angle
    double dd = F * std::sin(2.0*lat) * (1.0 - F/2.0 - alt/r0);

    double dbi = r0 + alt;

    // Position in geodetic D-frame (the local NED at this lat/lon).
    // In Zipfel's convention this has z-axis pointing inward along
    // the geodetic normal, so a point on the ellipsoid surface at
    // alt=0 has SBID = (-dbi*sin(dd), 0, -dbi*cos(dd)).
    double sbid1 = -dbi * std::sin(dd);
    double sbid3 = -dbi * std::cos(dd);

    // Rotate D -> I via inverse of TDI (which is ~TDI for orthogonal).
    // Zipfel inlines the rotation here as a fixed expression.
    double lon_cel = GW_CLONG + WGS84_WEII3 * time + lon;
    double slat = std::sin(lat);
    double clat = std::cos(lat);
    double slon = std::sin(lon_cel);
    double clon = std::cos(lon_cel);

    double sbii1 = -slat * clon * sbid1 - clat * clon * sbid3;
    double sbii2 = -slat * slon * sbid1 - clat * slon * sbid3;
    double sbii3 =  clat        * sbid1 - slat        * sbid3;

    return osk::Vec(sbii1, sbii2, sbii3);
}

//  cad_geo84_in:  Iterative geodetic lon/lat/alt from inertial pos.
//  Converges in ~3-4 iterations to CAD_SMALL.
void cad_geo84_in(double& lon, double& lat, double& alt,
                  const osk::Vec& SBII, double time) {
    const double F = WGS84_FLATTENING;
    const double A = WGS84_SMAJOR_AXIS;

    double dbi  = SBII.mag();
    double latg = std::asin(SBII.z / dbi);     // geocentric latitude as seed
    lat = latg;

    int count = 0;
    double lat0;
    do {
        lat0 = lat;
        double r0 = A * (1.0 - F * (1.0 - std::cos(2.0*lat0)) / 2.0
                             + 5.0 * F*F * (1.0 - std::cos(4.0*lat0)) / 16.0);
        alt = dbi - r0;
        double dd = F * std::sin(2.0*lat0) * (1.0 - F/2.0 - alt/r0);
        lat = latg + dd;
        if (++count > 100) {
            std::fprintf(stderr,
                         "*** cad_geo84_in: latitude failed to converge ***\n");
            std::exit(1);
        }
    } while (std::fabs(lat - lat0) > CAD_SMALL);

    // Longitude in inertial frame, then subtract Earth rotation
    double dum4 = std::asin(SBII.y / std::sqrt(SBII.x*SBII.x + SBII.y*SBII.y));
    double alamda;
    if      (SBII.x >= 0.0 && SBII.y >= 0.0) alamda = dum4;
    else if (SBII.x <  0.0 && SBII.y >= 0.0) alamda = osk::PI - dum4;
    else if (SBII.x <  0.0 && SBII.y <  0.0) alamda = osk::PI - dum4;
    else                                     alamda = 2.0*osk::PI + dum4;

    lon = alamda - WGS84_WEII3 * time - GW_CLONG;
    if (lon > osk::PI) lon = -(2.0*osk::PI - lon);
}

//  cad_geoc_in:  Spherical-Earth lon/lat/alt from inertial pos.
//  Used internally by cad_grav84.
void cad_geoc_in(double& lonc, double& latc, double& altc,
                 const osk::Vec& SBII, double time) {
    double dbi = SBII.mag();
    latc = std::asin(SBII.z / dbi);
    altc = dbi - REARTH_SPHERICAL;

    double dum4 = std::asin(SBII.y / std::sqrt(SBII.x*SBII.x + SBII.y*SBII.y));
    double lon_cel;
    if      (SBII.x >= 0.0 && SBII.y >= 0.0) lon_cel = dum4;
    else if (SBII.x <  0.0 && SBII.y >= 0.0) lon_cel = osk::PI - dum4;
    else if (SBII.x <  0.0 && SBII.y <  0.0) lon_cel = osk::PI - dum4;
    else                                     lon_cel = 2.0*osk::PI + dum4;

    lonc = lon_cel - WGS84_WEII3 * time - GW_CLONG;
    if (lonc > osk::PI) lonc = -(2.0*osk::PI - lonc);
}

//  cad_grav84:  WGS84 gravity in geocentric frame.
//  GRAVG = (lateral_J2, 0, radial_inv_sq + J2_correction).
osk::Vec cad_grav84(const osk::Vec& SBII, double time) {
    double lonc, latc, altc;
    cad_geoc_in(lonc, latc, altc, SBII, time);

    double dbi  = SBII.mag();
    double dum1 = WGS84_GM / (dbi * dbi);
    double dum2 = 3.0 * std::sqrt(5.0);
    double dum3 = (WGS84_SMAJOR_AXIS / dbi) * (WGS84_SMAJOR_AXIS / dbi);

    double gravg1 = -dum1 * dum2 * WGS84_C20 * dum3 * std::sin(latc)*std::cos(latc);
    double gravg2 =  0.0;
    double gravg3 =  dum1 * (1.0 + dum2/2.0 * WGS84_C20 * dum3
                              * (3.0 * std::sin(latc)*std::sin(latc) - 1.0));

    return osk::Vec(gravg1, gravg2, gravg3);
}

//  cad_tdi84:  DCM from I to D (geodetic NED).
osk::Mat cad_tdi84(double lon, double lat, double /*alt*/, double time) {
    double lon_cel = GW_CLONG + WGS84_WEII3 * time + lon;

    double tdi13 =  std::cos(lat);
    double tdi33 = -std::sin(lat);
    double tdi22 =  std::cos(lon_cel);
    double tdi21 = -std::sin(lon_cel);

    return osk::Mat(
        tdi33*tdi22,   -tdi33*tdi21,   tdi13,
        tdi21,          tdi22,         0.0,
        -tdi13*tdi22,   tdi13*tdi21,   tdi33
    );
}

//  cad_tgi84:  DCM from I to G (geocentric).
//  G is rotated from D by the deflection-of-normal angle 'dd'.
osk::Mat cad_tgi84(double lon, double lat, double alt, double time) {
    const double F = WGS84_FLATTENING;
    const double A = WGS84_SMAJOR_AXIS;

    osk::Mat TDI = cad_tdi84(lon, lat, alt, time);

    double r0 = A * (1.0 - F * (1.0 - std::cos(2.0*lat)) / 2.0
                         + 5.0 * F*F * (1.0 - std::cos(4.0*lat)) / 16.0);
    double dd = F * std::sin(2.0*lat) * (1.0 - F/2.0 - alt/r0);

    // TGD: rotation by 'dd' about the local east axis (y in NED)
    double cdd = std::cos(dd), sdd = std::sin(dd);
    osk::Mat TGD(
         cdd,  0.0,  -sdd,
         0.0,  1.0,   0.0,
         sdd,  0.0,   cdd
    );

    return TGD * TDI;
}

//  cad_tei:  DCM from I to E (Earth-fixed).  Pure rotation about z.
osk::Mat cad_tei(double time) {
    double xi = WGS84_WEII3 * time + GW_CLONG;
    double sxi = std::sin(xi), cxi = std::cos(xi);
    return osk::Mat(
         cxi,  sxi,  0.0,
        -sxi,  cxi,  0.0,
         0.0,  0.0,  1.0
    );
}

//  cad_geo84vel_in:  Geographic velocity magnitude + flight path angles
//  from inertial states.  Uses cad_tdi84 + Earth-rotation correction.
void cad_geo84vel_in(double& dvbe, double& psivdx, double& thtvdx,
                     const osk::Vec& SBII, const osk::Vec& VBII,
                     double time) {
    double lon, lat, alt;
    cad_geo84_in(lon, lat, alt, SBII, time);

    osk::Mat TDI = cad_tdi84(lon, lat, alt, time);

    // VBEI = VBII - WEII x SBII (subtract Earth rotation contribution)
    osk::Vec WEII(0.0, 0.0, WGS84_WEII3);
    osk::Vec VBEI = VBII - WEII.cross(SBII);

    osk::Vec VBED = TDI * VBEI;

    // Convert to polar: magnitude, heading (atan2(E,N)), flight path (-atan(D/horiz))
    dvbe = VBED.mag();
    double vN = VBED.x, vE = VBED.y, vD = VBED.z;
    psivdx = std::atan2(vE, vN) * DEG;
    double vh = std::sqrt(vN*vN + vE*vE);
    thtvdx = std::atan2(-vD, vh) * DEG;
}

//  cad_kepler  --  Project a two-body inertial state forward in time
//
//  Iteratively solves the universal Kepler equation (Bate, Mueller,
//  White, "Fundamentals of Astrodynamics", Dover 1971) to propagate
//  (SBII, VBII) along a Keplerian orbit for tgo seconds.  Used by LTG
//  guidance to predict the inertial state at burnout when only gravity
//  acts between now and then.
//
//  Returns 0 on success, 1 if iteration fails to converge.  On failure
//  the outputs are left unmodified.
//
//  Ported verbatim from Zipfel's utility_functions.cpp:cad_kepler().
int cad_kepler(osk::Vec& SPII, osk::Vec& VPII,
               const osk::Vec& SBII, const osk::Vec& VBII,
               double tgo) {
    const double SMALL = 1.0e-6;
    double sqrt_GM = std::sqrt(WGS84_GM);
    double ro      = SBII.mag();
    double vo      = VBII.mag();
    double rvo     = SBII.dot(VBII);
    double a1      = vo * vo / WGS84_GM;
    double sa      = ro / (2.0 - ro * a1);
    if (sa < 0.0) {
        return 1;     // hyperbolic / parabolic; not handled
    }
    double smua = sqrt_GM * std::sqrt(sa);
    double mdot = smua / (sa * sa);

    double dm = mdot * tgo;
    double de = dm;
    double a11 = rvo / smua;
    double a21 = (sa - ro) / sa;

    double sde = 0.0, cde = 0.0;
    int count = 0;
    double adm = 0.0;
    do {
        cde = 1.0 - std::cos(de);
        sde = std::sin(de);
        double dmn   = de + a11 * cde - a21 * sde;
        double dmerr = dm - dmn;
        adm = std::fabs(dmerr) / mdot;
        double dmde = 1.0 + a11 * sde - a21 * (1.0 - cde);
        de += dmerr / dmde;
        ++count;
        if (count > 20) {
            return 1;
        }
    } while (adm > SMALL);

    // Projected position
    double fk = (ro - sa * cde) / ro;
    double gk = (dm + sde - de) / mdot;
    SPII = SBII * fk + VBII * gk;

    // Projected velocity
    double rp  = SPII.mag();
    double fdk = -smua * sde / ro;
    double gdk = rp - sa * cde;
    VPII = SBII * (fdk / rp) + VBII * (gdk / rp);

    return 0;
}

//  angle  --  Angle between two 3-vectors, [0, pi] radians
//
//  Returns 0 if either vector has zero magnitude.  Equivalent to
//  acos(v1.v2 / (|v1|*|v2|)) with safe clamping.
double angle(const osk::Vec& v1, const osk::Vec& v2) {
    double a1 = v1.mag();
    double a2 = v2.mag();
    if (a1 < 1.0e-12 || a2 < 1.0e-12) return 0.0;
    double c = v1.dot(v2) / (a1 * a2);
    if (c > 1.0)  c = 1.0;
    if (c < -1.0) c = -1.0;
    return std::acos(c);
}

} // namespace rocket6dof
