//  cad_test.cpp  --  Standalone verification of the CAD utility functions.
//
//  Three tests:
//    (1) Round-trip lon/lat/alt -> SBII -> lon/lat/alt should preserve
//        all three to ~1e-10 (geodetic latitude convergence threshold).
//    (2) Gravity at the equator at WGS84 surface should be ~9.7803 m/s^2
//        (the published value); at the pole should be ~9.8322 m/s^2.
//    (3) TDI should be orthogonal (||I - T*T^T|| < 1e-14).  TGI too.
//        Earth rotation: TEI at t=0 should be identity; at t=86164.09 s
//        (1 sidereal day) should be identity again.

#include "cad_utils.h"
#include <cmath>
#include <cstdio>

using namespace rocket6dof;

namespace {
constexpr double DEG = 180.0 / osk::PI;
constexpr double RAD = osk::PI / 180.0;

double mat_inf_norm(const osk::Mat& M) {
    double m = 0.0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            m = std::max(m, std::fabs(M[i][j]));
    return m;
}

int test_roundtrip() {
    std::printf("\n=== Test 1: lon/lat/alt round-trip ===\n");
    struct Pt { double lon_deg, lat_deg, alt_m; const char* name; };
    Pt pts[] = {
        {   0.0,   0.0,      0.0, "Equator, prime meridian, sea level" },
        {   0.0,  90.0,      0.0, "North pole, sea level" },
        { 139.7,  35.7,    100.0, "Tokyo, 100 m alt" },
        { -74.0,  40.7,  10000.0, "New York, 10 km alt" },
        {  45.0,  45.0, 500000.0, "Mid-lat, 500 km alt" },
    };

    int fails = 0;
    for (const Pt& p : pts) {
        double lon = p.lon_deg * RAD;
        double lat = p.lat_deg * RAD;
        double alt = p.alt_m;
        double t   = 0.0;

        osk::Vec SBII = cad_in_geo84(lon, lat, alt, t);

        double lon2, lat2, alt2;
        cad_geo84_in(lon2, lat2, alt2, SBII, t);

        double dlon = (lon2 - lon) * DEG;
        double dlat = (lat2 - lat) * DEG;
        double dalt = alt2 - alt;

        bool ok = std::fabs(dlon) < 1e-5
               && std::fabs(dlat) < 1e-5
               && std::fabs(dalt) < 1e-3;

        std::printf("  %-40s  dlon=%+.2e  dlat=%+.2e  dalt=%+.2e  %s\n",
                    p.name, dlon, dlat, dalt, ok ? "OK" : "FAIL");
        if (!ok) ++fails;
    }
    return fails;
}

int test_gravity() {
    std::printf("\n=== Test 2: gravity at canonical points ===\n");
    int fails = 0;

    // Equator, sea level: published value ~9.7803 m/s^2.  The inverse-
    // square component is GM/a^2 = 3.986e14/6378137^2 = 9.79828 m/s^2.
    // With C20 oblateness correction we should land near 9.8144 m/s^2 in
    // the geocentric frame, which is what 'grav' magnitude reports here.
    // (Surface gravity also includes centrifugal subtraction, which we
    // do NOT include since this function returns true gravitation only.)
    {
        osk::Vec S = cad_in_geo84(0.0, 0.0, 0.0, 0.0);
        osk::Vec G = cad_grav84(S, 0.0);
        double mag = G.mag();
        double expected = 9.8144;   // GM/a^2 + J2 correction at equator
        bool ok = std::fabs(mag - expected) < 0.05;
        std::printf("  Equator surface gravity   = %.4f m/s^2 (expect ~%.4f)  %s\n",
                    mag, expected, ok ? "OK" : "FAIL");
        if (!ok) ++fails;
    }
    {
        osk::Vec S = cad_in_geo84(0.0, 90.0*RAD, 0.0, 0.0);
        osk::Vec G = cad_grav84(S, 0.0);
        double mag = G.mag();
        double expected = 9.8322;   // pole gravity, published
        bool ok = std::fabs(mag - expected) < 0.05;
        std::printf("  Pole surface gravity      = %.4f m/s^2 (expect ~%.4f)  %s\n",
                    mag, expected, ok ? "OK" : "FAIL");
        if (!ok) ++fails;
    }
    return fails;
}

int test_dcms() {
    std::printf("\n=== Test 3: DCM orthogonality and Earth rotation ===\n");
    int fails = 0;
    osk::Mat I3(1,0,0, 0,1,0, 0,0,1);

    // Test orthogonality of TDI at a mid-latitude
    {
        osk::Mat T = cad_tdi84(45.0*RAD, 30.0*RAD, 1000.0, 0.0);
        osk::Mat E = I3 - T * T.transpose();
        double n = mat_inf_norm(E);
        bool ok = n < 1e-14;
        std::printf("  ||I - TDI*TDI^T||_inf       = %.2e   %s\n",
                    n, ok ? "OK" : "FAIL");
        if (!ok) ++fails;
    }
    // Test orthogonality of TGI
    {
        osk::Mat T = cad_tgi84(45.0*RAD, 30.0*RAD, 1000.0, 0.0);
        osk::Mat E = I3 - T * T.transpose();
        double n = mat_inf_norm(E);
        bool ok = n < 1e-14;
        std::printf("  ||I - TGI*TGI^T||_inf       = %.2e   %s\n",
                    n, ok ? "OK" : "FAIL");
        if (!ok) ++fails;
    }
    // Test TEI(0) = I
    {
        osk::Mat T = cad_tei(0.0);
        osk::Mat E = I3 - T;
        double n = mat_inf_norm(E);
        bool ok = n < 1e-14;
        std::printf("  ||TEI(0) - I||_inf          = %.2e   %s\n",
                    n, ok ? "OK" : "FAIL");
        if (!ok) ++fails;
    }
    // Test TEI(sidereal day) = I
    {
        double t_sidereal = 2.0 * osk::PI / WGS84_WEII3;   // exact
        osk::Mat T = cad_tei(t_sidereal);
        osk::Mat E = I3 - T;
        double n = mat_inf_norm(E);
        bool ok = n < 1e-10;
        std::printf("  ||TEI(siderealday) - I||_inf = %.2e   %s\n",
                    n, ok ? "OK" : "FAIL");
        if (!ok) ++fails;
    }
    return fails;
}

} // anon

int main() {
    int fails = 0;
    fails += test_roundtrip();
    fails += test_gravity();
    fails += test_dcms();
    std::printf("\n=== Total failures: %d ===\n", fails);
    return fails == 0 ? 0 : 1;
}
