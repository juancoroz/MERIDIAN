//  aero_deck_test.cpp  --  Verify the aero deck loader against
//  Zipfel's aero_deck_SLV.asc and our own aero.txt

#include "aero_deck.h"

#include <cmath>
#include <cstdio>
#include <string>

using namespace rocket6dof;

namespace {

int fails = 0;

bool check(bool cond, const char* desc) {
    std::printf("  %s %s\n", cond ? "OK  " : "FAIL", desc);
    if (!cond) fails++;
    return cond;
}

bool close(double a, double b, double tol) {
    return std::fabs(a - b) <= tol;
}

// ---- Test 1: OSK format - load and interpolate ----
int test_osk_1d() {
    std::printf("\n=== Test 1: OSK 1D format (aero.txt ca0slv3_vs_mach) ===\n");

    AeroTable1 t("aero.txt");
    bool ok = t.read("ca0slv3_vs_mach");
    check(ok, "read OK");
    check(t.size() == 12, "size == 12 entries");

    // Expected values from aero.txt:
    //  0.05 -> 0.2295
    //  0.95 -> 0.7780
    //  1.20 -> 1.0428  (peak)
    //  25.00 -> 0.3636
    check(close(t.interp(0.05), 0.2295, 1e-4), "interp at 0.05 == 0.2295");
    check(close(t.interp(0.95), 0.7780, 1e-4), "interp at 0.95 == 0.7780");
    check(close(t.interp(1.20), 1.0428, 1e-4), "interp at 1.20 == 1.0428");
    check(close(t.interp(25.0), 0.3636, 1e-4), "interp at 25 == 0.3636");

    // Linear interpolation at M=0.5: between (0.05, 0.2295) and (0.6, 0.2475).
    // t = (0.5 - 0.05) / (0.6 - 0.05) = 0.45/0.55 = 0.81818
    // y = 0.2295 + 0.8182*(0.2475 - 0.2295) = 0.2295 + 0.0147 = 0.2442
    check(close(t.interp(0.5), 0.2442, 1e-3), "linear interp at 0.5 ~ 0.2442");

    // Extrapolation should clamp
    check(close(t.interp(-10.0), 0.2295, 1e-4),  "below range -> first value");
    check(close(t.interp(100.0), 0.3636, 1e-4),  "above range -> last value");

    return 0;
}

// ---- Test 2: OSK format - 2D ----
int test_osk_2d() {
    std::printf("\n=== Test 2: OSK 2D format (aero.txt cn0slv3_vs_mach_alpha) ===\n");

    AeroTable2 t("aero.txt");
    bool ok = t.read("cn0slv3_vs_mach_alpha");
    check(ok, "read OK");
    check(t.n1() == 12, "n1 (Mach) == 12");
    check(t.n2() == 6, "n2 (alpha) == 6");

    // Corner: M=0.05, alpha=0 -> 0.0000
    check(close(t.interp(0.05, 0.0), 0.0000, 1e-4), "(0.05, 0) == 0.0000");

    // At M=0.05, alpha=10: from our aero.txt the last column should be
    // monotonically growing.  Don't pin to exact value (depends on our
    // truncated table) -- just verify it's positive and reasonable.
    double v = t.interp(0.05, 10.0);
    check(v > 0.3 && v < 1.0, "(0.05, 10) is reasonable normal force");

    return 0;
}

// ---- Test 3: DATCOM format - 1D ----
int test_datcom_1d() {
    std::printf("\n=== Test 3: DATCOM 1D format (aero_deck_SLV.asc) ===\n");

    // Use the Zipfel reference deck if available
    AeroTable1 t("../../ssr_unpack/aero_deck_SLV.asc");
    bool ok = t.read("ca0slv3_vs_mach");
    check(ok, "read ca0slv3_vs_mach from .asc");
    if (!ok) return 0;

    check(t.size() == 17, "size == 17 entries (full Zipfel deck)");

    // Zipfel values:
    //   0.05 -> 0.2295
    //   0.60 -> 0.2475
    //   1.20 -> 1.0428
    //  25.00 -> 0.3636
    check(close(t.interp(0.05), 0.2295, 1e-4), "M=0.05 -> 0.2295");
    check(close(t.interp(0.60), 0.2475, 1e-4), "M=0.60 -> 0.2475");
    check(close(t.interp(1.20), 1.0428, 1e-4), "M=1.20 -> 1.0428");
    check(close(t.interp(25.0), 0.3636, 1e-4), "M=25  -> 0.3636");

    // 2nd order verification: M=0.30 -> 0.2080 (this point is in the
    // .asc deck but not in our truncated aero.txt)
    check(close(t.interp(0.30), 0.2080, 1e-4), "M=0.30 -> 0.2080 (Zipfel-only point)");

    return 0;
}

// ---- Test 4: DATCOM format - 2D ----
int test_datcom_2d() {
    std::printf("\n=== Test 4: DATCOM 2D format (aero_deck_SLV.asc) ===\n");

    AeroTable2 t("../../ssr_unpack/aero_deck_SLV.asc");
    bool ok = t.read("cn0slv3_vs_mach_alpha");
    check(ok, "read cn0slv3_vs_mach_alpha from .asc");
    if (!ok) return 0;

    check(t.n1() == 17, "n1 (Mach) == 17");
    check(t.n2() == 11, "n2 (alpha) == 11");

    // From Zipfel's deck, row M=0.05 starts: 0.05 0.00 0.0000 0.0715 0.1546 ...
    // So z(M=0.05, alpha=0) = 0.0000 and z(M=0.05, alpha=2) = 0.0715
    check(close(t.interp(0.05, 0.0),  0.0000, 1e-4), "(M=0.05, a=0)  == 0.0000");
    check(close(t.interp(0.05, 2.0),  0.0715, 1e-4), "(M=0.05, a=2)  == 0.0715");
    check(close(t.interp(0.05, 4.0),  0.1546, 1e-4), "(M=0.05, a=4)  == 0.1546");
    check(close(t.interp(0.05, 30.0), 2.0214, 1e-4), "(M=0.05, a=30) == 2.0214");

    // Row M=1.05 alpha=2:  0.0659
    check(close(t.interp(1.05, 2.0),  0.0659, 1e-4), "(M=1.05, a=2)  == 0.0659");

    // Row M=25.0 last (alpha=30): 3.5823
    check(close(t.interp(25.0, 30.0), 3.5823, 1e-4), "(M=25,   a=30) == 3.5823");

    return 0;
}

// ---- Test 5: DATCOM - read damping table later in same file ----
int test_datcom_seek() {
    std::printf("\n=== Test 5: DATCOM format - tags later in file ===\n");

    AeroTable1 t("../../ssr_unpack/aero_deck_SLV.asc");
    bool ok = t.read("clmqslv3_vs_mach");
    check(ok, "read clmqslv3_vs_mach (5th table in the file)");
    if (!ok) return 0;

    check(t.size() == 17, "size == 17");
    // From Zipfel:  0.05 -> -0.2197;  1.20 -> -0.6539
    check(close(t.interp(0.05), -0.2197, 1e-4), "M=0.05 -> -0.2197");
    check(close(t.interp(1.20), -0.6539, 1e-4), "M=1.20 -> -0.6539");

    return 0;
}

// ---- Test 6: format auto-detection ----
int test_format_detection() {
    std::printf("\n=== Test 6: format auto-detection ===\n");

    check(detect_aero_format("aero.txt")          == AERO_FORMAT_OSK,
          "aero.txt -> OSK");
    check(detect_aero_format("foo.dat")           == AERO_FORMAT_OSK,
          "foo.dat -> OSK");
    check(detect_aero_format("aero_deck_SLV.asc") == AERO_FORMAT_DATCOM,
          "aero_deck_SLV.asc -> DATCOM");
    check(detect_aero_format("FOO.ASC")           == AERO_FORMAT_DATCOM,
          "FOO.ASC -> DATCOM (case-insensitive)");

    return 0;
}

} // anon

int main() {
    test_format_detection();
    test_osk_1d();
    test_osk_2d();
    test_datcom_1d();
    test_datcom_2d();
    test_datcom_seek();
    std::printf("\n=== Total failures: %d ===\n", fails);
    return fails == 0 ? 0 : 1;
}
