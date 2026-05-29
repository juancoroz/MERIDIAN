// =====================================================================
//  util_demo.cpp  --  exercises Filer + Table1/2/3 utilities
//
//  Verifies:
//    - Filer tag scoping: same key name appears at top of file with one
//      value and inside the Airframe block with another; we should
//      read the Airframe value when setLine0 is used.
//    - Filer handles inline trailing comments, equals signs, multi-
//      space whitespace, free-form prose lines.
//    - Filer getInt / getString work.
//    - Table1 interpolation: linear, end-clamped.
//    - Table2 bilinear interpolation.
//    - Table3 trilinear interpolation.
//    - All three tables can live in one file.
// =====================================================================

#include "../../osk/osk.h"
#include <cstdio>
#include <cmath>
#include <iostream>

using namespace osk;

static int tests_run = 0, tests_failed = 0;

static void check(const char* what, double got, double want, double tol = 1e-9) {
    ++tests_run;
    bool ok = std::fabs(got - want) <= tol;
    std::printf("  %s %-40s got=%-12g want=%-12g\n",
                ok ? "OK  " : "FAIL", what, got, want);
    if (!ok) ++tests_failed;
}
static void check_eq(const char* what, int got, int want) {
    ++tests_run;
    bool ok = got == want;
    std::printf("  %s %-40s got=%-12d want=%-12d\n",
                ok ? "OK  " : "FAIL", what, got, want);
    if (!ok) ++tests_failed;
}
static void check_str(const char* what, const std::string& got, const std::string& want) {
    ++tests_run;
    bool ok = got == want;
    std::printf("  %s %-40s got='%s' want='%s'\n",
                ok ? "OK  " : "FAIL", what, got.c_str(), want.c_str());
    if (!ok) ++tests_failed;
}

int main() {
    std::printf("== Filer ==\n");
    Filer ff("data.txt");

    // Top-level: xcp should be the value from line 4 (2.45).
    check("xcp at top  ",   ff.getDouble("xcp"), 2.45);
    check("xka at top  ",   ff.getDouble("xka"), -44.486);

    // Scope to Airframe.  Now xcp should be 1.448 even though 2.45
    // appears earlier in the file.
    ff.setLine0("Airframe");
    check("xcp in Airframe",  ff.getDouble("xcp"),  1.448);
    check("xcg in Airframe",  ff.getDouble("xcg"),  1.422);
    check("dia in Airframe",  ff.getDouble("dia"),  0.886315);
    check("cmq with leading 0", ff.getDouble("cmq"), 0.0001);
    check_eq ("n_targets (int)",  ff.getInt("n_targets"), 3);
    check_str("run_title (str)",  ff.getString("run_title"), "baseline");

    std::printf("\n== Table1 (thrust_profile) ==\n");
    Table1 tab1("data.txt");
    tab1.read("thrust_profile");
    // Exact node values
    check("interp at node t=0",   tab1(0.0),   0.0);
    check("interp at node t=2",   tab1(2.0),   2000.0);
    check("interp at node t=4",   tab1(4.0),   4000.0);
    // Linear in the interior: y = 1000*t
    check("interp at t=1.5",      tab1(1.5),   1500.0);
    check("interp at t=3.25",     tab1(3.25),  3250.0);
    // End-clamping
    check("clamp below",          tab1(-1.0),  0.0);
    check("clamp above",          tab1(10.0),  4000.0);

    std::printf("\n== Table2 (cd_table) ==\n");
    Table2 tab2("data.txt");
    tab2.read("cd_table");
    // Spot-check the four corners from the data file:
    //  mach=1, alph=5  -> 5
    //  mach=1, alph=15 -> 15
    //  mach=4, alph=5  -> 20
    //  mach=4, alph=15 -> 60
    check("corner (1,5)",   tab2(1.0,  5.0),  5.0);
    check("corner (1,15)",  tab2(1.0, 15.0), 15.0);
    check("corner (4,5)",   tab2(4.0,  5.0), 20.0);
    check("corner (4,15)",  tab2(4.0, 15.0), 60.0);
    // Interior point: at mach=2.5, alph=7.5
    //  row 2 (mach=2): 10, 20, 30 at alph=5,10,15  -> bilinear in (mach,alph)
    //  By inspection cd = 5*mach*alph/5 = mach*alph; check that pattern:
    //  cd_table[i][j] = mach[i] * alph[j]  -- verify:
    //   (mach=1,alph=5) -> 5 ok; (mach=2,alph=10) -> 20 ok; (mach=4,alph=15) -> 60 ok.
    // So at (2.5, 7.5) the analytical bilinear answer is 2.5 * 7.5 = 18.75.
    check("interior (2.5,7.5)",  tab2(2.5, 7.5),  18.75);
    // Mid-cell: (3.5, 12.5) -> 3.5 * 12.5 = 43.75
    check("interior (3.5,12.5)", tab2(3.5, 12.5), 43.75);
    // Clamping
    check("clamp below both",   tab2(0.0, 0.0),   5.0);
    check("clamp above both",   tab2(9.0, 99.0), 60.0);

    std::printf("\n== Table3 (cd_table_thrust) ==\n");
    Table3 tab3("data.txt");
    tab3.read("cd_table_thrust");
    // Block 0 is the original 2-D table; block 1 is 10x that.
    // Verify corners in block 0:
    check("block0 corner (1,5,1)",   tab3(1.0,  5.0,  1.0),   5.0);
    check("block0 corner (4,15,1)",  tab3(4.0, 15.0,  1.0),  60.0);
    // Verify corners in block 1:
    check("block1 corner (1,5,10)",  tab3(1.0,  5.0, 10.0),  50.0);
    check("block1 corner (4,15,10)", tab3(4.0, 15.0, 10.0), 600.0);
    // Interpolated across thrust: at thrust=5.5 (midpoint of 1..10),
    // result should be the mean of block0 and block1 at same (mach,alph).
    // At (mach=2, alph=10): block0=20, block1=200, midpoint -> 110.
    check("thrust interp (2,10,5.5)", tab3(2.0, 10.0, 5.5), 110.0);
    // Trilinear interior: (2.5, 7.5, 1) should be 2.5*7.5 = 18.75 in block 0.
    check("trilinear at thrust=1",  tab3(2.5, 7.5, 1.0),  18.75);
    // At thrust=10: 10x. So (2.5, 7.5, 10) = 187.5.
    check("trilinear at thrust=10", tab3(2.5, 7.5, 10.0), 187.5);

    std::printf("\n%d tests, %d failed\n", tests_run, tests_failed);
    return tests_failed == 0 ? 0 : 1;
}
