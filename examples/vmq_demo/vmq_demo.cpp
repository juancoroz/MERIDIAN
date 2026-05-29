// =====================================================================
//  vmq_demo.cpp  --  Vec / Mat / Quat unit tests.
//
//  Verifies every API documented in Appendix 1.2 of the CMD User's
//  Guide, plus the math identities that the 6DOF model code depends
//  on (round-trip Euler -> DCM -> Euler, etc).
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
    std::printf("  %s %-44s got=%-13g want=%-13g\n",
                ok ? "OK  " : "FAIL", what, got, want);
    if (!ok) ++tests_failed;
}

static void check_vec(const char* what, const Vec& got, const Vec& want,
                      double tol = 1e-9) {
    ++tests_run;
    bool ok = std::fabs(got.x - want.x) <= tol &&
              std::fabs(got.y - want.y) <= tol &&
              std::fabs(got.z - want.z) <= tol;
    std::printf("  %s %-44s got=(%g,%g,%g) want=(%g,%g,%g)\n",
                ok ? "OK  " : "FAIL", what,
                got.x, got.y, got.z, want.x, want.y, want.z);
    if (!ok) ++tests_failed;
}

static void check_mat(const char* what, const Mat& got, const Mat& want,
                      double tol = 1e-9) {
    ++tests_run;
    bool ok = true;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            if (std::fabs(got[i][j] - want[i][j]) > tol) ok = false;
    if (ok) {
        std::printf("  OK   %-44s\n", what);
    } else {
        std::printf("  FAIL %s\n    got:  %g %g %g | %g %g %g | %g %g %g\n"
                                "    want: %g %g %g | %g %g %g | %g %g %g\n",
                    what,
                    got[0][0], got[0][1], got[0][2],
                    got[1][0], got[1][1], got[1][2],
                    got[2][0], got[2][1], got[2][2],
                    want[0][0], want[0][1], want[0][2],
                    want[1][0], want[1][1], want[1][2],
                    want[2][0], want[2][1], want[2][2]);
        ++tests_failed;
    }
}

// A small user-supplied transform for testing .apply()
static double add1(double v) { return v + 1.0; }

int main() {
    // ==================== Vec ============================
    std::printf("== Vec ==\n");

    // Constructor + .x/.y/.z access
    Vec v0;
    check("default constructor zero-inits .x", v0.x, 0.0);
    check("default constructor zero-inits .y", v0.y, 0.0);
    check("default constructor zero-inits .z", v0.z, 0.0);

    Vec v1(10., 20., 30.);
    check_vec("Vec(10,20,30) constructor", v1, Vec(10, 20, 30));

    // Indexed access
    Vec v2;
    v2[0] = 100.0; v2[1] = 200.0; v2[2] = 300.0;
    check_vec("indexed write [0..2]", v2, Vec(100, 200, 300));
    check("indexed read v2[1]", v2[1], 200.0);

    // .x/.y/.z access
    Vec v3;
    v3.x = -1.0; v3.y = -2.0; v3.z = -3.0;
    check_vec(".x/.y/.z assignment", v3, Vec(-1, -2, -3));

    // operator() reassignment
    v3(7.0, 8.0, 9.0);
    check_vec("operator() reassignment", v3, Vec(7, 8, 9));

    // Chained reassignment per Appendix 1.2.1.2
    Vec v4 = v3(-100, -200, -300);
    check_vec("chained reassignment in expression", v4, Vec(-100, -200, -300));
    check_vec("source modified by chained reassign", v3, Vec(-100, -200, -300));

    // Scalar broadcast
    Vec v6;
    v6 = 1.0;
    check_vec("scalar broadcast v = 1.0", v6, Vec(1, 1, 1));
    v6 = 0.0;
    check_vec("scalar broadcast v = 0.0", v6, Vec(0, 0, 0));

    // extract()
    Vec v7(1.5, 2.5, 3.5);
    double a, b, c;
    v7.extract(a, b, c);
    check("extract a", a, 1.5);
    check("extract b", b, 2.5);
    check("extract c", c, 3.5);

    // Arithmetic
    Vec va(1, 2, 3), vb(10, 20, 30);
    check_vec("v + v",       va + vb,   Vec(11, 22, 33));
    check_vec("v - v",       vb - va,   Vec(9, 18, 27));
    check_vec("v * scalar",  va * 2.0,  Vec(2, 4, 6));
    check_vec("v / scalar",  vb / 10.0, Vec(1, 2, 3));
    check_vec("scalar * v",  2.0 * va,  Vec(2, 4, 6));

    // In-place scale
    Vec vc(1, 2, 3);
    vc *= 3.0;
    check_vec("v *= scalar (in place)", vc, Vec(3, 6, 9));

    // scale() returns explicit, doesn't modify in place
    Vec vorig(1, 2, 3);
    Vec vscaled = vorig.scale(10.0);
    check_vec("scale() returns scaled",     vscaled, Vec(10, 20, 30));
    check_vec("scale() does NOT modify src", vorig,   Vec(1, 2, 3));

    // mag, .m, unit
    Vec vmag(3, 4, 0);
    check("mag()", vmag.mag(), 5.0);
    vmag.m = vmag.mag();
    check(".m after explicit set", vmag.m, 5.0);
    Vec u = vmag.unit();
    check_vec("unit()", u, Vec(0.6, 0.8, 0.0));
    check("unit() doesn't auto-scale .m", vmag.x, 3.0);  // mag itself unchanged

    // apply
    Vec vapp(-1.0, -2.0, -3.0);
    Vec vabs = vapp.apply(std::fabs);
    check_vec("apply(fabs)", vabs, Vec(1, 2, 3));
    Vec vp1 = vapp.apply(add1);
    check_vec("apply(custom fn)", vp1, Vec(0, -1, -2));

    // dot, cross
    Vec u1(1, 0, 0), u2(0, 1, 0), u3(0, 0, 1);
    check("dot(u1, u2) = 0",  u1.dot(u2), 0.0);
    check("dot(u1, u1) = 1",  u1.dot(u1), 1.0);
    check_vec("cross(u1, u2) = u3", u1.cross(u2), u3);
    check_vec("cross(u2, u3) = u1", u2.cross(u3), u1);
    check_vec("cross(u3, u1) = u2", u3.cross(u1), u2);
    // Anti-commutativity
    check_vec("cross(u2, u1) = -u3", u2.cross(u1), Vec(0, 0, -1));

    // ==================== Mat ============================
    std::printf("\n== Mat ==\n");

    // Constructors
    Mat m0;
    check_mat("default constructor zero-inits", m0,
              Mat(0,0,0, 0,0,0, 0,0,0));
    Mat m3(10., 20., 30., 40., 50., 60., 70., 80., 90.);
    check_mat("Mat(9 doubles row-major)", m3,
              Mat(10,20,30, 40,50,60, 70,80,90));
    Mat m4(Vec(1, 2, 3), Vec(4, 5, 6), Vec(7, 8, 9));
    check_mat("Mat(v1, v2, v3) row constructor", m4,
              Mat(1,2,3, 4,5,6, 7,8,9));

    // Element access
    check("m[0][1]", m3[0][1], 20.0);
    check("m[2][2]", m3[2][2], 90.0);
    check_vec("m[0] returns row 0",    m4[0], Vec(1, 2, 3));
    check_vec("m[1] returns row 1",    m4[1], Vec(4, 5, 6));

    // Reassignment
    Mat m5;
    m5(1, 2, 3, 4, 5, 6, 7, 8, 9);
    check_mat("m(...) reassignment", m5, m4);
    m5(Vec(9,8,7), Vec(6,5,4), Vec(3,2,1));
    check_mat("m(v1,v2,v3) reassignment", m5,
              Mat(9,8,7, 6,5,4, 3,2,1));

    // Matrix arithmetic
    Mat I(1,0,0, 0,1,0, 0,0,1);   // identity
    check_mat("I + I",
              I + I,
              Mat(2,0,0, 0,2,0, 0,0,2));
    check_mat("m - m = 0",
              m4 - m4,
              Mat(0,0,0, 0,0,0, 0,0,0));

    // Matrix * vector
    Vec result_v = m4 * Vec(1, 1, 1);
    check_vec("m * (1,1,1) = sum of each row", result_v,
              Vec(6, 15, 24));
    check_vec("I * v = v",
              I * Vec(7, 8, 9),
              Vec(7, 8, 9));

    // Matrix * matrix
    check_mat("I * m = m", I * m4, m4);
    check_mat("m * I = m", m4 * I, m4);

    // Scalar scaling
    check_mat("m * 2",  m4 * 2.0,
              Mat(2,4,6, 8,10,12, 14,16,18));
    check_mat("m / 2", Mat(2,4,6, 8,10,12, 14,16,18) / 2.0, m4);
    check_mat("2.0 * m", 2.0 * m4,
              Mat(2,4,6, 8,10,12, 14,16,18));

    Mat mm(1,2,3, 4,5,6, 7,8,9);
    Mat mm_orig = mm;
    Mat mm_scaled = mm.scale(0.5);
    check_mat("scale() returns scaled",     mm_scaled,
              Mat(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5));
    check_mat("scale() does NOT modify src", mm, mm_orig);

    // Transpose
    Mat mt(1, 2, 3, 4, 5, 6, 7, 8, 9);
    check_mat("transpose()", mt.transpose(),
              Mat(1, 4, 7, 2, 5, 8, 3, 6, 9));

    // Determinant
    check("det(I) = 1", I.det(), 1.0);
    check("det(singular) = 0", Mat(1,2,3, 4,5,6, 7,8,9).det(), 0.0);
    check("det(diag(2,3,4)) = 24",
          Mat(2,0,0, 0,3,0, 0,0,4).det(), 24.0);
    // det of rotation matrix should be 1
    Vec euler(0.3, -0.4, 0.7);
    Mat R = euler.getDCM();
    check("det(rotation matrix) = 1", R.det(), 1.0, 1e-10);

    // Inverse
    Mat A(2, 0, 0,  0, 3, 0,  0, 0, 4);
    Mat Ainv = A.inv();
    check_mat("inv(diag(2,3,4)) = diag(1/2,1/3,1/4)",
              Ainv,
              Mat(0.5, 0, 0,  0, 1.0/3.0, 0,  0, 0, 0.25));
    check_mat("A * A.inv() = I", A * Ainv, I, 1e-10);
    // Rotation: inverse is transpose
    Mat Rinv = R.inv();
    check_mat("inv(R) == transpose(R) for rotation", Rinv,
              R.transpose(), 1e-10);

    // apply
    Mat mapp(-1, -2, -3, -4, -5, -6, -7, -8, -9);
    check_mat("apply(fabs)", mapp.apply(std::fabs),
              Mat(1, 2, 3, 4, 5, 6, 7, 8, 9));

    // ==================== Quat ============================
    std::printf("\n== Quat ==\n");

    Quat q0;
    check("Quat default .s", q0.s, 0.0);
    check("Quat default .x", q0.x, 0.0);
    check("Quat default .y", q0.y, 0.0);
    check("Quat default .z", q0.z, 0.0);

    Quat q1(0.962, -0.023, 0.084, 0.258);
    check("Quat .s init", q1.s, 0.962);
    check("Quat .x init", q1.x, -0.023);
    check("Quat .y init", q1.y, 0.084);
    check("Quat .z init", q1.z, 0.258);

    // Indexed access
    check("q[0]", q1[0], 0.962);
    check("q[3]", q1[3], 0.258);

    // Reassignment
    q1(1.0, 0.0, 0.0, 0.0);
    check("Quat reassign .s", q1.s, 1.0);
    check("Quat reassign .x", q1.x, 0.0);

    // Normalize
    Quat qn(2.0, 0.0, 0.0, 0.0);
    qn.normalize();
    check("normalize() in-place affects .s", qn.s, 1.0);
    Quat qn2(1.0, 1.0, 1.0, 1.0);
    qn2.normalize();
    check("normalize() yields unit quaternion", qn2.mag(), 1.0, 1e-12);

    // ==================== Round-trip identity ============
    std::printf("\n== Euler <-> DCM <-> Quat round-trips ==\n");

    // Test several Euler triples, including some near edge cases
    double triples[][3] = {
        { 0.0, 0.0, 0.0 },
        { 0.1, 0.2, 0.3 },
        { -0.5, 0.4, -1.0 },
        { 0.3, -0.4, 0.7 },
        { 1.5, 0.6, -0.8 }
    };
    for (auto& t : triples) {
        Vec  e(t[0], t[1], t[2]);
        Mat  D = e.getDCM();
        Vec  e_back = D.getEuler();
        char buf[80];
        std::snprintf(buf, sizeof buf, "Vec(%.2f,%.2f,%.2f).getDCM().getEuler()",
                      t[0], t[1], t[2]);
        check_vec(buf, e_back, e, 1e-10);

        Quat Q = e.getQuat();
        Mat  D2 = Q.getDCM();
        std::snprintf(buf, sizeof buf, "Vec(%.2f,%.2f,%.2f).getQuat().getDCM()",
                      t[0], t[1], t[2]);
        check_mat(buf, D2, D, 1e-10);

        Quat Q2 = D.getQuat();
        Mat  D3 = Q2.getDCM();
        std::snprintf(buf, sizeof buf, "Mat(...).getQuat().getDCM() matches");
        check_mat(buf, D3, D, 1e-10);
    }

    // ==================== Physical sanity ================
    std::printf("\n== Rotation physics ==\n");

    // 90-degree yaw (psi = pi/2) about +z: the inertial-x vector,
    // expressed in body frame after the yaw rotation, becomes -y_body
    // (since the body now points along the original +y direction, so
    // the old +x is now to the body's right, which is -y in body).
    Vec yaw90(0.0, 0.0, PI / 2.0);
    Mat Dyaw = yaw90.getDCM();
    Vec xi(1, 0, 0);
    Vec xb = Dyaw * xi;
    check_vec("DCM(yaw=90) applied to inertial-x", xb,
              Vec(0.0, -1.0, 0.0), 1e-10);

    // 90-degree pitch (theta = pi/2) about +y in NED convention
    // (+z down): the inertial-x vector becomes +z in body (nose
    // points up, so the original "forward" direction is now "up" in
    // body, and with +z down, "up" = +z_body.
    // Let's just verify what the formula gives.)  With theta=pi/2,
    // phi=psi=0, DCM = [0 0 -1; 0 1 0; 1 0 0].  So DCM * (1,0,0) =
    // (0, 0, 1).
    Vec pitch90(0.0, PI / 2.0, 0.0);
    Mat Dpitch = pitch90.getDCM();
    Vec xb2 = Dpitch * xi;
    check_vec("DCM(pitch=90) applied to inertial-x",
              xb2, Vec(0.0, 0.0, 1.0), 1e-10);

    // Quaternion of zero rotation = (1, 0, 0, 0)
    Vec zero_euler(0, 0, 0);
    Quat qid = zero_euler.getQuat();
    check_vec("Quat(0,0,0).x,y,z = 0", Vec(qid.x, qid.y, qid.z),
              Vec(0, 0, 0), 1e-15);
    check("Quat(0,0,0).s = 1", qid.s, 1.0, 1e-15);

    // Identity DCM should produce zero Euler angles
    Vec id_euler = I.getEuler();
    check_vec("I.getEuler() = (0,0,0)", id_euler, Vec(0,0,0), 1e-15);

    std::printf("\n%d tests, %d failed\n", tests_run, tests_failed);
    return tests_failed == 0 ? 0 : 1;
}
