//  kinematics.cpp  --  Rotational kinematics for the 6DOF rocket
//
//  Cross-reference to Zipfel "Modeling and Simulation of Aerospace
//  Vehicle Dynamics" 3rd ed., Section 10.3.x (round6[100-149]):
//
//      Zipfel name      OSK member       Meaning
//      -----------      ----------       -------
//      TBI              TBI              body-from-inertial DCM
//      TBD              TBD              body-from-geodetic DCM
//      ortho_error      ortho_error      |I - T*T^T| residual
//      psibd,thtbd,     psibd,thtbd,     Euler angles (3-2-1) [rad]
//        phibd            phibd
//      psibdx,thtbdx,   psibdx,...       same [deg]
//        phibdx
//      alphax, betax    alphax, betax    incidence angles [deg]
//      alppx, phipx     alppx, phipx     aeroballistic angles [deg]
//
//  Key design choices, with rationale:
//
//  (1) DCM integration vs quaternion: Zipfel integrates TBI as a 3x3
//      matrix and orthonormalizes via Bar-Itzhack each step.  We do
//      the same -- 9 scalar states + orthonormalize -- to keep the
//      book cross-reference exact.  Quaternion integration would be
//      4 states + normalize, and OSK supports it via the Quat utility,
//      but porting Zipfel's downstream blocks would then need a
//      quaternion-to-DCM step every update.  Not worth the complexity.
//
//  (2) Orthonormalization placement: the correction must apply to the
//      POST-integration TBI, not the pre-integration value.  In OSK,
//      integration happens in propagate() AFTER update().  So the
//      orthonormalization runs at the start of update(), gated on
//      State::stepstart && !State::tickfirst -- meaning "a full step
//      has just finished, and we're at the boundary".  This runs once
//      per integration step, not once per RK4 sub-stage.
//
//  (3) Non-rotating Earth: TDI = identity, so TBD = TBI.  Zipfel's
//      code has TBD = TBI * ~TDI; the simplification falls out
//      automatically when TDI = I.

#include "kinematics.h"
#include "environment.h"
#include "newton.h"
#include "euler.h"
#include "cad_utils.h"
#include <cmath>
#include <cstdio>

namespace rocket6dof {

namespace {
constexpr double DEG = 180.0 / osk::PI;
constexpr double RAD = osk::PI / 180.0;
constexpr double EPS = 1.0e-10;

// Build a 3-2-1 (yaw-pitch-roll) DCM: body from reference.
// Same convention as Zipfel's mat3tr(psi, theta, phi).
osk::Mat mat3tr(double psi, double theta, double phi) {
    double cpsi = std::cos(psi),   spsi = std::sin(psi);
    double cthe = std::cos(theta), sthe = std::sin(theta);
    double cphi = std::cos(phi),   sphi = std::sin(phi);
    return osk::Mat(
        cthe*cpsi,                     cthe*spsi,                    -sthe,
        sphi*sthe*cpsi - cphi*spsi,    sphi*sthe*spsi + cphi*cpsi,    sphi*cthe,
        cphi*sthe*cpsi + sphi*spsi,    cphi*sthe*spsi - sphi*cpsi,    cphi*cthe
    );
}

// Skew-symmetric matrix [w]_x such that [w]_x * v = w cross v.
osk::Mat skew(const osk::Vec& w) {
    return osk::Mat(
         0.0, -w.z,  w.y,
         w.z,  0.0, -w.x,
        -w.y,  w.x,  0.0
    );
}

inline double sign(double x) { return (x > 0) - (x < 0); }
} // anon

Kinematics::Kinematics()
    : env(nullptr), newton(nullptr), euler(nullptr)
{
    mearth = 0;
    psibdx0 = thtbdx0 = phibdx0 = 0.0;

    // Initialise TBI to identity in raw scalars
    t00 = 1; t01 = 0; t02 = 0;
    t10 = 0; t11 = 1; t12 = 0;
    t20 = 0; t21 = 0; t22 = 1;
    t00d = t01d = t02d = 0.0;
    t10d = t11d = t12d = 0.0;
    t20d = t21d = t22d = 0.0;

    ortho_error = 0.0;
    psibd = thtbd = phibd = 0.0;
    psibdx = thtbdx = phibdx = 0.0;
    alphax = betax = 0.0;
    alppx = phipx = 0.0;
    alphaix = betaix = 0.0;

    // Register all 9 DCM scalars
    addIntegrator(t00, t00d);
    addIntegrator(t01, t01d);
    addIntegrator(t02, t02d);
    addIntegrator(t10, t10d);
    addIntegrator(t11, t11d);
    addIntegrator(t12, t12d);
    addIntegrator(t20, t20d);
    addIntegrator(t21, t21d);
    addIntegrator(t22, t22d);
}

void Kinematics::init() {
    if (initCount == 0) {
        // Initial body-from-geodetic DCM from Euler angles (3-2-1 sequence)
        osk::Mat TBD0 = mat3tr(psibdx0*RAD, thtbdx0*RAD, phibdx0*RAD);

        // TBI(0) = TBD(0) * TDI(0).  Even in spherical, non-rotating
        // mode, TDI is generally not identity (it depends on lat/lon
        // at the launch site).  We use cad_tdi84 to compute it in BOTH
        // modes; in mearth=0 the time argument is held at 0 so TDI
        // remains constant for the run.
        osk::Mat TDI0(1,0,0, 0,1,0, 0,0,1);
        if (newton) {
            double lon = newton->lonx0 * RAD;
            double lat = newton->latx0 * RAD;
            double alt = newton->alt0;
            // In mearth=1 mode, time is the actual sim time (0 at init);
            // in mearth=0 mode, 0 is passed so Earth rotation is frozen.
            double t0 = 0.0;
            TDI0 = cad_tdi84(lon, lat, alt, t0);
        }
        osk::Mat T0 = TBD0 * TDI0;

        t00 = T0[0][0]; t01 = T0[0][1]; t02 = T0[0][2];
        t10 = T0[1][0]; t11 = T0[1][1]; t12 = T0[1][2];
        t20 = T0[2][0]; t21 = T0[2][1]; t22 = T0[2][2];
        TBI = T0;
        TBD = TBD0;
    }
}

void Kinematics::update() {
    // ---- (1) Pull the integrator scalars into TBI ----
    TBI = osk::Mat(t00, t01, t02,
                   t10, t11, t12,
                   t20, t21, t22);

    // ---- (2) Orthonormalize ONCE PER STEP at the step boundary ----
    // Bar-Itzhack first-order: T <- T + 0.5*(I - T*T^T)*T
    // This must not run inside RK4 sub-stages, only at step boundaries
    // after a full integration step completes.
    if (osk::State::stepstart && !osk::State::tickfirst) {
        osk::Mat I_(1,0,0, 0,1,0, 0,0,1);
        osk::Mat E  = I_ - TBI * TBI.transpose();
        TBI = TBI + (E * TBI) * 0.5;

        // Track residual for diagnostics (sum of |diag(E)|^2 sqrt)
        double e1 = E[0][0], e2 = E[1][1], e3 = E[2][2];
        ortho_error = std::sqrt(e1*e1 + e2*e2 + e3*e3);

        // Write back to scalars so the integrator continues from the
        // orthonormalized state, not the drifted one.
        t00 = TBI[0][0]; t01 = TBI[0][1]; t02 = TBI[0][2];
        t10 = TBI[1][0]; t11 = TBI[1][1]; t12 = TBI[1][2];
        t20 = TBI[2][0]; t21 = TBI[2][1]; t22 = TBI[2][2];
    }

    // ---- (3) Compute DCM derivative from angular velocity ----
    osk::Vec WBIB(0, 0, 0);
    if (euler) WBIB = euler->WBIB_();

    // d/dt TBI = -[WBIB]_x * TBI
    osk::Mat TBID = skew(WBIB) * TBI * (-1.0);
    t00d = TBID[0][0]; t01d = TBID[0][1]; t02d = TBID[0][2];
    t10d = TBID[1][0]; t11d = TBID[1][1]; t12d = TBID[1][2];
    t20d = TBID[2][0]; t21d = TBID[2][1]; t22d = TBID[2][2];

    // ---- (4) TBD = TBI * ~TDI ----
    // In mearth=1 mode, TDI varies with time (Earth rotation).  In
    // mearth=0 mode, TDI is held at its t=0 value (non-rotating Earth).
    if (newton) {
        double t_for_tdi = (mearth == 1) ? osk::State::t : 0.0;
        osk::Mat TDI = cad_tdi84(newton->lonx * RAD, newton->latx * RAD,
                                 newton->alt,        t_for_tdi);
        TBD = TBI * TDI.transpose();
    } else {
        TBD = TBI;
    }

    // ---- (5) Euler angle extraction (3-2-1) ----
    double tbd13 = TBD[0][2];
    double tbd11 = TBD[0][0];
    double tbd33 = TBD[2][2];
    double tbd12 = TBD[0][1];
    double tbd23 = TBD[1][2];

    double cthtbd;
    if (std::fabs(tbd13) < 1.0) {
        thtbd  = std::asin(-tbd13);
        cthtbd = std::cos(thtbd);
    } else {
        thtbd  = (osk::PI/2.0) * sign(-tbd13);
        cthtbd = EPS;
    }

    double cpsi = tbd11 / cthtbd;
    if (std::fabs(cpsi) > 1.0) cpsi = sign(cpsi);
    psibd = std::acos(cpsi) * sign(tbd12);

    double cphi = tbd33 / cthtbd;
    if (std::fabs(cphi) > 1.0) cphi = sign(cphi);
    phibd = std::acos(cphi) * sign(tbd23);

    psibdx = psibd * DEG;
    thtbdx = thtbd * DEG;
    phibdx = phibd * DEG;

    // ---- (6) Incidence angles using wind vector VAED ----
    if (newton && env) {
        osk::Vec VBED = newton->VBED_();
        osk::Vec VAED = env->VAED_();
        // Recompute dvba locally from VBED-VAED rather than reading
        // env->dvba.  Env runs AFTER Kine in our block order, so
        // env->dvba would be the previous substage's value, leading to
        // numerical issues when VBAB.x marginally exceeds the stale
        // dvba and the acos clamp pushes alpp to 0 spuriously.
        osk::Vec VBAD = VBED - VAED;
        double   dvba = VBAD.mag();
        osk::Vec VBAB = TBD * VBAD;

        double alpha = std::atan2(VBAB.z, VBAB.x);
        double beta  = (dvba > EPS) ? std::asin(VBAB.y / dvba) : 0.0;
        alphax = alpha * DEG;
        betax  = beta  * DEG;

        // Aeroballistic angles (total alpha, aero roll)
        double dum = (dvba > EPS) ? VBAB.x / dvba : 1.0;
        if (std::fabs(dum) > 1.0) dum = sign(dum);
        double alpp = std::acos(dum);

        double phip;
        if (VBAB.y == 0.0 && VBAB.z == 0.0) {
            phip = 0.0;
        } else if (std::fabs(VBAB.y) < EPS) {
            phip = (VBAB.z >= 0) ? 0.0 : osk::PI;
        } else {
            phip = std::atan2(VBAB.y, VBAB.z);
        }
        alppx = alpp * DEG;
        phipx = phip * DEG;

        // Inertial incidence diagnostics
        osk::Vec VBII = newton->VBII_();
        osk::Vec VBIB = TBI * VBII;
        double   dvbi = VBII.mag();
        alphaix = std::atan2(VBIB.z, VBIB.x) * DEG;
        betaix  = (dvbi > EPS) ? std::asin(VBIB.y / dvbi) * DEG : 0.0;
    }
}

void Kinematics::rpt() {
    if (osk::State::sample(1.0)) {
        std::printf("Kine t=%7.3f  psi=%8.3f  tht=%8.3f  phi=%8.3f  "
                    "alpha=%6.2f  beta=%6.2f  ortho=%.2e\n",
                    osk::State::t, psibdx, thtbdx, phibdx,
                    alphax, betax, ortho_error);
    }
}

} // namespace rocket6dof
