//  euler.cpp  --  Rotational EOM for the 6DOF rocket
//
//  Cross-reference to Zipfel "Modeling and Simulation of Aerospace
//  Vehicle Dynamics" 3rd ed., Section 10.3.x (round6[150-199]):
//
//      Zipfel name      OSK member        Meaning
//      -----------      ----------        -------
//      WBIB             WBIB              ang vel body wrt I, body coords  [rad/s]
//      WBIBD            WBIBD             angular acceleration             [rad/s^2]
//      WBEB             WBEB              ang vel body wrt Earth, body co. [rad/s]
//      WBII             WBII              ang vel body wrt I, inertial co. [rad/s]
//      ppx, qqx, rrx    ppx, qqx, rrx     roll/pitch/yaw rates wrt Earth   [deg/s]
//
//  Equation (Zipfel Eq. 4.7, also called Euler's rotational EOM):
//
//      I_BBB * d/dt(WBIB) = FMB - WBIB x (I_BBB * WBIB)
//
//  Solved for the derivative:
//
//      WBIBD = inv(I_BBB) * (FMB - skew(WBIB) * I_BBB * WBIB)
//
//  where skew(w) is the 3x3 skew-symmetric matrix such that
//  skew(w) * v = w cross v.
//
//  Note on initial conditions:  Zipfel's user input is ppx0/qqx0/rrx0,
//  which are body rates WBEB wrt the rotating Earth frame.  Init
//  converts these to WBIB (inertial-relative) via:
//      WBIB = WBEB + TBI * WEII
//  This means Kinematics must be initialized BEFORE Euler so that
//  TBI(0) is available.  In mearth=0 (non-rotating) mode, WEII=0 and
//  WBIB = WBEB directly.

#include "euler.h"
#include "forces.h"
#include "mass_props.h"
#include "kinematics.h"
#include "propulsion.h"
#include "cad_utils.h"
#include <cmath>
#include <cstdio>

namespace rocket6dof {

namespace {
constexpr double DEG = 180.0 / osk::PI;
constexpr double RAD = osk::PI / 180.0;

// Skew-symmetric matrix [w]_x such that [w]_x * v = w x v.
osk::Mat skew(const osk::Vec& w) {
    return osk::Mat(
         0.0, -w.z,  w.y,
         w.z,  0.0, -w.x,
        -w.y,  w.x,  0.0
    );
}
} // anon

Euler::Euler()
    : forces(nullptr), mass(nullptr), kin(nullptr), prop(nullptr)
{
    mearth = 0;
    ppx0 = qqx0 = rrx0 = 0.0;

    wx = wy = wz = 0.0;
    wxd = wyd = wzd = 0.0;

    WBIB  = osk::Vec(0, 0, 0);
    WBIBD = osk::Vec(0, 0, 0);
    WBEB  = osk::Vec(0, 0, 0);
    WBII  = osk::Vec(0, 0, 0);
    ppx = qqx = rrx = 0.0;

    // Register the three angular-velocity integrator states
    addIntegrator(wx, wxd);
    addIntegrator(wy, wyd);
    addIntegrator(wz, wzd);
}

void Euler::init() {
    if (initCount == 0) {
        // User-supplied body rates wrt Earth, converted to rad/s
        WBEB = osk::Vec(ppx0 * RAD, qqx0 * RAD, rrx0 * RAD);

        // Convert WBEB -> WBIB by adding the Earth-rotation contribution
        // expressed in body frame: WBIB = WBEB + TBI * WEII.
        if (mearth == 1 && kin) {
            osk::Vec WEII(0.0, 0.0, WGS84_WEII3);
            WBIB = WBEB + kin->TBI * WEII;
        } else {
            // Non-rotating Earth: WEII = 0, so WBIB == WBEB
            WBIB = WBEB;
        }

        wx = WBIB.x;  wy = WBIB.y;  wz = WBIB.z;
    }
    // On subsequent stage entries, keep WBIB as it was -- multi-stage
    // continuity pattern.
}

void Euler::update() {
    // ---- Pull integrator scalars into WBIB ----
    WBIB = osk::Vec(wx, wy, wz);

    // ---- Read inputs from upstream blocks (with safe defaults) ----
    osk::Vec FMB(0, 0, 0);
    osk::Mat IBBB(1,0,0, 0,1,0, 0,0,1);   // identity default
    osk::Mat TBI (1,0,0, 0,1,0, 0,0,1);   // identity default

    if (forces) FMB  = forces->FMB_();
    if (prop)        IBBB = prop->IBBB_();
    else if (mass)   IBBB = mass->IBBB_();
    if (kin)    TBI  = kin->TBI_();

    // ---- Euler's rotational EOM ----
    //   WBIBD = inv(IBBB) * (FMB - skew(WBIB) * IBBB * WBIB)
    osk::Vec gyro = skew(WBIB) * (IBBB * WBIB);     // = WBIB x (IBBB * WBIB)
    osk::Vec rhs  = FMB - gyro;
    WBIBD = IBBB.inv() * rhs;

    // ---- Push derivatives to the integrator scalars ----
    wxd = WBIBD.x;
    wyd = WBIBD.y;
    wzd = WBIBD.z;

    // ---- Diagnostics ----
    // Angular velocity wrt inertial frame in inertial coordinates
    WBII = TBI.transpose() * WBIB;

    // Angular velocity wrt Earth in body coordinates
    if (mearth == 1) {
        osk::Vec WEII(0.0, 0.0, WGS84_WEII3);
        WBEB = WBIB - TBI * WEII;
    } else {
        WBEB = WBIB;
    }

    // Body rates in deg/s
    ppx = WBEB.x * DEG;
    qqx = WBEB.y * DEG;
    rrx = WBEB.z * DEG;
}

void Euler::rpt() {
    if (osk::State::sample(1.0)) {
        std::printf("Eulr t=%7.3f  ppx=%8.3f  qqx=%8.3f  rrx=%8.3f  |WBIB|=%8.4f rad/s\n",
                    osk::State::t, ppx, qqx, rrx, WBIB.mag());
    }
}

} // namespace rocket6dof
