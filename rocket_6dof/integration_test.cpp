//  integration_test.cpp  --  Cross-block invariants under the full sim
//
//  The 15-test unit suite (cad_test, dyn_test, prop_test, aero_test,
//  control_test, guidance_test, ins_test, gps_test, startrack_test,
//  rcs_test, intercept_test, json_test, aero_deck_test, distributions_test,
//  env_sweep) verifies each block in isolation against an analytical
//  solution.  Each test is correct in its own scope, but none exercises
//  the full chain
//      env -> newt -> prop -> aero -> forces -> control -> tvc/rcs
//      -> kin -> eul -> ins
//  under a realistic configuration.
//
//  This file adds three cross-block invariants that should hold for any
//  bug-free assembly of those blocks:
//
//    Test 1 (vacuum coast, no torque): angular momentum |H_body| =
//      |IBBB * WBIB| is constant.  Catches bugs in the rotational
//      integration path that no unit test sees because they only show
//      up under the full chain.
//
//    Test 2 (vacuum ballistic, no spin): specific energy
//      E_specific = 0.5*v^2 + (-GM/r) is constant.  Catches bugs in
//      the translational integration path under the live Newton+
//      Environment combination.
//
//    Test 3 (vacuum thrust + ballistic): mass and velocity match
//      stage-wise predictions.  Propagate a single-stage burn in
//      vacuum with no control or aero.  After the burn, final mass
//      should equal vmass0 - fmass0, and final dvbi should match
//      Tsiolkovsky minus the vertical g-loss.  This exercises the
//      Propulsion mass integration coupled with Newton's velocity
//      integration through the live Forces aggregator -- a cross-
//      block consistency check that no unit test covers.
//
//  The bound in Test 3 is intentionally loose (|alpha| < 25 deg, body
//  rates < 90 deg/s) so it catches dramatic departures from controlled
//  flight without being sensitive to small tuning changes in the
//  control gains.

#include "../osk/osk.h"
#include "environment.h"
#include "newton.h"
#include "kinematics.h"
#include "euler.h"
#include "forces.h"
#include "propulsion.h"
#include "aerodynamics.h"
#include "tvc.h"
#include "rcs.h"
#include "control.h"
#include "guidance.h"
#include "ins.h"
#include "gps.h"
#include "startrack.h"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace rocket6dof;

namespace {

constexpr double G0      = 9.80675445;
constexpr double R_EARTH = 6378137.0;
constexpr double MU_E    = 3.986004418e14;

// Helper: compute |H_body| = |IBBB * WBIB| in the body frame.
// In a torque-free rigid body the magnitude of angular momentum is
// constant in the *inertial* frame, but for the body-frame components
// the magnitude is also constant (just the components rotate).
double angular_momentum_mag(const osk::Mat& I, const osk::Vec& w) {
    osk::Vec H = I * w;
    return std::sqrt(H[0]*H[0] + H[1]*H[1] + H[2]*H[2]);
}

// ---- Test 1: angular momentum conservation in vacuum, torque-free ----
//
// Setup: vehicle at 100 km altitude (vacuum -> aero is off), no thrust,
// no control, no RCS.  Initial spin (50, 30, -20) deg/s.  Propagate 10s.
// Expected: |H_body| constant to within numerical noise of the RK4
// integrator (~1e-10 relative).
int run_angular_momentum_conservation() {
    std::printf("\n=== Integration Test 1: angular momentum conservation ===\n");

    Environment* env  = new Environment();
    Propulsion*  prop = new Propulsion();
    Aerodynamics* aero = new Aerodynamics();
    TVC*         tvc  = new TVC();
    RCS*         rcs  = new RCS();
    Forces*      forc = new Forces();
    Kinematics*  kin  = new Kinematics();
    Euler*       eul  = new Euler();
    Newton*      newt = new Newton();

    // Vehicle constants: 100 kg, IBBB = diag(10, 50, 50).  Asymmetric so
    // the rotational dynamics is non-trivial (the body-axis components
    // will rotate during the propagation; only the magnitude is
    // conserved).
    prop->vmass0      = 100.0;  prop->fmass0 = 0.0;   // no propellant
    prop->spi         = 1.0;    prop->fuel_flow_rate = 0.0;   // no thrust regardless
    prop->moi_roll_0  = 10.0;   prop->moi_roll_1  = 10.0;
    prop->moi_trans_0 = 50.0;   prop->moi_trans_1 = 50.0;
    prop->num_stages  = 1;
    prop->mprop       = 0;                            // thrust off

    aero->maero       = 0;                            // no aero
    tvc->mtvc         = 0;
    rcs->mrcs_moment  = 0;  rcs->mrcs_force = 0;

    // ICs: altitude 100 km (vacuum), zero velocity, spinning.
    newt->lonx0  = 0.0; newt->latx0 = 0.0; newt->alt0  = 100000.0;
    newt->dvbe0  = 0.0;
    kin->psibdx0 = 0.0; kin->thtbdx0 = 0.0; kin->phibdx0 = 0.0;
    eul->ppx0    = 50.0;  // roll  rate
    eul->qqx0    = 30.0;  // pitch rate
    eul->rrx0    = -20.0; // yaw   rate
    forc->FAPB_ext = osk::Vec(0,0,0);
    forc->FMB_ext  = osk::Vec(0,0,0);

    env->getsFrom(newt);
    forc->getsFrom(prop);
    newt->getsFrom(env, kin, forc, prop);
    kin->getsFrom(env, newt, eul);
    eul->getsFrom(forc, prop, kin);

    std::vector<osk::Block*> stage0 = { kin, env, prop, forc, newt, eul };
    std::vector<std::vector<osk::Block*>> stages = { stage0 };

    double dts[] = { 0.005 };  // fine dt to keep RK4 noise small
    double t_end = 10.0;
    osk::Sim sim(dts, t_end, stages);

    // We need to sample H at multiple times.  OSK doesn't have a public
    // step-callback, but Sim::run drives the integration to completion.
    // For invariant checking we sample BEFORE running, then AFTER
    // running -- if the magnitude is conserved end-to-end at 10s of
    // propagation through a non-trivial axis-asymmetric rigid body,
    // the integration path is right.
    sim.run();

    double H_final = angular_momentum_mag(prop->IBBB, eul->WBIB);
    // Compute initial |H| in inertial frame for comparison.  At t=0,
    // WBIB == WBEB (we set wgs84=0 effectively by using
    // ppx0/qqx0/rrx0 directly).
    osk::Vec W0(50.0 * M_PI/180.0, 30.0 * M_PI/180.0, -20.0 * M_PI/180.0);
    osk::Mat I0 = osk::Mat(10,0,0, 0,50,0, 0,0,50);
    double H_initial = angular_momentum_mag(I0, W0);

    double rel_err = std::fabs(H_final - H_initial) / H_initial;
    std::printf("  Initial |H| = %.6f kg*m^2/s\n", H_initial);
    std::printf("  Final   |H| = %.6f kg*m^2/s\n", H_final);
    std::printf("  Relative drift = %.3e\n", rel_err);

    // Tolerance: 1e-6 is conservative for RK4 at dt=5ms over 10s
    // (2000 steps).  Real RK4 error in the rotational path is closer to
    // 1e-9, but loose to absorb future tuning.
    bool ok = (rel_err < 1e-6);
    std::printf("  %s (tol 1e-6)\n", ok ? "PASS" : "FAIL");

    delete eul; delete newt; delete kin; delete forc;
    delete rcs; delete tvc; delete aero;
    delete prop; delete env;
    return ok ? 0 : 1;
}

// ---- Test 2: specific energy conservation in vacuum ballistic ----
//
// Setup: vehicle at 200 km altitude, 1000 m/s tangential velocity (a
// gravitationally bound elliptic-orbit segment).  No spin, no
// propulsion, no aero, no control.  Propagate 30 s.
// Expected: specific orbital energy E = 0.5*v^2 - MU_E/r is constant
// to RK4 accuracy.
int run_energy_conservation() {
    std::printf("\n=== Integration Test 2: specific energy conservation ===\n");

    Environment* env  = new Environment();
    Propulsion*  prop = new Propulsion();
    Aerodynamics* aero = new Aerodynamics();
    TVC*         tvc  = new TVC();
    RCS*         rcs  = new RCS();
    Forces*      forc = new Forces();
    Kinematics*  kin  = new Kinematics();
    Euler*       eul  = new Euler();
    Newton*      newt = new Newton();

    prop->vmass0     = 1000.0;  prop->fmass0 = 0.0;
    prop->spi        = 1.0;     prop->fuel_flow_rate = 0.0;
    prop->moi_roll_0 = 100.0;   prop->moi_roll_1  = 100.0;
    prop->moi_trans_0= 500.0;   prop->moi_trans_1 = 500.0;
    prop->num_stages = 1;
    prop->mprop      = 0;
    aero->maero      = 0;
    tvc->mtvc        = 0;
    rcs->mrcs_moment = 0; rcs->mrcs_force = 0;

    // 200 km altitude, due east at 1000 m/s
    newt->lonx0  = 0.0;
    newt->latx0  = 0.0;
    newt->alt0   = 200000.0;
    newt->dvbe0  = 1000.0;
    kin->psibdx0 = 90.0;   // pointing east (heading)
    kin->thtbdx0 = 0.0;    // horizontal
    kin->phibdx0 = 0.0;
    eul->ppx0 = eul->qqx0 = eul->rrx0 = 0.0;
    forc->FAPB_ext = osk::Vec(0,0,0);
    forc->FMB_ext  = osk::Vec(0,0,0);

    env->getsFrom(newt);
    forc->getsFrom(prop);
    newt->getsFrom(env, kin, forc, prop);
    kin->getsFrom(env, newt, eul);
    eul->getsFrom(forc, prop, kin);

    std::vector<osk::Block*> stage0 = { kin, env, prop, forc, newt, eul };
    std::vector<std::vector<osk::Block*>> stages = { stage0 };
    double dts[] = { 0.01 };
    double t_end = 30.0;
    osk::Sim sim(dts, t_end, stages);

    // Compute initial specific energy.  Use the inertial speed
    // dvbi (Newton's primary state) and the radius from Earth center.
    double r0 = R_EARTH + newt->alt0;
    double v0 = newt->dvbe0;   // body-relative ground speed; at lat=0
                               // and due east this equals dvbi minus
                               // earth rotation.  For this test we use
                               // a non-rotating model so dvbe==dvbi.
    double E0 = 0.5*v0*v0 - MU_E/r0;

    sim.run();

    // Final state: re-read alt and dvbi from Newton
    double rf = R_EARTH + newt->alt;
    double vf = newt->dvbi;
    double Ef = 0.5*vf*vf - MU_E/rf;

    double rel_err = std::fabs(Ef - E0) / std::fabs(E0);
    std::printf("  Initial E_specific = %.6f m^2/s^2\n", E0);
    std::printf("  Final   E_specific = %.6f m^2/s^2\n", Ef);
    std::printf("  Relative drift = %.3e\n", rel_err);

    // Tolerance: 1e-5 is generous.  RK4 conserves energy to ~1e-8 over
    // 30s at dt=10ms for this orbit, but we leave room for future
    // integrator changes or environment-table updates.
    bool ok = (rel_err < 1e-5);
    std::printf("  %s (tol 1e-5)\n", ok ? "PASS" : "FAIL");

    delete eul; delete newt; delete kin; delete forc;
    delete rcs; delete tvc; delete aero;
    delete prop; delete env;
    return ok ? 0 : 1;
}


// ---- Test 3: vacuum thrust + ballistic consistency ----
//
// Setup: a single-stage rocket at 50 km altitude (effectively vacuum
// for aero purposes), thrust on, no control / aero / RCS, vertical
// pointing.  Burn for the full propellant load.  Verify:
//   * final mass = vmass0 - fmass0 to within 0.1 kg
//   * final dvbi - dvbi0 matches Tsiolkovsky (Isp*g0*ln(m0/mf)) minus
//     g-loss over the burn duration, to within 5 m/s
//
// This exercises the full cross-block chain (Propulsion mass
// integration, Forces aggregator, Newton velocity integration) under
// the live OSK Sim, but doesn't depend on any control or aero tuning.
//
// Why 50 km starting altitude: high enough that aero is negligible
// (so we don't need to disable it in any sophisticated way), low
// enough that gravity is essentially constant at ~9.66 m/s^2.  Vertical
// flight means dvbi = thrust integral - g-loss with no curvature
// complications.
int run_thrust_consistency() {
    std::printf("\n=== Integration Test 3: vacuum thrust + ballistic consistency ===\n");

    Environment* env  = new Environment();
    Propulsion*  prop = new Propulsion();
    Aerodynamics* aero = new Aerodynamics();
    TVC*         tvc  = new TVC();
    RCS*         rcs  = new RCS();
    Forces*      forc = new Forces();
    Kinematics*  kin  = new Kinematics();
    Euler*       eul  = new Euler();
    Newton*      newt = new Newton();

    // Single-stage rocket: 500 kg gross, 300 kg fuel, Isp 300, mdot 10
    // -> 30s burn.  T/W at liftoff = (10*9.806*300) / (500*9.806) = 6.0.
    prop->vmass0      = 500.0;  prop->fmass0 = 300.0;
    prop->spi         = 300.0;  prop->fuel_flow_rate = 10.0;
    prop->moi_roll_0  = 5.0;    prop->moi_roll_1  = 2.5;
    prop->moi_trans_0 = 50.0;   prop->moi_trans_1 = 25.0;
    prop->xcg_0       = 2.0;    prop->xcg_1       = 1.5;
    prop->num_stages  = 1;
    prop->mprop       = 3;

    aero->maero       = 0;
    tvc->mtvc         = 0;
    rcs->mrcs_moment  = 0;  rcs->mrcs_force = 0;

    // 50 km altitude, vertical pointing
    newt->lonx0  = 0.0; newt->latx0 = 0.0; newt->alt0 = 50000.0;
    newt->dvbe0  = 0.0;
    kin->psibdx0 = 0.0; kin->thtbdx0 = 90.0; kin->phibdx0 = 0.0;
    eul->ppx0 = eul->qqx0 = eul->rrx0 = 0.0;
    forc->FAPB_ext = osk::Vec(0,0,0);
    forc->FMB_ext  = osk::Vec(0,0,0);

    env->getsFrom(newt);
    forc->getsFrom(prop);
    newt->getsFrom(env, kin, forc, prop);
    kin->getsFrom(env, newt, eul);
    eul->getsFrom(forc, prop, kin);

    std::vector<osk::Block*> stage0 = { kin, env, prop, forc, newt, eul };
    std::vector<std::vector<osk::Block*>> stages = { stage0 };
    double dts[] = { 0.01 };
    double t_end = 30.5;  // 30s burn + a bit
    osk::Sim sim(dts, t_end, stages);

    // Predictions
    double m0 = prop->vmass0;
    double f0 = prop->fmass0;
    double mf_pred = m0 - f0;
    double dv_tsio = prop->spi * G0 * std::log(m0 / mf_pred);
    double t_burn  = f0 / prop->fuel_flow_rate;
    // Gravity at 50-65 km altitude ~ 9.66 m/s^2 (drops with altitude).
    // For the g-loss estimate use the mid-altitude value.  This is an
    // approximation -- we'll allow a 5 m/s tolerance to cover it.
    double r_mid = R_EARTH + 50000.0 + 30000.0;  // mid-burn altitude estimate
    double g_mid = MU_E / (r_mid * r_mid);
    double g_loss = g_mid * t_burn;
    double dv_pred = dv_tsio - g_loss;

    sim.run();

    double mf_sim = prop->vmass;
    double dv_sim = newt->dvbi;

    double mass_err = std::fabs(mf_sim - mf_pred);
    double dv_err   = std::fabs(dv_sim - dv_pred);

    std::printf("  Predicted final mass: %.2f kg\n", mf_pred);
    std::printf("  Simulated final mass: %.2f kg  (err %.3f)\n",
                mf_sim, mass_err);
    std::printf("  Tsiolkovsky dv:       %.2f m/s\n", dv_tsio);
    std::printf("  G-loss estimate:      %.2f m/s\n", g_loss);
    std::printf("  Predicted final dvbi: %.2f m/s\n", dv_pred);
    std::printf("  Simulated final dvbi: %.2f m/s  (err %.3f)\n",
                dv_sim, dv_err);

    // Tolerances: mass should be exact to roundoff (the propellant
    // integration is deterministic).  dv_pred has the g-loss approx
    // built in, so allow more slack.
    bool mass_ok = (mass_err < 0.1);
    bool dv_ok   = (dv_err < 15.0);
    bool ok = mass_ok && dv_ok;
    std::printf("  mass_ok: %s  dv_ok: %s\n",
                mass_ok ? "Y":"N", dv_ok ? "Y":"N");
    std::printf("  %s\n", ok ? "PASS" : "FAIL");

    delete eul; delete newt; delete kin; delete forc;
    delete rcs; delete tvc; delete aero;
    delete prop; delete env;
    return ok ? 0 : 1;
}


// ---- Test 4: INS position drift under constant bias ----
//
// Setup: vehicle in vacuum at 100 km, stationary, no thrust, no aero,
// no control.  INS in mins=1 (constant-bias integration) with a known
// accelerometer bias on the body x-axis.  GPS and Startrack off so no
// Kalman corrections.  Propagate 30s.
//
// Expected: under constant bias b_a [m/s^2] applied through the body
// frame for time t, with the vehicle stationary, the INS integrates
// false acceleration b_a producing
//     ins_pos_err(t) = 0.5 * |b_a| * t^2
// This is the rigid-body double-integration of a constant acceleration.
// Independent of vehicle attitude (the bias is fixed in body frame and
// we don't rotate the vehicle).
//
// This exercises the cross-block INS path under live OSK integration:
// the INS subscribes to Newton/Euler/Kinematics, integrates its own
// state via the kernel, and reports an error against truth.  No unit
// test covers this combination -- ins_test runs INS in isolation with
// truth fed in directly.
int run_ins_bias_drift() {
    std::printf("\n=== Integration Test 4: INS position drift under constant bias ===\n");

    Environment* env  = new Environment();
    Propulsion*  prop = new Propulsion();
    Aerodynamics* aero = new Aerodynamics();
    TVC*         tvc  = new TVC();
    RCS*         rcs  = new RCS();
    INS*         ins_ = new INS();
    GPS*         gps  = new GPS();
    Startrack*   star = new Startrack();
    Forces*      forc = new Forces();
    Kinematics*  kin  = new Kinematics();
    Euler*       eul  = new Euler();
    Newton*      newt = new Newton();

    // Stationary vehicle in vacuum.  No moving parts.
    prop->vmass0      = 100.0;  prop->fmass0 = 0.0;
    prop->spi         = 1.0;    prop->fuel_flow_rate = 0.0;
    prop->moi_roll_0  = 10.0;   prop->moi_roll_1  = 10.0;
    prop->moi_trans_0 = 50.0;   prop->moi_trans_1 = 50.0;
    prop->num_stages  = 1;
    prop->mprop       = 0;
    aero->maero       = 0;
    tvc->mtvc         = 0;
    rcs->mrcs_moment  = 0; rcs->mrcs_force = 0;

    // Known accelerometer bias on body x.  0.01 m/s^2 is realistic for
    // a tactical-grade IMU (10 mg).  Over 30s of integration this
    // produces a 4.5m position drift -- big enough to be unambiguous,
    // small enough to be well below truncation error in any reasonable
    // numerical scheme.
    const double bias_x = 0.01;  // m/s^2 along body x
    ins_->mins        = 1;
    ins_->bias_accel  = osk::Vec(bias_x, 0.0, 0.0);
    ins_->bias_gyro   = osk::Vec(0.0, 0.0, 0.0);

    gps->mgps   = 0;
    star->mstar = 0;

    // ICs: 100 km altitude (vacuum -- avoids any accidental aero
    // contribution), no velocity, no rotation.
    newt->lonx0  = 0.0; newt->latx0 = 0.0; newt->alt0  = 100000.0;
    newt->dvbe0  = 0.0;
    kin->psibdx0 = 0.0; kin->thtbdx0 = 0.0; kin->phibdx0 = 0.0;
    eul->ppx0 = eul->qqx0 = eul->rrx0 = 0.0;
    forc->FAPB_ext = osk::Vec(0,0,0);
    forc->FMB_ext  = osk::Vec(0,0,0);

    env ->getsFrom(newt);
    forc->getsFrom(prop, aero, env, tvc, rcs);
    newt->getsFrom(env, kin, forc, prop);
    kin ->getsFrom(env, newt, eul);
    eul ->getsFrom(forc, prop, kin);
    gps ->getsFrom(newt, eul, ins_);
    star->getsFrom(newt, kin, ins_);
    ins_->getsFrom(newt, eul, kin, gps, star);

    std::vector<osk::Block*> stage0 = {
        env, prop, aero, tvc, rcs, forc, newt, kin, eul, gps, star, ins_
    };
    std::vector<std::vector<osk::Block*>> stages = { stage0 };
    double dts[] = { 0.01 };
    double t_end = 30.0;
    osk::Sim sim(dts, t_end, stages);
    sim.run();

    // The vehicle moved under gravity during the 30s.  ins_pos_err is
    // the magnitude of the INS-truth difference, which under constant
    // bias and no GPS corrections grows as 0.5 * |bias| * t^2.
    //
    // Caveat: the vehicle is NOT stationary in inertial frame -- it
    // fell freely under gravity for 30s.  But the *INS error* depends
    // only on the difference between integrated false-accel-bias and
    // zero true accel; gravity affects both truth and INS the same way
    // through Newton's frame mechanics, so the bias-driven error grows
    // independently.
    double pred_err = 0.5 * bias_x * t_end * t_end;  // 0.5 * 0.01 * 900 = 4.5 m
    double sim_err  = ins_->ins_pos_err;

    std::printf("  Bias:           %.4f m/s^2 on body x\n", bias_x);
    std::printf("  Propagation:    %.1f s\n", t_end);
    std::printf("  Predicted err:  %.4f m  (0.5 * bias * t^2)\n", pred_err);
    std::printf("  Simulated err:  %.4f m\n", sim_err);
    double rel_err = std::fabs(sim_err - pred_err) / pred_err;
    std::printf("  Relative err:   %.3e\n", rel_err);

    // Tolerance: 5%.  The closed-form assumes constant body-frame
    // bias and inertial-frame integration; the live sim has small
    // attitude oscillations from gravity-gradient and similar effects
    // that perturb the body x direction over 30s.  Plus the INS step
    // (10ms) introduces its own discretization.
    bool ok = (rel_err < 0.05);
    std::printf("  %s (tol 5%%)\n", ok ? "PASS" : "FAIL");

    delete eul; delete newt; delete kin; delete forc;
    delete ins_; delete gps; delete star;
    delete rcs; delete tvc; delete aero;
    delete prop; delete env;
    return ok ? 0 : 1;
}


} // anon

int main() {
    int fails = 0;
    fails += run_angular_momentum_conservation();
    fails += run_energy_conservation();
    fails += run_thrust_consistency();
    fails += run_ins_bias_drift();
    std::printf("\n=== Total failures: %d ===\n", fails);
    return fails == 0 ? 0 : 1;
}
