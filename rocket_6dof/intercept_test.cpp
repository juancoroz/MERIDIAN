//  intercept_test.cpp  --  Verify Intercept fires on each condition
//
//  Three tests:
//    (1) Ground impact: drop a body from 1 km altitude; verify Intercept
//        fires at the right time (analytic free-fall: t = sqrt(2h/g)).
//    (2) Apogee: launch a vertical rocket with limited burn; verify
//        Intercept fires at apogee, with vertical velocity ~ 0.
//    (3) Max time: set t_max to a value smaller than the natural run;
//        verify the sim terminates at that time.
//
//  After each test, verify Sim::stop was set (the sim terminated early)
//  and trcond holds the correct condition code.

#include "../osk/osk.h"
#include "environment.h"
#include "newton.h"
#include "kinematics.h"
#include "euler.h"
#include "forces.h"
#include "propulsion.h"
#include "intercept.h"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace rocket6dof;

namespace {

constexpr double G0 = 9.80675445;

int run_ground_impact() {
    std::printf("\n=== Test 1: Ground impact detection ===\n");

    Environment* env   = new Environment();
    Propulsion*  prop  = new Propulsion();
    Forces*      forc  = new Forces();
    Euler*       eul   = new Euler();
    Kinematics*  kin   = new Kinematics();
    Newton*      newt  = new Newton();
    Intercept*   icpt  = new Intercept();

    // Inert vehicle: no thrust, no aero, just gravity
    prop->mprop  = 0;
    prop->vmass0 = 100.0;
    prop->fmass0 = 0.0;
    forc->FAPB_ext = osk::Vec(0, 0, 0);
    forc->FMB_ext  = osk::Vec(0, 0, 0);
    eul->ppx0 = eul->qqx0 = eul->rrx0 = 0.0;

    // Drop from 1 km altitude with zero velocity
    newt->lonx0 = 0.0; newt->latx0 = 0.0; newt->alt0 = 1000.0;
    newt->dvbe0 = 0.0;
    kin->psibdx0 = 0.0; kin->thtbdx0 = 0.0; kin->phibdx0 = 0.0;

    // Intercept: only ground-impact check enabled
    icpt->check_ground_impact = true;
    icpt->check_apogee        = false;
    icpt->check_time_max      = false;

    forc->getsFrom(prop);
    newt->getsFrom(env, kin, forc, prop);
    kin->getsFrom(env, newt, eul);
    eul->getsFrom(forc, prop, kin);
    env->getsFrom(newt);
    icpt->getsFrom(env, newt);

    // Stage order: kin, env, prop, forc, newt, eul, icpt
    std::vector<osk::Block*> stage0 = { kin, env, prop, forc, newt, eul, icpt };
    std::vector<std::vector<osk::Block*>> stages = { stage0 };

    // Run for 30 sec max -- should terminate well before that
    double dts[] = { 0.01 };
    osk::Sim sim(dts, 30.0, stages);
    sim.run();

    // Expected: free fall from h=1000m takes t = sqrt(2*h/g).
    // The vehicle starts at the equator at h=1000 m; gravity weakens
    // slightly with altitude.  Use the average gravity over the fall.
    const double MU = 3.986004418e14;
    const double Re = 6378137.0;
    double g_avg     = MU / ((Re + 500.0) * (Re + 500.0));
    double t_predict = std::sqrt(2.0 * newt->alt0 / g_avg);

    std::printf("  Predicted impact time (free fall): t = %.3f s\n", t_predict);
    std::printf("  Sim terminated at t = %.3f s\n", icpt->t_terminate);
    std::printf("  trcond = %d (1 = ground impact)\n", icpt->trcond);
    std::printf("  final alt = %.2f m\n", newt->alt);

    bool ok =    icpt->trcond == 1
              && std::fabs(icpt->t_terminate - t_predict) < 0.05
              && newt->alt <= 0.1
              && newt->alt > -10.0;     // didn't overshoot far below ground

    std::printf("  %s\n", ok ? "PASS" : "FAIL");

    delete icpt; delete newt; delete kin; delete eul; delete forc; delete prop; delete env;
    return ok ? 0 : 1;
}

int run_apogee() {
    std::printf("\n=== Test 2: Apogee detection ===\n");

    Environment* env   = new Environment();
    Propulsion*  prop  = new Propulsion();
    Forces*      forc  = new Forces();
    Euler*       eul   = new Euler();
    Kinematics*  kin   = new Kinematics();
    Newton*      newt  = new Newton();
    Intercept*   icpt  = new Intercept();

    // Brief thrust to give the rocket some altitude, then coast to apogee.
    prop->mprop          = 3;
    prop->vmass0         = 1000.0;
    prop->fmass0         = 100.0;     // small fuel = short burn
    prop->spi            = 300.0;
    prop->fuel_flow_rate = 50.0;      // 2 s burn
    prop->moi_roll_0  = 10.0;  prop->moi_roll_1  = 10.0;
    prop->moi_trans_0 = 100.0; prop->moi_trans_1 = 100.0;

    forc->FAPB_ext = osk::Vec(0, 0, 0);
    forc->FMB_ext  = osk::Vec(0, 0, 0);
    eul->ppx0 = eul->qqx0 = eul->rrx0 = 0.0;

    // Vertical launch (body pointing up: theta = 90)
    newt->lonx0 = 0.0; newt->latx0 = 0.0; newt->alt0 = 0.0;
    newt->dvbe0 = 0.0;
    kin->psibdx0 = 0.0; kin->thtbdx0 = 90.0; kin->phibdx0 = 0.0;

    // Only apogee detection
    icpt->check_ground_impact = false;
    icpt->check_apogee        = true;
    icpt->check_time_max      = false;
    icpt->t_no_check_until    = 0.5;   // wait until rocket has clearly lifted

    forc->getsFrom(prop);
    newt->getsFrom(env, kin, forc, prop);
    kin->getsFrom(env, newt, eul);
    eul->getsFrom(forc, prop, kin);
    env->getsFrom(newt);
    icpt->getsFrom(env, newt);

    std::vector<osk::Block*> stage0 = { kin, env, prop, forc, newt, eul, icpt };
    std::vector<std::vector<osk::Block*>> stages = { stage0 };

    double dts[] = { 0.01 };
    osk::Sim sim(dts, 200.0, stages);
    sim.run();

    // Analytic prediction (constant-g approximation, gravity ~9.80):
    //   t_burn = 2.0 s, m0=1000, mf=900, ratio = 1.111, ln(1.111) = 0.1054
    //   dv_ideal = 300 * 9.80675 * 0.1054 = 310.0 m/s
    //   gravity loss over burn = 9.80675 * 2 = 19.6 m/s
    //   v_burnout = 310.0 - 19.6 = 290.4 m/s
    //   apogee after burn: t_coast = v_burnout / g = 290.4 / 9.80 = 29.6 s
    //   total time to apogee = 2.0 + 29.6 = 31.6 s
    double t_burn = prop->fmass0 / prop->fuel_flow_rate;
    double dv_ideal = prop->spi * G0
                    * std::log(prop->vmass0 / (prop->vmass0 - prop->fmass0));
    double v_burnout = dv_ideal - G0 * t_burn;
    double t_apogee_predict = t_burn + v_burnout / G0;

    std::printf("  Predicted apogee time: %.3f s\n", t_apogee_predict);
    std::printf("  Sim terminated at t = %.3f s\n", icpt->t_terminate);
    std::printf("  trcond = %d (2 = apogee)\n", icpt->trcond);
    std::printf("  final alt = %.1f m  vz_up = %+.4f m/s (~0 at apogee)\n",
                newt->alt, -newt->VBED.z);

    double vz_at_apogee = -newt->VBED.z;

    bool ok =    icpt->trcond == 2
              && std::fabs(icpt->t_terminate - t_apogee_predict) < 0.5
              && std::fabs(vz_at_apogee) < 0.5;    // very close to zero vertical

    std::printf("  %s\n", ok ? "PASS" : "FAIL");

    delete icpt; delete newt; delete kin; delete eul; delete forc; delete prop; delete env;
    return ok ? 0 : 1;
}

int run_time_max() {
    std::printf("\n=== Test 3: Time max termination ===\n");

    Environment* env   = new Environment();
    Propulsion*  prop  = new Propulsion();
    Forces*      forc  = new Forces();
    Euler*       eul   = new Euler();
    Kinematics*  kin   = new Kinematics();
    Newton*      newt  = new Newton();
    Intercept*   icpt  = new Intercept();

    prop->mprop  = 0;
    prop->vmass0 = 100.0;
    forc->FAPB_ext = osk::Vec(0, 0, 0);
    forc->FMB_ext  = osk::Vec(0, 0, 0);
    eul->ppx0 = eul->qqx0 = eul->rrx0 = 0.0;

    // Just sit on the ground for a while
    newt->lonx0 = 0.0; newt->latx0 = 0.0; newt->alt0 = 100000.0; // high alt
    newt->dvbe0 = 0.0;
    kin->psibdx0 = 0.0; kin->thtbdx0 = 0.0; kin->phibdx0 = 0.0;

    // Only time-max check, set t_max = 5 sec
    icpt->check_ground_impact = false;
    icpt->check_apogee        = false;
    icpt->check_time_max      = true;
    icpt->t_max               = 5.0;

    forc->getsFrom(prop);
    newt->getsFrom(env, kin, forc, prop);
    kin->getsFrom(env, newt, eul);
    eul->getsFrom(forc, prop, kin);
    env->getsFrom(newt);
    icpt->getsFrom(env, newt);

    std::vector<osk::Block*> stage0 = { kin, env, prop, forc, newt, eul, icpt };
    std::vector<std::vector<osk::Block*>> stages = { stage0 };

    double dts[] = { 0.01 };
    osk::Sim sim(dts, 60.0, stages);    // Sim's own max is 60; Intercept's is 5
    sim.run();

    std::printf("  t_max set to 5.0 s, Sim's own tmax = 60.0 s\n");
    std::printf("  Sim terminated at t = %.3f s\n", icpt->t_terminate);
    std::printf("  trcond = %d (5 = time max)\n", icpt->trcond);

    bool ok =    icpt->trcond == 5
              && icpt->t_terminate >= 5.0
              && icpt->t_terminate <  5.02;       // within one step of t_max

    std::printf("  %s\n", ok ? "PASS" : "FAIL");

    delete icpt; delete newt; delete kin; delete eul; delete forc; delete prop; delete env;
    return ok ? 0 : 1;
}

} // anon

int main() {
    int fails = 0;
    fails += run_ground_impact();
    fails += run_apogee();
    fails += run_time_max();
    std::printf("\n=== Total failures: %d ===\n", fails);
    return fails == 0 ? 0 : 1;
}
