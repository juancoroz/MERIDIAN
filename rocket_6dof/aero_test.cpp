//  aero_test.cpp  --  Verify Aerodynamics against analytical predictions
//
//  Three tests:
//
//    (1) Static coefficient lookup.  Set up env (rho, Mach, q), set
//        kinematics outputs (alpha=0, phipx=0), confirm Aero produces
//        cx = -ca0(M), cy=cz=0, cll=clm=cln=0.  Then run Forces and
//        confirm FAPB.x = -pdynmc * refa * ca0.
//
//    (2) Terminal velocity.  Drop a 100 kg vehicle from 5 km altitude
//        with refa=1.0 m^2, no thrust, vertical orientation.  At
//        terminal velocity, drag balances weight:
//            v_term = sqrt(2*m*g / (rho*S*Ca))
//        For Ca~0.25 (subsonic), v_term should be ~80 m/s.  Run 30 sec
//        and verify the body has approached terminal velocity.
//
//    (3) Pitching moment sign check.  Set alpha = 4 deg at fixed Mach
//        with a configured vehicle, confirm clm is negative (stable;
//        moment opposes alpha) and matches the table-lookup value
//        within numerical precision.

#include "../osk/osk.h"
#include "environment.h"
#include "newton.h"
#include "kinematics.h"
#include "euler.h"
#include "forces.h"
#include "propulsion.h"
#include "aerodynamics.h"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace rocket6dof;

namespace {

constexpr double G0 = 9.80675445;

// ---- Test 1: static coefficient lookup + Forces aggregation ----
int run_static_lookup() {
    std::printf("\n=== Test 1: Static aero coefficient lookup ===\n");

    Aerodynamics* aero = new Aerodynamics();
    Environment*  env  = new Environment();
    Kinematics*   kin  = new Kinematics();
    Forces*       forc = new Forces();
    Propulsion*   prop = new Propulsion();   // off

    // Configure aero
    aero->maero    = 13;       // 3-stage launch config
    aero->refa     = 1.0;
    aero->refd     = 1.0;
    aero->xcg_ref  = 0.0;
    aero->aero_file = "aero.txt";
    aero->tag_ca0   = "ca0slv3_vs_mach";
    aero->tag_caa   = "caaslv3_vs_mach";
    aero->tag_ca0b  = "ca0bslv3_vs_mach";
    aero->tag_cn0   = "cn0slv3_vs_mach_alpha";
    aero->tag_clm0  = "clm0slv3_vs_mach_alpha";
    aero->tag_clmq  = "clmqslv3_vs_mach";

    aero->getsFrom(env, kin, prop);

    // Force the load via init()
    aero->initCount = 0;
    aero->init();

    // Manually set Env outputs (bypassing the chain since we want a fixed state)
    env->rho    = 1.225;
    env->pdynmc = 0.5 * 1.225 * 100.0 * 100.0;   // q at 100 m/s
    env->vmach  = 0.05;                           // very low Mach
    env->dvba   = 100.0;
    // Kinematics aero angles (zero alpha)
    kin->alppx = 0.0;
    kin->phipx = 0.0;
    kin->alphax = 0.0;
    kin->betax  = 0.0;
    // Propulsion off
    prop->mprop  = 0;
    prop->thrust = 0.0;
    prop->vmass  = 100.0;
    prop->xcg    = 0.0;

    // Tick Aero and Forces once
    aero->update();
    forc->getsFrom(prop, aero, env);
    forc->update();

    // Expected at M=0.05, alpha=0:
    //   ca = ca0(M=0.05) + caa(0.05)*0 + 0 (thrust off, no ca0b)
    //      = 0.2295
    //   cn = cn0(M=0.05, alpha=0) = 0
    //   cx = -ca = -0.2295
    //   cy = -cn*sin(phip) = 0
    //   cz = -cn*cos(phip) = 0
    //   cll = clm = cln = 0
    //
    //   FAPB.x = pdynmc * refa * cx = -0.2295 * (0.5*1.225*100^2) * 1.0
    //          = -0.2295 * 6125.0 = -1405.7 N

    double ca_expected = 0.2295;
    double cx_expected = -ca_expected;
    double FAPBx_expected = env->pdynmc * aero->refa * cx_expected;

    std::printf("  M=%.3f, alpha=%.2f deg\n", env->vmach, kin->alppx);
    std::printf("    aero->ca0  = %.4f  (expect %.4f)\n", aero->ca0, ca_expected);
    std::printf("    aero->ca   = %.4f  (expect %.4f)\n", aero->ca,  ca_expected);
    std::printf("    aero->cna  = %.4f  (expect 0)\n", aero->cna);
    std::printf("    aero->cx   = %.4f  (expect %.4f)\n", aero->cx,  cx_expected);
    std::printf("    aero->cy   = %.4f  (expect 0)\n", aero->cy);
    std::printf("    aero->cz   = %.4f  (expect 0)\n", aero->cz);
    std::printf("    aero->clm  = %.4f  (expect 0)\n", aero->clm);
    std::printf("    FAPB.x     = %.2f N (expect %.2f)\n",
                forc->FAPB.x, FAPBx_expected);

    bool ok =    std::fabs(aero->ca0 - ca_expected)    < 1e-6
              && std::fabs(aero->ca  - ca_expected)    < 1e-6
              && std::fabs(aero->cna)                  < 1e-9
              && std::fabs(aero->cx  - cx_expected)    < 1e-6
              && std::fabs(aero->cy)                   < 1e-9
              && std::fabs(aero->cz)                   < 1e-9
              && std::fabs(forc->FAPB.x - FAPBx_expected) < 1e-3;

    std::printf("  %s\n", ok ? "PASS" : "FAIL");

    delete prop; delete forc; delete kin; delete env; delete aero;
    return ok ? 0 : 1;
}

// ---- Test 2: aero drag deceleration of a coasting body ----
// A simpler test of aero drag than terminal velocity (which puts us
// at gimbal lock for a vertical fall).  Launch a body horizontally at
// 200 m/s at low altitude, no thrust, body x-axis pointing along the
// velocity (alpha = 0).  Aero drag decelerates it.  Compare the
// integrated deceleration over 10 sec against the analytic prediction
// from F = -q*S*Ca with the realized Mach-dependent Ca.
int run_aero_drag() {
    std::printf("\n=== Test 2: Aero drag deceleration ===\n");

    Environment*  env  = new Environment();
    Aerodynamics* aero = new Aerodynamics();
    Propulsion*   prop = new Propulsion();
    Forces*       forc = new Forces();
    Euler*        eul  = new Euler();
    Kinematics*   kin  = new Kinematics();
    Newton*       newt = new Newton();

    aero->maero    = 13;
    aero->refa     = 1.0;
    aero->refd     = 1.0;
    aero->xcg_ref  = 0.0;
    aero->aero_file = "aero.txt";
    aero->tag_ca0   = "ca0slv3_vs_mach";
    aero->tag_caa   = "caaslv3_vs_mach";
    aero->tag_cn0   = "cn0slv3_vs_mach_alpha";
    aero->tag_clm0  = "clm0slv3_vs_mach_alpha";
    aero->tag_clmq  = "clmqslv3_vs_mach";

    prop->mprop  = 0;
    prop->vmass0 = 1000.0;
    prop->fmass0 = 0.0;
    prop->spi    = 0.0;
    prop->fuel_flow_rate = 0.0;
    prop->moi_roll_0  = 100.0;  prop->moi_roll_1  = 100.0;
    prop->moi_trans_0 = 1000.0; prop->moi_trans_1 = 1000.0;

    // Launch horizontally (theta=0) with initial velocity dvbe0 = 200 m/s
    // heading due North (psivdx0 = 0).  Body is aligned with velocity:
    // psi=0, theta=0, phi=0 -> body-x points North in NED, aligned with
    // VBED.  Aero acts as pure drag in the -body-x direction.  Gravity
    // pulls the body down but the body's body-x velocity is preserved
    // (just decays from drag).
    // Use low altitude so air is dense (rho=1.225 at sea level).
    newt->lonx0 = 0.0; newt->latx0 = 0.0; newt->alt0 = 100.0;
    newt->dvbe0 = 200.0; newt->psivdx0 = 0.0; newt->thtvdx0 = 0.0;
    kin->psibdx0 = 0.0; kin->thtbdx0 = 0.0; kin->phibdx0 = 0.0;

    forc->FAPB_ext = osk::Vec(0, 0, 0);
    forc->FMB_ext  = osk::Vec(0, 0, 0);
    eul->ppx0 = eul->qqx0 = eul->rrx0 = 0.0;

    aero->getsFrom(env, kin, prop);
    forc->getsFrom(prop, aero, env);
    newt->getsFrom(env, kin, forc, prop);
    kin->getsFrom(env, newt, eul);
    eul->getsFrom(forc, prop, kin);
    env->getsFrom(newt);

    std::vector<osk::Block*> stage0 = { kin, env, prop, aero, forc, newt, eul };
    std::vector<std::vector<osk::Block*>> stages = { stage0 };

    double t_end = 10.0;
    double dts[] = { 0.01 };
    osk::Sim sim(dts, t_end, stages);
    sim.run();

    // Compute predicted final velocity from numerical integration of
    // dv/dt = -q*S*Ca/m = -(0.5*rho*v^2)*S*Ca/m
    // with Ca varying with Mach.  We use a simple Euler integration here
    // for the analytic prediction (small enough dt that error is < 1%):
    double m = prop->vmass;
    double v_pred = 200.0;
    double dt_pred = 0.001;
    int n_pred = static_cast<int>(t_end / dt_pred);
    for (int i = 0; i < n_pred; ++i) {
        double rho_pred = 1.225;     // approximately constant at low altitudes
        double a_sound  = std::sqrt(1.4 * 287.053 * 288.15);
        double mach_pred = v_pred / a_sound;
        // Use Ca0 lookup approximated by linear interpolation in Mach
        // (matches our 12-point table).  At v=200 m/s, Mach~0.59 -> Ca~0.248.
        // Quick analytic: Ca = Ca0(Mach) (no alpha effect, alpha=0)
        double Ca_pred;
        if      (mach_pred < 0.05) Ca_pred = 0.2295;
        else if (mach_pred < 0.60) Ca_pred = 0.2295 + (0.2475 - 0.2295) * (mach_pred - 0.05) / (0.60 - 0.05);
        else if (mach_pred < 0.95) Ca_pred = 0.2475 + (0.7780 - 0.2475) * (mach_pred - 0.60) / (0.95 - 0.60);
        else                       Ca_pred = 0.7780;
        double q = 0.5 * rho_pred * v_pred * v_pred;
        double drag = q * aero->refa * Ca_pred;
        double accel = -drag / m;
        v_pred += accel * dt_pred;
    }

    std::printf("  Launched at 200 m/s horizontal, drag-only deceleration\n");
    std::printf("  After %g s:\n", t_end);
    std::printf("    total speed dvbe  = %.3f m/s\n", newt->dvbe);
    std::printf("    VBED N component  = %.3f m/s (horizontal, drag-decayed)\n",
                newt->VBED.x);
    std::printf("    VBED D component  = %.3f m/s (vertical, gravity-driven)\n",
                newt->VBED.z);
    std::printf("    Mach              = %.4f\n", env->vmach);
    std::printf("    aero ca           = %.4f\n", aero->ca);
    std::printf("    altitude          = %.1f m  (started 100, gravity pulled down)\n",
                newt->alt);
    std::printf("    body alpha        = %.3f deg (small -> body tracks wind)\n",
                kin->alppx);
    std::printf("  Analytic prediction (pure-drag, no gravity): v_horizontal = %.3f m/s\n",
                v_pred);

    // Compare to horizontal component.  The vehicle's body is statically
    // stable, so the body x-axis follows the velocity vector, keeping
    // alpha small.  The drag deceleration along the (now-tilted) body-x
    // direction is captured by the changing velocity magnitude AND
    // direction.  Decomposing: VBED.x is the surviving northward velocity
    // and should track the pure-drag prediction within RK4 truncation
    // (and the small modulation from non-zero alpha as the body rotates).
    double v_horiz = newt->VBED.x;
    double abs_err = v_horiz - v_pred;
    double rel_err = 100.0 * abs_err / v_pred;
    std::printf("  horizontal error: %+.3f m/s (%+.3f%% rel)\n", abs_err, rel_err);

    // Tolerance: 3% allows for (a) coarse forward-Euler in the analytic
    // baseline, (b) the body's pitch rotation introducing a small
    // dependence on alpha that the analytic ignores.
    bool ok = std::fabs(rel_err) < 3.0;

    std::printf("  %s\n", ok ? "PASS" : "FAIL");

    delete newt; delete kin; delete eul; delete forc; delete prop; delete aero; delete env;
    return ok ? 0 : 1;
}

// ---- Test 3: pitching moment sign and magnitude at alpha = 4 deg ----
int run_pitching_moment() {
    std::printf("\n=== Test 3: Pitching moment at alpha = 4 deg ===\n");

    Aerodynamics* aero = new Aerodynamics();
    Environment*  env  = new Environment();
    Kinematics*   kin  = new Kinematics();
    Propulsion*   prop = new Propulsion();
    Forces*       forc = new Forces();

    aero->maero    = 13;
    aero->refa     = 1.0;
    aero->refd     = 5.0;
    aero->xcg_ref  = 2.5;
    aero->aero_file = "aero.txt";
    aero->tag_ca0   = "ca0slv3_vs_mach";
    aero->tag_caa   = "caaslv3_vs_mach";
    aero->tag_cn0   = "cn0slv3_vs_mach_alpha";
    aero->tag_clm0  = "clm0slv3_vs_mach_alpha";
    aero->tag_clmq  = "clmqslv3_vs_mach";

    aero->getsFrom(env, kin, prop);
    aero->initCount = 0;
    aero->init();

    // Force a state: Mach=0.95 (transonic), alpha=4 deg, no aero roll
    env->rho    = 0.5;
    env->dvba   = 320.0;       // ~ Mach 0.95 in sea-level air
    env->vmach  = 0.95;
    env->pdynmc = 0.5 * 0.5 * 320.0 * 320.0;

    kin->alppx = 4.0;
    kin->phipx = 0.0;     // pitch plane is body-x/body-z

    prop->mprop = 0;
    prop->thrust = 0.0;
    prop->vmass  = 100.0;
    prop->xcg    = 2.5;   // same as xcg_ref so no CG-shift correction

    aero->update();

    // Expected from table (M=0.95, alpha=4 deg from cn0slv3 and clm0slv3):
    // cn0 row for Mach=0.95 (index 2 in our 12-row table, 0-based)
    //   columns: 0.0000 0.0648 0.1429 0.2349 0.3414 0.4626
    // alpha=4 deg is the column index 2 (0,2,4,6,8,10 -> indices 0..5)
    //   cn0(M=0.95, alpha=4) = 0.1429
    // clm0(M=0.95, alpha=4) = -0.79
    double cn0_exp = 0.1429;
    double clm0_exp = -0.79;

    // With xcg = xcg_ref, no CG correction:
    //   cna = cn0 = 0.1429
    //   clma = clm0 = -0.79
    // Aero roll phip=0 -> sphip=0, cphip=1:
    //   cx = -ca   (depends on ca0(0.95) + caa(0.95)*4)
    //   cy = -cn*0 = 0
    //   cz = -cn*1 = -0.1429
    //   clm = clma*1 = -0.79
    //   cln = -clma*0 = 0

    // Also predict ca from table values:
    //   ca0(0.95) = 0.7780, caa(0.95) = -0.0053, alpha=4
    //   ca = 0.7780 + (-0.0053)*4 = 0.7780 - 0.0212 = 0.7568
    double ca_exp = 0.7780 + (-0.0053) * 4.0;

    std::printf("  M=%.2f, alpha=%.1f deg\n", env->vmach, kin->alppx);
    std::printf("    ca0  = %.4f  (expect %.4f)\n", aero->ca0, 0.7780);
    std::printf("    ca   = %.4f  (expect %.4f)\n", aero->ca,  ca_exp);
    std::printf("    cn0  = %.4f  (expect %.4f)\n", aero->cn0, cn0_exp);
    std::printf("    clm0 = %.4f  (expect %.4f)\n", aero->clm0_v, clm0_exp);
    std::printf("    clm  = %.4f  (expect %.4f, statically stable < 0)\n",
                aero->clm, clm0_exp);

    bool ok =    std::fabs(aero->ca0 - 0.7780)    < 1e-4
              && std::fabs(aero->ca  - ca_exp)    < 1e-4
              && std::fabs(aero->cn0 - cn0_exp)   < 1e-4
              && std::fabs(aero->clm0_v - clm0_exp) < 1e-4
              && std::fabs(aero->clm  - clm0_exp) < 1e-4   // since phip=0
              && aero->clm < 0;                              // sign: stable

    std::printf("  %s\n", ok ? "PASS" : "FAIL");

    delete forc; delete prop; delete kin; delete env; delete aero;
    return ok ? 0 : 1;
}

} // anon

int main() {
    int fails = 0;
    fails += run_static_lookup();
    fails += run_aero_drag();
    fails += run_pitching_moment();
    std::printf("\n=== Total failures: %d ===\n", fails);
    return fails == 0 ? 0 : 1;
}
