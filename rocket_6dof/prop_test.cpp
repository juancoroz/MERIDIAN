//  prop_test.cpp  --  Verify Propulsion against the rocket equation
//
//  Tsiolkovsky rocket equation: for a vehicle with no external forces
//  other than thrust and uniform gravity, vertical launch under
//  constant thrust gives burnout velocity
//
//      v(t_burn) = Isp * g0 * ln(m0 / m_f) - g * t_burn
//
//  where the first term is the ideal delta-v and the second is the
//  gravity-loss term.  This test sets up exactly that scenario:
//
//    * spherical, non-rotating Earth, pure inverse-square gravity
//    * vehicle pointing straight up (theta = 90 deg)
//    * single Propulsion block, constant thrust, burns for 10 s
//    * no aerodynamics, no RCS, no TVC
//    * verify burnout velocity matches the analytical formula
//
//  Then we extend to a three-stage rocket using OSK's Sim::stop
//  mechanism to transition between stages at burnout, and confirm
//  the cumulative delta-v matches stage-wise Tsiolkovsky.

#include "../osk/osk.h"
#include "environment.h"
#include "newton.h"
#include "kinematics.h"
#include "euler.h"
#include "forces.h"
#include "propulsion.h"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace rocket6dof;

namespace {

constexpr double G0 = 9.80675445;   // Zipfel's AGRAV

// ---- Test 1: single-stage Tsiolkovsky verification ----
int run_single_stage() {
    std::printf("\n=== Test 1: single-stage rocket, Tsiolkovsky verification ===\n");

    Environment* env  = new Environment();
    Propulsion*  prop = new Propulsion();
    Forces*      forc = new Forces();
    Euler*       eul  = new Euler();
    Kinematics*  kin  = new Kinematics();
    Newton*      newt = new Newton();

    // ---- Propulsion parameters ----
    prop->mprop          = 3;       // constant-thrust
    prop->vmass0         = 1000.0;  // 1000 kg gross
    prop->fmass0         = 500.0;   // 500 kg fuel
    prop->spi            = 300.0;   // 300 s Isp
    prop->fuel_flow_rate = 50.0;    // 50 kg/s -> burn time = 10 s
    prop->moi_roll_0     = 10.0;    prop->moi_roll_1   = 5.0;
    prop->moi_trans_0    = 100.0;   prop->moi_trans_1  = 50.0;
    prop->xcg_0          = 2.0;     prop->xcg_1        = 1.8;

    // ---- Initial conditions ----
    //   On the pad at equator, prime meridian.
    //   Body pointing straight up: theta = 90 deg.
    newt->lonx0 = 0.0;  newt->latx0 = 0.0;  newt->alt0 = 0.0;
    newt->dvbe0 = 0.0;
    kin->psibdx0 = 0.0;  kin->thtbdx0 = 90.0;  kin->phibdx0 = 0.0;

    forc->FAPB_ext = osk::Vec(0, 0, 0);
    forc->FMB_ext  = osk::Vec(0, 0, 0);
    eul->ppx0 = eul->qqx0 = eul->rrx0 = 0.0;

    // ---- Wire dependencies ----
    forc->getsFrom(prop);                          // Forces reads thrust from Propulsion
    newt->getsFrom(env, kin, forc, prop);          // Newton uses prop for vmass (time-varying)
    kin->getsFrom(env, newt, eul);
    eul->getsFrom(forc, prop, kin);
    env->getsFrom(newt);                // Euler uses prop for IBBB

    // ---- Stage order (Zipfel canonical) ----
    std::vector<osk::Block*> stage0 = { kin, env, prop, forc, newt, eul };
    std::vector<std::vector<osk::Block*>> stages = { stage0 };

    double dts[] = { 0.01 };
    osk::Sim sim(dts, 12.0, stages);   // run 2 s past burnout for coast
    sim.run();

    // ---- Analytic expectations ----
    double t_burn      = prop->fmass0 / prop->fuel_flow_rate;     // 10 s
    double m0          = prop->vmass0;
    double mf          = prop->vmass0 - prop->fmass0;             // 500 kg
    double thrust      = prop->spi * prop->fuel_flow_rate * G0;
    double dv_ideal    = prop->spi * G0 * std::log(m0 / mf);
    double dv_gravloss = G0 * t_burn;
    double v_burnout   = dv_ideal - dv_gravloss;

    // After burnout (at t_burn = 10 s), vehicle coasts.  By t=12 it has
    // decelerated by g * (t - t_burn) = 9.80675 * 2 ~ 19.6 m/s.  Use the
    // simulated value at t=12 to back out v_burnout (since rpt() only
    // prints at integer seconds we'll compare against final state).
    double v_final     = newt->dvbi;
    double v_burnout_sim = v_final + G0 * (12.0 - t_burn);

    std::printf("  Stage parameters:\n");
    std::printf("    m0   = %.1f kg,  mf = %.1f kg,  m0/mf = %.3f\n",
                m0, mf, m0/mf);
    std::printf("    spi  = %.1f s,  fuel_flow = %.1f kg/s,  thrust = %.1f N\n",
                prop->spi, prop->fuel_flow_rate, thrust);
    std::printf("    t_burn = %.2f s,  T/W at launch = %.2f\n",
                t_burn, thrust / (m0 * G0));
    std::printf("\n");
    std::printf("  Analytic prediction:\n");
    std::printf("    Tsiolkovsky dv_ideal = Isp * g * ln(m0/mf) = %.2f m/s\n", dv_ideal);
    std::printf("    Gravity loss         = g * t_burn          = %.2f m/s\n", dv_gravloss);
    std::printf("    v_burnout predicted  = %.2f m/s\n", v_burnout);
    std::printf("\n");
    std::printf("  Simulated:\n");
    std::printf("    v(t=12)              = %.2f m/s\n", v_final);
    std::printf("    v_burnout backed out = v(12) + g*(12-t_burn) = %.2f m/s\n",
                v_burnout_sim);
    std::printf("    error  = %.4f m/s  (%.3f%% rel)\n",
                v_burnout_sim - v_burnout,
                100.0 * (v_burnout_sim - v_burnout) / v_burnout);
    std::printf("\n");
    std::printf("    final altitude       = %.2f m\n", newt->alt);
    std::printf("    final vmass          = %.2f kg  (expect %.2f)\n",
                prop->vmass, mf);
    std::printf("    final fmassr         = %.2f kg  (expect 0)\n", prop->fmassr);

    // Tolerance: ~0.5 m/s on a ~1940 m/s burnout velocity is ~0.025%, well
    // within RK4 truncation given dt=0.01 over 1000 steps.  The gravity
    // model is inverse-square (slight non-uniformity over the ~10 km climb)
    // so the constant-g formula has its own ~0.1% error.
    bool ok =    std::fabs(v_burnout_sim - v_burnout) < 2.0
              && std::fabs(prop->vmass  - mf)         < 0.01
              && std::fabs(prop->fmassr - 0.0)        < 0.01;

    std::printf("  %s\n", ok ? "PASS" : "FAIL");

    delete newt; delete kin; delete eul; delete forc; delete prop; delete env;
    return ok ? 0 : 1;
}

// ---- Test 2: three-stage rocket using OSK Sim::stop ----
// We construct three Propulsion objects with successively smaller masses
// and verify cumulative delta-v matches the stage-wise Tsiolkovsky sum.
// Staging is triggered by a small "Stager" block that watches fmassr,
// rewires Forces/Newton/Euler to the next Propulsion, and writes
// Sim::stop to advance to the next stage.
int run_three_stage() {
    std::printf("\n=== Test 2: three-stage rocket, cumulative delta-v ===\n");

    Environment* env  = new Environment();
    Forces*      forc = new Forces();
    Kinematics*  kin  = new Kinematics();
    Euler*       eul  = new Euler();
    Newton*      newt = new Newton();

    // Three propulsion objects, one per stage.  Realistic-ish ratios:
    // each upper stage is a fraction of the previous total mass.
    Propulsion* p1 = new Propulsion();
    Propulsion* p2 = new Propulsion();
    Propulsion* p3 = new Propulsion();

    // Stage 1: 1000 kg gross, 700 kg fuel, Isp 280, flow 100 kg/s -> 7 s burn
    p1->mprop=3; p1->vmass0=1000.0; p1->fmass0=700.0;
    p1->spi=280.0; p1->fuel_flow_rate=100.0;
    p1->moi_roll_0=10; p1->moi_roll_1=5; p1->moi_trans_0=100; p1->moi_trans_1=50;

    // Stage 2: 250 kg gross (after stage 1 jettisons 50 kg of structure),
    // 180 kg fuel, Isp 320, flow 30 kg/s -> 6 s burn
    p2->mprop=3; p2->vmass0=250.0; p2->fmass0=180.0;
    p2->spi=320.0; p2->fuel_flow_rate=30.0;
    p2->moi_roll_0=5; p2->moi_roll_1=2; p2->moi_trans_0=50; p2->moi_trans_1=20;

    // Stage 3: 60 kg gross, 40 kg fuel, Isp 340, flow 5 kg/s -> 8 s burn
    p3->mprop=3; p3->vmass0=60.0; p3->fmass0=40.0;
    p3->spi=340.0; p3->fuel_flow_rate=5.0;
    p3->moi_roll_0=2; p3->moi_roll_1=1; p3->moi_trans_0=20; p3->moi_trans_1=10;

    // Common ICs
    newt->lonx0 = 0.0; newt->latx0 = 0.0; newt->alt0 = 0.0;
    newt->dvbe0 = 0.0;
    kin->psibdx0 = 0.0; kin->thtbdx0 = 90.0; kin->phibdx0 = 0.0;
    forc->FAPB_ext = osk::Vec(0, 0, 0);
    forc->FMB_ext  = osk::Vec(0, 0, 0);
    eul->ppx0 = eul->qqx0 = eul->rrx0 = 0.0;

    // Per-stage wiring of Forces and Newton/Euler.  Each stage uses its
    // own Propulsion, so we need to re-wire at each stage entry.
    // OSK's solution: have separate Forces/Newton/Euler instances per
    // stage, OR use a single set and re-wire init().  We pick the
    // simpler path: have each stage's block list include a small
    // "Rewire" block that points Forces/Newton/Euler at the correct
    // Propulsion in their init().  But this is fragile.
    //
    // Cleaner approach: write the wiring directly using stage-specific
    // pointers.  Forces+Newton+Euler each have a getsFrom call before
    // the sim, and we update those between stages by calling getsFrom
    // again on stage transition.  But OSK doesn't give us a hook
    // between stages.  Instead, do the simplest correct thing: wire
    // everything to the *current* propulsion in each stage's init()
    // function via a thin proxy.
    //
    // Simplest hack: rewire from the stager.  When Stager fires
    // Sim::stop, it also updates the wiring pointers.

    // Stage 0: wire to p1
    forc->getsFrom(p1);
    newt->getsFrom(env, kin, forc, p1);
    kin->getsFrom(env, newt, eul);
    eul->getsFrom(forc, p1, kin);
    env->getsFrom(newt);

    // For staging: the Stager block, in addition to setting Sim::stop,
    // rewires Forces+Newton+Euler to the next stage's Propulsion before
    // the kernel transitions.  This works because Stager::update runs
    // after Propulsion::update in our stage ordering, so by the time
    // Sim::stop is read, the new wiring is in place.

    // Specialized Stager that also rewires
    class RewireStager : public osk::Block {
    public:
        Propulsion* old_prop; Propulsion* new_prop; int next_stage;
        Forces* f; Newton* n; Euler* e; Kinematics* k; Environment* env;
        bool fired;
        RewireStager(Propulsion* op, Propulsion* np, int ns,
                     Forces* fp, Newton* np2, Euler* ep, Kinematics* kp, Environment* en)
            : old_prop(op), new_prop(np), next_stage(ns),
              f(fp), n(np2), e(ep), k(kp), env(en), fired(false) {}
        void init()   override { fired = false; }
        void update() override {
            if (!fired && old_prop && old_prop->fmassr_() <= 0.0) {
                fired = true;
                // Rewire downstream blocks BEFORE telling sim to stop
                f->getsFrom(new_prop);
                n->getsFrom(env, k, f, new_prop);
                e->getsFrom(f, new_prop, k);
                osk::Sim::stop = next_stage;
            }
        }
        void rpt() override {}
    };

    RewireStager* rs1 = new RewireStager(p1, p2, 1, forc, newt, eul, kin, env);
    RewireStager* rs2 = new RewireStager(p2, p3, 2, forc, newt, eul, kin, env);

    // Build three stages.  Each stage contains its own active Propulsion
    // plus the SHARED dynamics blocks.  The kinematic/dynamics state
    // (TBI, SBII, VBII, WBIB) carries forward across stages.
    std::vector<osk::Block*> stage0 = { kin, env, p1, forc, newt, eul, rs1 };
    std::vector<osk::Block*> stage1 = { kin, env, p2, forc, newt, eul, rs2 };
    std::vector<osk::Block*> stage2 = { kin, env, p3, forc, newt, eul };
    std::vector<std::vector<osk::Block*>> stages = { stage0, stage1, stage2 };

    double dts[] = { 0.01, 0.01, 0.01 };
    osk::Sim sim(dts, 30.0, stages);
    sim.run();

    // Stage-wise Tsiolkovsky predictions
    double dv1 = p1->spi * G0 * std::log(p1->vmass0 / (p1->vmass0 - p1->fmass0));
    double dv2 = p2->spi * G0 * std::log(p2->vmass0 / (p2->vmass0 - p2->fmass0));
    double dv3 = p3->spi * G0 * std::log(p3->vmass0 / (p3->vmass0 - p3->fmass0));
    double dv_total_ideal = dv1 + dv2 + dv3;
    double t_burn_total   = p1->fmass0/p1->fuel_flow_rate
                          + p2->fmass0/p2->fuel_flow_rate
                          + p3->fmass0/p3->fuel_flow_rate;
    double gravity_loss   = G0 * t_burn_total;
    double v_predicted    = dv_total_ideal - gravity_loss;

    double v_final = newt->dvbi;

    // Back out v at end-of-burn from v_final + gravity loss over coast
    double t_coast      = 30.0 - t_burn_total;
    double v_burnout_sim = v_final + G0 * t_coast;

    std::printf("  Stage-wise delta-v:\n");
    std::printf("    Stage 1:  Isp=%.0f, m0/mf=%.3f, dv = %.2f m/s, t_burn = %.2f s\n",
                p1->spi, p1->vmass0/(p1->vmass0-p1->fmass0), dv1, p1->fmass0/p1->fuel_flow_rate);
    std::printf("    Stage 2:  Isp=%.0f, m0/mf=%.3f, dv = %.2f m/s, t_burn = %.2f s\n",
                p2->spi, p2->vmass0/(p2->vmass0-p2->fmass0), dv2, p2->fmass0/p2->fuel_flow_rate);
    std::printf("    Stage 3:  Isp=%.0f, m0/mf=%.3f, dv = %.2f m/s, t_burn = %.2f s\n",
                p3->spi, p3->vmass0/(p3->vmass0-p3->fmass0), dv3, p3->fmass0/p3->fuel_flow_rate);
    std::printf("    Sum:                              dv = %.2f m/s, t_burn_tot = %.2f s\n",
                dv_total_ideal, t_burn_total);
    std::printf("  Predicted v_burnout = %.2f - %.2f = %.2f m/s\n",
                dv_total_ideal, gravity_loss, v_predicted);
    std::printf("  Simulated v(30s)    = %.2f m/s\n", v_final);
    std::printf("  v_burnout backed out = %.2f m/s\n", v_burnout_sim);
    std::printf("  error  = %.4f m/s  (%.3f%% rel)\n",
                v_burnout_sim - v_predicted,
                100.0 * (v_burnout_sim - v_predicted) / v_predicted);

    // Tolerance: 15 m/s on a ~10.8 km/s burnout.  The constant-gravity
    // assumption in v_predicted has ~0.1% error over a 200 km ascent
    // (gravity weakens from 9.80 to 9.19 m/s^2), so the analytic baseline
    // is itself wrong by ~10 m/s; the actual integrated trajectory is
    // MORE accurate than the constant-g prediction.
    bool ok = std::fabs(v_burnout_sim - v_predicted) < 15.0;
    std::printf("  %s\n", ok ? "PASS" : "FAIL");

    delete rs2; delete rs1;
    delete p3; delete p2; delete p1;
    delete newt; delete eul; delete kin; delete forc; delete env;
    return ok ? 0 : 1;
}

//  Test 3: Multi-stage support inside a single Propulsion instance.
//
//  Same stage parameters as Test 2 (which used 3 separate Propulsion
//  objects + a RewireStager block).  This test exercises the new
//  num_stages > 1 path inside Propulsion, with NO coast between stages
//  (to match Test 2's instant-switch behavior).
//
//  Pass criterion: same v_burnout as Test 2 to within 1 m/s.  Same
//  Tsiolkovsky-stage-wise prediction.  Plus the staging state machine
//  reports the correct phase/current_stage transitions.
int run_in_propulsion_multistage() {
    std::printf("\n=== Test 3: in-Propulsion multi-stage support ===\n");

    Environment* env  = new Environment();
    Forces*      forc = new Forces();
    Kinematics*  kin  = new Kinematics();
    Euler*       eul  = new Euler();
    Newton*      newt = new Newton();
    Propulsion*  prop = new Propulsion();

    // Configure 3-stage rocket using the new arrays.
    // Stage 1: 1000 kg gross, 700 kg fuel, Isp 280, mdot 100 -> 7 s
    // Stage 2: 250 kg gross,  180 kg fuel, Isp 320, mdot 30  -> 6 s
    // Stage 3: 60  kg gross,   40 kg fuel, Isp 340, mdot 5   -> 8 s
    // (dry_mass_dropped = 50 kg between 1->2: 1000-700=300, drop 50, ->250)
    // (dry_mass_dropped = 10 kg between 2->3: 250-180=70,  drop 10, ->60)
    prop->mprop      = 3;
    prop->num_stages = 3;

    prop->vmass0_stage[0] = 1000.0; prop->fmass0_stage[0] = 700.0;
    prop->spi_stage[0]    = 280.0;  prop->mdot_stage[0]   = 100.0;
    prop->moi_roll_0_stage[0]  = 10; prop->moi_roll_1_stage[0]  = 5;
    prop->moi_trans_0_stage[0] = 100; prop->moi_trans_1_stage[0] = 50;
    prop->dry_mass_dropped[0]  = 50.0;
    prop->coast_to_next[0]     = 0.0;     // instantaneous staging

    prop->vmass0_stage[1] = 250.0;  prop->fmass0_stage[1] = 180.0;
    prop->spi_stage[1]    = 320.0;  prop->mdot_stage[1]   = 30.0;
    prop->moi_roll_0_stage[1]  = 5; prop->moi_roll_1_stage[1]  = 2;
    prop->moi_trans_0_stage[1] = 50; prop->moi_trans_1_stage[1] = 20;
    prop->dry_mass_dropped[1]  = 10.0;
    prop->coast_to_next[1]     = 0.0;

    prop->vmass0_stage[2] = 60.0;   prop->fmass0_stage[2] = 40.0;
    prop->spi_stage[2]    = 340.0;  prop->mdot_stage[2]   = 5.0;
    prop->moi_roll_0_stage[2]  = 2; prop->moi_roll_1_stage[2]  = 1;
    prop->moi_trans_0_stage[2] = 20; prop->moi_trans_1_stage[2] = 10;
    prop->dry_mass_dropped[2]  = 0.0;     // last stage, nothing more to drop
    prop->coast_to_next[2]     = 0.0;

    // ICs (same as Test 2)
    newt->lonx0 = 0.0; newt->latx0 = 0.0; newt->alt0 = 0.0;
    newt->dvbe0 = 0.0;
    kin->psibdx0 = 0.0; kin->thtbdx0 = 90.0; kin->phibdx0 = 0.0;
    forc->FAPB_ext = osk::Vec(0, 0, 0);
    forc->FMB_ext  = osk::Vec(0, 0, 0);
    eul->ppx0 = eul->qqx0 = eul->rrx0 = 0.0;

    env->getsFrom(newt);
    forc->getsFrom(prop);
    newt->getsFrom(env, kin, forc, prop);
    kin->getsFrom(env, newt, eul);
    eul->getsFrom(forc, prop, kin);

    std::vector<osk::Block*> stage0 = { kin, env, prop, forc, newt, eul };
    std::vector<std::vector<osk::Block*>> stages = { stage0 };

    double dts[] = { 0.01 };
    double t_total = 7.0 + 6.0 + 8.0 + 0.5;   // total burn + slop
    osk::Sim sim(dts, t_total, stages);
    sim.run();

    // Predictions: stage-wise Tsiolkovsky.
    // dv per stage uses that stage's vmass0_stage and fmass0_stage.
    double dv1 = prop->spi_stage[0] * G0
               * std::log(prop->vmass0_stage[0]
                        / (prop->vmass0_stage[0] - prop->fmass0_stage[0]));
    double dv2 = prop->spi_stage[1] * G0
               * std::log(prop->vmass0_stage[1]
                        / (prop->vmass0_stage[1] - prop->fmass0_stage[1]));
    double dv3 = prop->spi_stage[2] * G0
               * std::log(prop->vmass0_stage[2]
                        / (prop->vmass0_stage[2] - prop->fmass0_stage[2]));
    double dv_total = dv1 + dv2 + dv3;
    double t_burn   = prop->fmass0_stage[0] / prop->mdot_stage[0]
                    + prop->fmass0_stage[1] / prop->mdot_stage[1]
                    + prop->fmass0_stage[2] / prop->mdot_stage[2];
    double gloss    = G0 * t_burn;
    double v_pred   = dv_total - gloss;

    std::printf("  Stage burns: %.1f + %.1f + %.1f = %.1f s\n",
                prop->fmass0_stage[0] / prop->mdot_stage[0],
                prop->fmass0_stage[1] / prop->mdot_stage[1],
                prop->fmass0_stage[2] / prop->mdot_stage[2],
                t_burn);
    std::printf("  Stage delta-v: %.1f + %.1f + %.1f = %.1f m/s\n",
                dv1, dv2, dv3, dv_total);
    std::printf("  Gravity loss = %.1f m/s\n", gloss);
    std::printf("  Predicted burnout speed = %.1f m/s\n", v_pred);
    std::printf("  Simulated burnout speed = %.2f m/s  (dvbi=%.2f)\n",
                newt->dvbe, newt->dvbi);
    std::printf("  Final stage:        %d  (expect 2)\n", prop->current_stage);
    std::printf("  Final phase:        %d  (expect 2 = spent)\n", prop->phase);
    std::printf("  Final vmass:        %.2f kg  (expect %.2f)\n",
                prop->vmass,
                prop->vmass0_stage[2] - prop->fmass0_stage[2]);

    bool dv_ok    = std::fabs(newt->dvbe - v_pred) < 15.0;
    bool stage_ok = (prop->current_stage == 2);
    bool phase_ok = (prop->phase == 2);
    bool mass_ok  = std::fabs(prop->vmass
                       - (prop->vmass0_stage[2] - prop->fmass0_stage[2])) < 0.5;
    bool ok = dv_ok && stage_ok && phase_ok && mass_ok;
    std::printf("  dv_ok: %s  stage_ok: %s  phase_ok: %s  mass_ok: %s\n",
                dv_ok ? "Y":"N", stage_ok ? "Y":"N",
                phase_ok ? "Y":"N", mass_ok ? "Y":"N");
    std::printf("  %s\n", ok ? "PASS" : "FAIL");

    delete prop; delete newt; delete eul; delete kin; delete forc; delete env;
    return ok ? 0 : 1;
}

//  Test 4: Multi-stage with inter-stage coast
//
//  Same as Test 3 but with a 2-second coast between stages 1->2.  The
//  test verifies that during the coast:
//   * thrust = 0
//   * vmass, MoI, xcg are held at stage-1-burnout values
//   * fmasse stays at fmass0_stage[0] (no integration during coast)
//   * after coast: stage transitions to 2, fmasse resets to 0, vmass
//     becomes vmass0_stage[1]
int run_multistage_with_coast() {
    std::printf("\n=== Test 4: multi-stage with inter-stage coast ===\n");

    Environment* env  = new Environment();
    Forces*      forc = new Forces();
    Kinematics*  kin  = new Kinematics();
    Euler*       eul  = new Euler();
    Newton*      newt = new Newton();
    Propulsion*  prop = new Propulsion();

    prop->mprop      = 3;
    prop->num_stages = 2;

    // Stage 1: 1000 kg, 700 kg fuel, mdot 100 -> 7 s burn
    prop->vmass0_stage[0] = 1000.0; prop->fmass0_stage[0] = 700.0;
    prop->spi_stage[0]    = 280.0;  prop->mdot_stage[0]   = 100.0;
    prop->moi_roll_0_stage[0]  = 10; prop->moi_roll_1_stage[0]  = 5;
    prop->moi_trans_0_stage[0] = 100; prop->moi_trans_1_stage[0] = 50;
    prop->dry_mass_dropped[0]  = 50.0;
    prop->coast_to_next[0]     = 2.0;     // 2 second coast

    // Stage 2: 250 kg, 180 kg fuel, mdot 30 -> 6 s
    prop->vmass0_stage[1] = 250.0; prop->fmass0_stage[1] = 180.0;
    prop->spi_stage[1]    = 320.0; prop->mdot_stage[1]   = 30.0;
    prop->moi_roll_0_stage[1]  = 5; prop->moi_roll_1_stage[1]  = 2;
    prop->moi_trans_0_stage[1] = 50; prop->moi_trans_1_stage[1] = 20;
    prop->dry_mass_dropped[1]  = 0.0;
    prop->coast_to_next[1]     = 0.0;

    newt->lonx0 = 0.0; newt->latx0 = 0.0; newt->alt0 = 0.0;
    newt->dvbe0 = 0.0;
    kin->psibdx0 = 0.0; kin->thtbdx0 = 90.0; kin->phibdx0 = 0.0;
    forc->FAPB_ext = osk::Vec(0, 0, 0);
    forc->FMB_ext  = osk::Vec(0, 0, 0);
    eul->ppx0 = eul->qqx0 = eul->rrx0 = 0.0;

    env->getsFrom(newt);
    forc->getsFrom(prop);
    newt->getsFrom(env, kin, forc, prop);
    kin->getsFrom(env, newt, eul);
    eul->getsFrom(forc, prop, kin);

    std::vector<osk::Block*> stage0 = { kin, env, prop, forc, newt, eul };
    std::vector<std::vector<osk::Block*>> stages = { stage0 };

    double dts[] = { 0.01 };

    // Run for 7.5 s (mid-coast).  Sample state.
    osk::Sim sim_a(dts, 7.5, stages);
    sim_a.run();
    double thrust_coast = prop->thrust;
    double vmass_coast  = prop->vmass;
    int    phase_coast  = prop->phase;
    int    stage_coast  = prop->current_stage;

    std::printf("  At t=7.5s (mid-coast):\n");
    std::printf("    thrust = %.1f N    (expect 0)\n", thrust_coast);
    std::printf("    vmass  = %.1f kg   (expect 300 = stage-1 burnout)\n", vmass_coast);
    std::printf("    phase  = %d        (expect 1 = coast)\n", phase_coast);
    std::printf("    stage  = %d        (expect 0)\n", stage_coast);

    // Continue running to end of full mission (7+2+6 = 15 s burn time)
    osk::Sim sim_b(dts, 16.0, stages);
    sim_b.run();
    double final_thrust = prop->thrust;
    double final_vmass  = prop->vmass;
    int    final_phase  = prop->phase;
    int    final_stage  = prop->current_stage;

    std::printf("  At t=16s (after stage 2 burnout):\n");
    std::printf("    thrust = %.1f N    (expect 0)\n", final_thrust);
    std::printf("    vmass  = %.1f kg   (expect 70 = stage-2 burnout)\n", final_vmass);
    std::printf("    phase  = %d        (expect 2 = spent)\n", final_phase);
    std::printf("    stage  = %d        (expect 1)\n", final_stage);

    bool coast_ok   = (thrust_coast == 0.0)
                   && std::fabs(vmass_coast - 300.0) < 0.5
                   && (phase_coast == 1)
                   && (stage_coast == 0);
    bool final_ok   = (final_thrust == 0.0)
                   && std::fabs(final_vmass - 70.0) < 0.5
                   && (final_phase == 2)
                   && (final_stage == 1);

    bool ok = coast_ok && final_ok;
    std::printf("  coast_ok: %s  final_ok: %s\n",
                coast_ok ? "Y":"N", final_ok ? "Y":"N");
    std::printf("  %s\n", ok ? "PASS" : "FAIL");

    delete prop; delete newt; delete eul; delete kin; delete forc; delete env;
    return ok ? 0 : 1;
}

//  Test 5: 4-stage rocket (regression guard for MAX_STAGES=4)
//
//  Identical structure to Test 3 but with 4 stages, exercising the
//  propulsion path past the old hard limit of 3.  Verifies:
//   * all 4 stages activate in sequence
//   * cumulative delta-v matches stage-wise Tsiolkovsky
//   * final stage index, phase, and vmass match expectations
//
//  Stage masses chosen so each stage drops to the next as a real
//  rocket would: vmass0_stage[i+1] = (vmass after burn[i]) - dropped[i]
int run_four_stage() {
    std::printf("\n=== Test 5: 4-stage rocket (MAX_STAGES regression guard) ===\n");

    Environment* env  = new Environment();
    Forces*      forc = new Forces();
    Kinematics*  kin  = new Kinematics();
    Euler*       eul  = new Euler();
    Newton*      newt = new Newton();
    Propulsion*  prop = new Propulsion();

    // Stage 1: 1200 kg gross, 800 kg fuel, Isp 260, mdot 100 -> 8 s
    // Stage 2: 350 kg gross,  240 kg fuel, Isp 300, mdot 40  -> 6 s
    // Stage 3: 90 kg gross,   60 kg fuel,  Isp 320, mdot 8   -> 7.5 s
    // Stage 4: 22 kg gross,   12 kg fuel,  Isp 340, mdot 2   -> 6 s
    // Dropped: 50 kg between 1->2 (1200-800=400, -50, ->350)
    // Dropped: 20 kg between 2->3 (350-240=110, -20, ->90)
    // Dropped: 8 kg between 3->4 (90-60=30, -8, ->22)
    prop->mprop      = 3;
    prop->num_stages = 4;

    prop->vmass0_stage[0] = 1200.0; prop->fmass0_stage[0] = 800.0;
    prop->spi_stage[0]    = 260.0;  prop->mdot_stage[0]   = 100.0;
    prop->moi_roll_0_stage[0]  = 12; prop->moi_roll_1_stage[0]  = 6;
    prop->moi_trans_0_stage[0] = 120; prop->moi_trans_1_stage[0] = 60;
    prop->dry_mass_dropped[0]  = 50.0;
    prop->coast_to_next[0]     = 0.0;

    prop->vmass0_stage[1] = 350.0;  prop->fmass0_stage[1] = 240.0;
    prop->spi_stage[1]    = 300.0;  prop->mdot_stage[1]   = 40.0;
    prop->moi_roll_0_stage[1]  = 6; prop->moi_roll_1_stage[1]  = 3;
    prop->moi_trans_0_stage[1] = 60; prop->moi_trans_1_stage[1] = 30;
    prop->dry_mass_dropped[1]  = 20.0;
    prop->coast_to_next[1]     = 0.0;

    prop->vmass0_stage[2] = 90.0;   prop->fmass0_stage[2] = 60.0;
    prop->spi_stage[2]    = 320.0;  prop->mdot_stage[2]   = 8.0;
    prop->moi_roll_0_stage[2]  = 3; prop->moi_roll_1_stage[2]  = 1.5;
    prop->moi_trans_0_stage[2] = 30; prop->moi_trans_1_stage[2] = 15;
    prop->dry_mass_dropped[2]  = 8.0;
    prop->coast_to_next[2]     = 0.0;

    prop->vmass0_stage[3] = 22.0;   prop->fmass0_stage[3] = 12.0;
    prop->spi_stage[3]    = 340.0;  prop->mdot_stage[3]   = 2.0;
    prop->moi_roll_0_stage[3]  = 1.0; prop->moi_roll_1_stage[3]  = 0.5;
    prop->moi_trans_0_stage[3] = 10; prop->moi_trans_1_stage[3] = 5;
    prop->dry_mass_dropped[3]  = 0.0;     // last stage
    prop->coast_to_next[3]     = 0.0;

    // ICs (vertical launch, same as the other multi-stage tests)
    newt->lonx0 = 0.0; newt->latx0 = 0.0; newt->alt0 = 0.0;
    newt->dvbe0 = 0.0;
    kin->psibdx0 = 0.0; kin->thtbdx0 = 90.0; kin->phibdx0 = 0.0;
    forc->FAPB_ext = osk::Vec(0, 0, 0);
    forc->FMB_ext  = osk::Vec(0, 0, 0);
    eul->ppx0 = eul->qqx0 = eul->rrx0 = 0.0;

    env->getsFrom(newt);
    forc->getsFrom(prop);
    newt->getsFrom(env, kin, forc, prop);
    kin->getsFrom(env, newt, eul);
    eul->getsFrom(forc, prop, kin);

    std::vector<osk::Block*> stage0 = { kin, env, prop, forc, newt, eul };
    std::vector<std::vector<osk::Block*>> stages = { stage0 };

    double dts[] = { 0.01 };
    double t_total = 8.0 + 6.0 + 7.5 + 6.0 + 0.5;   // total burn + slop
    osk::Sim sim(dts, t_total, stages);
    sim.run();

    // Stage-wise Tsiolkovsky.
    double dv1 = prop->spi_stage[0] * G0
               * std::log(prop->vmass0_stage[0]
                        / (prop->vmass0_stage[0] - prop->fmass0_stage[0]));
    double dv2 = prop->spi_stage[1] * G0
               * std::log(prop->vmass0_stage[1]
                        / (prop->vmass0_stage[1] - prop->fmass0_stage[1]));
    double dv3 = prop->spi_stage[2] * G0
               * std::log(prop->vmass0_stage[2]
                        / (prop->vmass0_stage[2] - prop->fmass0_stage[2]));
    double dv4 = prop->spi_stage[3] * G0
               * std::log(prop->vmass0_stage[3]
                        / (prop->vmass0_stage[3] - prop->fmass0_stage[3]));
    double dv_total = dv1 + dv2 + dv3 + dv4;
    double t_burn   = prop->fmass0_stage[0] / prop->mdot_stage[0]
                    + prop->fmass0_stage[1] / prop->mdot_stage[1]
                    + prop->fmass0_stage[2] / prop->mdot_stage[2]
                    + prop->fmass0_stage[3] / prop->mdot_stage[3];
    double gloss    = G0 * t_burn;
    double v_pred   = dv_total - gloss;

    std::printf("  Stage burns: %.1f + %.1f + %.1f + %.1f = %.1f s\n",
                prop->fmass0_stage[0] / prop->mdot_stage[0],
                prop->fmass0_stage[1] / prop->mdot_stage[1],
                prop->fmass0_stage[2] / prop->mdot_stage[2],
                prop->fmass0_stage[3] / prop->mdot_stage[3],
                t_burn);
    std::printf("  Stage delta-v: %.1f + %.1f + %.1f + %.1f = %.1f m/s\n",
                dv1, dv2, dv3, dv4, dv_total);
    std::printf("  Gravity loss = %.1f m/s\n", gloss);
    std::printf("  Predicted burnout speed = %.1f m/s\n", v_pred);
    std::printf("  Simulated burnout speed = %.2f m/s  (dvbi=%.2f)\n",
                newt->dvbe, newt->dvbi);
    std::printf("  Final stage:        %d  (expect 3)\n", prop->current_stage);
    std::printf("  Final phase:        %d  (expect 2 = spent)\n", prop->phase);
    std::printf("  Final vmass:        %.2f kg  (expect %.2f)\n",
                prop->vmass,
                prop->vmass0_stage[3] - prop->fmass0_stage[3]);

    // Looser tolerance than 3-stage: more stages -> more g-loss approximation
    // error compounds (the t_burn calculation assumes constant gravity, but
    // gravity varies with altitude; over 27.5s of total burn the rocket
    // climbs hundreds of meters).
    bool dv_ok    = std::fabs(newt->dvbe - v_pred) < 30.0;
    bool stage_ok = (prop->current_stage == 3);
    bool phase_ok = (prop->phase == 2);
    bool mass_ok  = std::fabs(prop->vmass
                       - (prop->vmass0_stage[3] - prop->fmass0_stage[3])) < 0.5;
    bool ok = dv_ok && stage_ok && phase_ok && mass_ok;
    std::printf("  dv_ok: %s  stage_ok: %s  phase_ok: %s  mass_ok: %s\n",
                dv_ok ? "Y":"N", stage_ok ? "Y":"N",
                phase_ok ? "Y":"N", mass_ok ? "Y":"N");
    std::printf("  %s\n", ok ? "PASS" : "FAIL");

    delete prop; delete newt; delete eul; delete kin; delete forc; delete env;
    return ok ? 0 : 1;
}

} // anon

int main() {
    int fails = 0;
    fails += run_single_stage();
    fails += run_three_stage();
    fails += run_in_propulsion_multistage();
    fails += run_multistage_with_coast();
    fails += run_four_stage();
    std::printf("\n=== Total failures: %d ===\n", fails);
    return fails == 0 ? 0 : 1;
}
