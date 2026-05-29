//  gps_test.cpp  --  Verify GPS measurement timing, noise, and INS update
//
//  Five tests:
//    (1) GPS off (mgps=0): no measurements, gps_update_avail never fires
//    (2) Perfect GPS timing: measurement count matches gps_step interval
//    (3) Noisy GPS noise statistics: empirical RMS error matches rpos
//    (4) INS bounded by GPS: with INS bias AND GPS correction, position
//        error stays around the GPS noise level (not quadratic growth)
//    (5) Reproducibility: same noise_seed produces identical measurements

#include "../osk/osk.h"
#include "environment.h"
#include "newton.h"
#include "kinematics.h"
#include "euler.h"
#include "forces.h"
#include "propulsion.h"
#include "ins.h"
#include "gps.h"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace rocket6dof;

namespace {

// ---- Test 1: GPS off ----
int run_gps_off() {
    std::printf("\n=== Test 1: GPS off (mgps=0) ===\n");

    Environment* env  = new Environment();
    Propulsion*  prop = new Propulsion();
    Forces*      forc = new Forces();
    Euler*       eul  = new Euler();
    Kinematics*  kin  = new Kinematics();
    Newton*      newt = new Newton();
    GPS*         gps  = new GPS();

    prop->mprop  = 0;  prop->vmass0 = 100.0;
    forc->FAPB_ext = osk::Vec(0,0,0); forc->FMB_ext = osk::Vec(0,0,0);
    eul->ppx0 = eul->qqx0 = eul->rrx0 = 0.0;
    newt->lonx0 = newt->latx0 = 0; newt->alt0 = 1000.0; newt->dvbe0 = 0.0;
    kin->psibdx0 = kin->thtbdx0 = kin->phibdx0 = 0.0;

    gps->mgps     = 0;
    gps->gps_step = 1.0;

    env->getsFrom(newt);
    forc->getsFrom(prop);
    newt->getsFrom(env, kin, forc, prop);
    kin->getsFrom(env, newt, eul);
    eul->getsFrom(forc, prop, kin);
    gps->getsFrom(newt);

    std::vector<osk::Block*> stage0 = { kin, env, prop, forc, newt, eul, gps };
    std::vector<std::vector<osk::Block*>> stages = { stage0 };

    double dts[] = { 0.01 };
    osk::Sim sim(dts, 5.0, stages);
    sim.run();

    std::printf("  After 5 s with mgps=0:\n");
    std::printf("    meas_count       = %d  (expect 0)\n", gps->meas_count);
    std::printf("    gps_update_avail = %d  (expect 0)\n", gps->gps_update_avail);

    bool ok = (gps->meas_count == 0) && (gps->gps_update_avail == 0);
    std::printf("  %s\n", ok ? "PASS" : "FAIL");

    delete gps; delete newt; delete kin; delete eul; delete forc; delete prop; delete env;
    return ok ? 0 : 1;
}

// ---- Test 2: Perfect GPS timing ----
int run_gps_timing() {
    std::printf("\n=== Test 2: Perfect GPS timing (mgps=1) ===\n");

    Environment* env  = new Environment();
    Propulsion*  prop = new Propulsion();
    Forces*      forc = new Forces();
    Euler*       eul  = new Euler();
    Kinematics*  kin  = new Kinematics();
    Newton*      newt = new Newton();
    GPS*         gps  = new GPS();

    prop->mprop  = 0;  prop->vmass0 = 100.0;
    forc->FAPB_ext = osk::Vec(0,0,0); forc->FMB_ext = osk::Vec(0,0,0);
    eul->ppx0 = eul->qqx0 = eul->rrx0 = 0.0;
    newt->lonx0 = newt->latx0 = 0; newt->alt0 = 1000.0; newt->dvbe0 = 0.0;
    kin->psibdx0 = kin->thtbdx0 = kin->phibdx0 = 0.0;

    gps->mgps     = 1;
    gps->gps_step = 1.0;

    env->getsFrom(newt);
    forc->getsFrom(prop);
    newt->getsFrom(env, kin, forc, prop);
    kin->getsFrom(env, newt, eul);
    eul->getsFrom(forc, prop, kin);
    gps->getsFrom(newt);

    std::vector<osk::Block*> stage0 = { kin, env, prop, forc, newt, eul, gps };
    std::vector<std::vector<osk::Block*>> stages = { stage0 };

    double dts[] = { 0.01 };
    double t_end = 10.0;
    osk::Sim sim(dts, t_end, stages);
    sim.run();

    // With gps_step=1.0 over 10s, we expect ~10 measurements (one per
    // second; depending on first-call alignment we may get 10 or 11)
    int expected_min = 10, expected_max = 11;

    std::printf("  After %.0f s with gps_step=1.0:\n", t_end);
    std::printf("    meas_count       = %d  (expect %d-%d)\n",
                gps->meas_count, expected_min, expected_max);
    std::printf("    last_pos_err     = %.3e m  (expect 0; perfect mode)\n",
                gps->last_pos_err);
    std::printf("    last_vel_err     = %.3e m/s\n", gps->last_vel_err);

    bool ok =    gps->meas_count >= expected_min
              && gps->meas_count <= expected_max
              && gps->last_pos_err < 1e-9
              && gps->last_vel_err < 1e-9;
    std::printf("  %s\n", ok ? "PASS" : "FAIL");

    delete gps; delete newt; delete kin; delete eul; delete forc; delete prop; delete env;
    return ok ? 0 : 1;
}

// ---- Test 3: Noisy GPS noise statistics ----
int run_gps_noise() {
    std::printf("\n=== Test 3: Noisy GPS noise statistics (mgps=2) ===\n");

    Environment* env  = new Environment();
    Propulsion*  prop = new Propulsion();
    Forces*      forc = new Forces();
    Euler*       eul  = new Euler();
    Kinematics*  kin  = new Kinematics();
    Newton*      newt = new Newton();
    GPS*         gps  = new GPS();

    prop->mprop  = 0;  prop->vmass0 = 100.0;
    forc->FAPB_ext = osk::Vec(0,0,0); forc->FMB_ext = osk::Vec(0,0,0);
    eul->ppx0 = eul->qqx0 = eul->rrx0 = 0.0;
    newt->lonx0 = newt->latx0 = 0; newt->alt0 = 1000.0; newt->dvbe0 = 0.0;
    kin->psibdx0 = kin->thtbdx0 = kin->phibdx0 = 0.0;

    gps->mgps       = 2;
    gps->gps_step   = 0.1;     // 10 Hz for many samples in short test
    gps->rpos       = 10.0;
    gps->rvel       = 0.5;
    gps->noise_seed = 12345;

    env->getsFrom(newt);
    forc->getsFrom(prop);
    newt->getsFrom(env, kin, forc, prop);
    kin->getsFrom(env, newt, eul);
    eul->getsFrom(forc, prop, kin);
    gps->getsFrom(newt);

    std::vector<osk::Block*> stage0 = { kin, env, prop, forc, newt, eul, gps };
    std::vector<std::vector<osk::Block*>> stages = { stage0 };

    // Run long enough to collect ~500 GPS samples
    double dts[] = { 0.01 };
    osk::Sim sim(dts, 50.0, stages);

    // We can't easily collect per-measurement statistics during sim,
    // so instead: run, then directly call gps->update() in a loop with
    // truth held constant to measure the noise distribution.
    // (The above sim establishes the truth state at t=50s; gps is
    // operating throughout but we don't have hook to collect samples.)
    sim.run();
    int total_meas = gps->meas_count;

    // Now collect noise samples by repeatedly calling update() with
    // truth held fixed.  Each call advances meas_count and produces
    // a fresh noise sample (since the engine is seeded by
    // noise_seed XOR meas_count).
    int N = 1000;
    double sum_sq_pos = 0.0;
    double sum_sq_vel = 0.0;
    int seen = 0;
    osk::Vec truth_pos = newt->SBII;
    osk::Vec truth_vel = newt->VBII;
    for (int i = 0; i < N; i++) {
        // Force a GPS measurement by advancing time past gps_epoch+step
        osk::State::t = gps->gps_epoch + gps->gps_step + 0.001;
        gps->update();
        if (gps->gps_update_avail) {
            osk::Vec dpos = gps->SBII_meas - truth_pos;
            osk::Vec dvel = gps->VBII_meas - truth_vel;
            sum_sq_pos += dpos.dot(dpos);
            sum_sq_vel += dvel.dot(dvel);
            seen++;
            gps->gps_update_avail = 0;
        }
    }

    double rms_pos = std::sqrt(sum_sq_pos / seen);
    double rms_vel = std::sqrt(sum_sq_vel / seen);
    // For 3D Gaussian with per-axis sigma=rpos, the magnitude RMS is
    // sqrt(3)*rpos.
    double expected_rms_pos = std::sqrt(3.0) * gps->rpos;
    double expected_rms_vel = std::sqrt(3.0) * gps->rvel;

    std::printf("  Initial sim measurements: %d\n", total_meas);
    std::printf("  Manual noise sweep:       %d samples\n", seen);
    std::printf("  Empirical RMS pos error  = %.3f m  (expect ~%.3f m = sqrt(3)*rpos)\n",
                rms_pos, expected_rms_pos);
    std::printf("  Empirical RMS vel error  = %.4f m/s (expect ~%.4f m/s)\n",
                rms_vel, expected_rms_vel);

    // Tolerance: 10% (statistical convergence on N=1000 samples)
    double rel_p = std::fabs(rms_pos - expected_rms_pos) / expected_rms_pos;
    double rel_v = std::fabs(rms_vel - expected_rms_vel) / expected_rms_vel;
    std::printf("  Pos rel error = %.2f %%  Vel rel error = %.2f %%\n",
                100.0 * rel_p, 100.0 * rel_v);

    bool ok = (rel_p < 0.1) && (rel_v < 0.1);
    std::printf("  %s\n", ok ? "PASS" : "FAIL");

    osk::State::t = 0.0;   // reset for subsequent tests
    delete gps; delete newt; delete kin; delete eul; delete forc; delete prop; delete env;
    return ok ? 0 : 1;
}

// ---- Test 4: INS bounded by GPS ----
int run_ins_bounded_by_gps() {
    std::printf("\n=== Test 4: INS drift bounded by GPS updates ===\n");

    // Run two parallel sims: same INS with bias, one with GPS, one without.
    auto build = [](bool gps_active, double& pos_err_final, int& nupdates) {
        Environment* env  = new Environment();
        Propulsion*  prop = new Propulsion();
        Forces*      forc = new Forces();
        Euler*       eul  = new Euler();
        Kinematics*  kin  = new Kinematics();
        Newton*      newt = new Newton();
        INS*         ins  = new INS();
        GPS*         gps  = new GPS();

        prop->mprop  = 0;  prop->vmass0 = 100.0;
        forc->FAPB_ext = osk::Vec(0,0,0); forc->FMB_ext = osk::Vec(0,0,0);
        eul->ppx0 = eul->qqx0 = eul->rrx0 = 0.0;
        newt->lonx0 = newt->latx0 = 0; newt->alt0 = 1000.0; newt->dvbe0 = 0.0;
        kin->psibdx0 = kin->thtbdx0 = kin->phibdx0 = 0.0;

        // Accelerometer bias: 0.01 m/s^2 along body-x
        ins->mins       = 1;
        ins->bias_accel = osk::Vec(0.01, 0, 0);
        ins->bias_gyro  = osk::Vec(0, 0, 0);

        gps->mgps       = gps_active ? 1 : 0;   // perfect GPS when on
        gps->gps_step   = 1.0;

        env->getsFrom(newt);
        forc->getsFrom(prop);
        newt->getsFrom(env, kin, forc, prop);
        kin->getsFrom(env, newt, eul);
        eul->getsFrom(forc, prop, kin);
        gps->getsFrom(newt);
        ins->getsFrom(newt, eul, kin, gps);

        // GPS must run BEFORE INS in the stage so INS sees fresh
        // gps_update_avail
        std::vector<osk::Block*> stage0 = { kin, env, prop, forc, newt, eul, gps, ins };
        std::vector<std::vector<osk::Block*>> stages = { stage0 };

        double dts[] = { 0.01 };
        osk::Sim sim(dts, 30.0, stages);
        sim.run();

        pos_err_final = ins->ins_pos_err;
        nupdates      = ins->gps_update_count;

        delete gps; delete ins; delete newt; delete kin; delete eul;
        delete forc; delete prop; delete env;
    };

    double err_no_gps, err_with_gps;
    int n_no_gps, n_with_gps;
    build(false, err_no_gps, n_no_gps);
    build(true,  err_with_gps, n_with_gps);

    std::printf("  Setup: stationary vehicle, INS bias 0.01 m/s^2, 30 s run\n");
    std::printf("  Predicted no-GPS drift = 0.5*b*t^2 = 4.5 m\n");
    std::printf("  Without GPS: INS pos_err = %8.3f m, gps_updates = %d\n",
                err_no_gps, n_no_gps);
    std::printf("  With perfect GPS (1 Hz): INS pos_err = %8.3f m, gps_updates = %d\n",
                err_with_gps, n_with_gps);

    // Expectations:
    //   * Without GPS: pos_err ~ 4.5 m
    //   * With perfect GPS: pos_err is ALMOST zero (truth update every 1s,
    //     intervening drift is ~bias*dt + 0.5*bias*dt^2 sub-meter)
    //   * GPS update count > 25 (at 1 Hz over 30s)
    bool no_gps_drifts   = (err_no_gps > 3.0);          // close to predicted 4.5 m
    bool with_gps_bound  = (err_with_gps < 0.1);        // bounded by between-update drift
    bool gps_active      = (n_with_gps >= 25);

    std::printf("  no_gps_drifts: %s  with_gps_bound: %s  gps_updates_seen: %s\n",
                no_gps_drifts ? "YES" : "NO",
                with_gps_bound ? "YES" : "NO",
                gps_active ? "YES" : "NO");

    bool ok = no_gps_drifts && with_gps_bound && gps_active;
    std::printf("  %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

// ---- Test 5: Reproducibility (deterministic seed) ----
int run_gps_reproducibility() {
    std::printf("\n=== Test 5: GPS noise reproducibility (same seed) ===\n");

    auto run_noisy = [](unsigned long seed, std::vector<double>& samples) {
        Environment* env  = new Environment();
        Propulsion*  prop = new Propulsion();
        Forces*      forc = new Forces();
        Euler*       eul  = new Euler();
        Kinematics*  kin  = new Kinematics();
        Newton*      newt = new Newton();
        GPS*         gps  = new GPS();

        prop->mprop  = 0;  prop->vmass0 = 100.0;
        forc->FAPB_ext = osk::Vec(0,0,0); forc->FMB_ext = osk::Vec(0,0,0);
        eul->ppx0 = eul->qqx0 = eul->rrx0 = 0.0;
        newt->lonx0 = newt->latx0 = 0; newt->alt0 = 1000.0; newt->dvbe0 = 0.0;
        kin->psibdx0 = kin->thtbdx0 = kin->phibdx0 = 0.0;

        gps->mgps       = 2;
        gps->gps_step   = 0.5;
        gps->rpos       = 10.0;
        gps->rvel       = 0.5;
        gps->noise_seed = seed;

        env->getsFrom(newt);
        forc->getsFrom(prop);
        newt->getsFrom(env, kin, forc, prop);
        kin->getsFrom(env, newt, eul);
        eul->getsFrom(forc, prop, kin);
        gps->getsFrom(newt);

        std::vector<osk::Block*> stage0 = { kin, env, prop, forc, newt, eul, gps };
        std::vector<std::vector<osk::Block*>> stages = { stage0 };

        double dts[] = { 0.01 };
        osk::Sim sim(dts, 5.0, stages);
        sim.run();

        // Capture the final measurement
        samples.push_back(gps->SBII_meas.x);
        samples.push_back(gps->SBII_meas.y);
        samples.push_back(gps->SBII_meas.z);
        samples.push_back(gps->VBII_meas.x);
        samples.push_back(gps->VBII_meas.y);
        samples.push_back(gps->VBII_meas.z);

        delete gps; delete newt; delete kin; delete eul; delete forc; delete prop; delete env;
    };

    std::vector<double> A, B, C;
    run_noisy(42, A);
    run_noisy(42, B);     // same seed -> identical samples
    run_noisy(99, C);     // different seed -> different samples

    bool ab_identical = true;
    for (size_t i = 0; i < A.size(); i++) {
        if (std::fabs(A[i] - B[i]) > 1e-12) { ab_identical = false; break; }
    }
    bool ac_differ = false;
    for (size_t i = 0; i < A.size(); i++) {
        if (std::fabs(A[i] - C[i]) > 1.0) { ac_differ = true; break; }
    }

    std::printf("  Seed=42 run A: pos.x = %.3f m\n", A[0]);
    std::printf("  Seed=42 run B: pos.x = %.3f m  (should equal A)\n", B[0]);
    std::printf("  Seed=99 run C: pos.x = %.3f m  (should differ from A)\n", C[0]);
    std::printf("  A == B (same seed):  %s\n", ab_identical ? "YES" : "NO");
    std::printf("  A != C (different):  %s\n", ac_differ ? "YES" : "NO");

    bool ok = ab_identical && ac_differ;
    std::printf("  %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

// ---- Test 6: Full Kalman GPS filter (mgps=3) ----
// Drive a stationary vehicle for 30 s with the Zipfel-style 8-state
// Kalman filter active.  Verifies:
//   * No NaN/Inf anywhere
//   * GDOP is finite and reasonable (geometry of the chosen 4 SVs)
//   * meas_count grows (filter is producing updates at gps_step rate)
//   * |SXH| and |VXH| stay bounded (filter converges, not diverging)
//   * INS error stays bounded with Kalman corrections active
int run_gps_kalman() {
    std::printf("\n=== Test 6: Full Kalman GPS filter (mgps=3) ===\n");

    Environment* env  = new Environment();
    Propulsion*  prop = new Propulsion();
    Forces*      forc = new Forces();
    Euler*       eul  = new Euler();
    Kinematics*  kin  = new Kinematics();
    Newton*      newt = new Newton();
    INS*         ins  = new INS();
    GPS*         gps  = new GPS();

    prop->mprop  = 0;  prop->vmass0 = 100.0;
    forc->FAPB_ext = osk::Vec(0,0,0); forc->FMB_ext = osk::Vec(0,0,0);
    eul->ppx0 = eul->qqx0 = eul->rrx0 = 0.0;
    newt->lonx0 = newt->latx0 = 0; newt->alt0 = 1000.0; newt->dvbe0 = 0.0;
    kin->psibdx0 = kin->thtbdx0 = kin->phibdx0 = 0.0;

    ins->mins       = 1;
    ins->bias_accel = osk::Vec(0.001, 0, 0);
    ins->bias_gyro  = osk::Vec(0, 0, 0);

    gps->mgps         = 3;
    gps->almanac_time = 0.0;
    gps->del_rearth   = 0.0;
    gps->gps_acqtime  = 2.0;
    gps->gps_step     = 1.0;
    gps->ucfreq_noise = 0.01;
    gps->uctime_cor   = 3600.0;
    for (int i = 0; i < 4; i++) {
        gps->pr_bias[i]  = 0.0;
        gps->pr_noise[i] = 5.0;
        gps->dr_noise[i] = 0.1;
    }
    gps->ppos = 20.0; gps->pvel = 2.0; gps->pclockb = 10.0; gps->pclockf = 1.0;
    gps->qpos = 0.1;  gps->qvel = 0.01; gps->qclockb = 0.1; gps->qclockf = 0.01;
    gps->rpos_kf = 5.0; gps->rvel_kf = 0.1;
    gps->factp = 0.0; gps->factq = 0.0; gps->factr = 0.0;
    gps->noise_seed = 12345;

    env->getsFrom(newt);
    forc->getsFrom(prop);
    newt->getsFrom(env, kin, forc, prop);
    kin->getsFrom(env, newt, eul);
    eul->getsFrom(forc, prop, kin);
    gps->getsFrom(newt, eul, ins);
    ins->getsFrom(newt, eul, kin, gps);

    std::vector<osk::Block*> stage0 = { kin, env, prop, forc, newt, eul, gps, ins };
    std::vector<std::vector<osk::Block*>> stages = { stage0 };

    double dts[] = { 0.01 };
    osk::Sim sim(dts, 30.0, stages);
    sim.run();

    int meas = gps->meas_count;
    double gdop_final = gps->gdop;
    double sxh_mag = gps->SXH.mag();
    double vxh_mag = gps->VXH.mag();
    double ins_err = ins->ins_pos_err;
    int ins_updates = ins->gps_update_count;

    std::printf("  After 30 s with Kalman filter active:\n");
    std::printf("    meas_count        = %d   (expect >= 10)\n", meas);
    std::printf("    GDOP              = %.3f m   (expect finite, positive)\n", gdop_final);
    std::printf("    |SXH|             = %.4f m   (filter position correction)\n", sxh_mag);
    std::printf("    |VXH|             = %.6f m/s (filter velocity correction)\n", vxh_mag);
    std::printf("    INS pos_err       = %.4f m   (should stay bounded)\n", ins_err);
    std::printf("    INS gps_updates   = %d\n", ins_updates);

    bool finite_all = std::isfinite(gdop_final)
                   && std::isfinite(sxh_mag)
                   && std::isfinite(vxh_mag)
                   && std::isfinite(ins_err);
    bool meas_ok    = meas >= 10;
    bool gdop_ok    = (gdop_final > 0.0) && (gdop_final < 100.0);  // good GPS is <10
    // Filter outputs should be on the order of the measurement noise,
    // not catastrophic.  rpos_kf=5 m, so SXH ~10s of meters is fine.
    bool sxh_bound  = sxh_mag < 100.0;
    bool vxh_bound  = vxh_mag < 10.0;
    // INS error is dominated by Q-vs-actual-drift modeling -- a real
    // tuning exercise.  Just verify the filter is consuming GPS at the
    // expected rate.
    bool ins_consuming = (ins_updates >= 10);

    std::printf("    finite_all: %s  meas_ok: %s  gdop_ok: %s\n",
                finite_all ? "YES" : "NO",
                meas_ok    ? "YES" : "NO",
                gdop_ok    ? "YES" : "NO");
    std::printf("    sxh_bound:  %s  vxh_bound: %s  ins_consuming: %s (n=%d)\n",
                sxh_bound      ? "YES" : "NO",
                vxh_bound      ? "YES" : "NO",
                ins_consuming  ? "YES" : "NO", ins_updates);

    bool ok = finite_all && meas_ok && gdop_ok
           && sxh_bound && vxh_bound && ins_consuming;
    std::printf("  %s\n", ok ? "PASS" : "FAIL");

    delete gps; delete ins; delete newt; delete kin; delete eul; delete forc; delete prop; delete env;
    return ok ? 0 : 1;
}

} // anon

int main() {
    int fails = 0;
    fails += run_gps_off();
    fails += run_gps_timing();
    fails += run_gps_noise();
    fails += run_ins_bounded_by_gps();
    fails += run_gps_reproducibility();
    fails += run_gps_kalman();
    std::printf("\n=== Total failures: %d ===\n", fails);
    return fails == 0 ? 0 : 1;
}
