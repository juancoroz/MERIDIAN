//  startrack_test.cpp -- Verify Startrack timing, noise, and INS update
//
//  Five tests, parallel to gps_test:
//    (1) mstar=0: no measurements, no INS update
//    (2) mstar=1 timing: perfect-mode measurements at startrack_step rate
//    (3) mstar=2 noise statistics: empirical RMS of URIC matches
//        sqrt(3)*tilt_noise
//    (4) INS attitude bounded: gyro bias produces drift; startrack
//        corrects to a small residual vs unbounded growth
//    (5) mstar=3 full triad: 3 stars selected, volume > 0.1

#include "../osk/osk.h"
#include "environment.h"
#include "newton.h"
#include "kinematics.h"
#include "euler.h"
#include "forces.h"
#include "propulsion.h"
#include "ins.h"
#include "startrack.h"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace rocket6dof;

namespace {

struct VehicleBundle {
    Environment* env;
    Propulsion*  prop;
    Forces*      forc;
    Euler*       eul;
    Kinematics*  kin;
    Newton*      newt;
    INS*         ins;
    Startrack*   star;

    VehicleBundle() {
        env  = new Environment();
        prop = new Propulsion();
        forc = new Forces();
        eul  = new Euler();
        kin  = new Kinematics();
        newt = new Newton();
        ins  = new INS();
        star = new Startrack();

        prop->mprop  = 0;  prop->vmass0 = 100.0;
        forc->FAPB_ext = osk::Vec(0,0,0); forc->FMB_ext = osk::Vec(0,0,0);
        eul->ppx0 = eul->qqx0 = eul->rrx0 = 0.0;
        newt->lonx0 = 0; newt->latx0 = 0; newt->alt0 = 500000.0;
        newt->dvbe0 = 0.0;
        kin->psibdx0 = 0; kin->thtbdx0 = 0; kin->phibdx0 = 0;
    }
    void wire() {
        env->getsFrom(newt);
        forc->getsFrom(prop);
        newt->getsFrom(env, kin, forc, prop);
        kin->getsFrom(env, newt, eul);
        eul->getsFrom(forc, prop, kin);
        star->getsFrom(newt, kin, ins);
        ins->getsFrom(newt, eul, kin, nullptr, star);
    }
    std::vector<osk::Block*> stage() {
        return { kin, env, prop, forc, newt, eul, star, ins };
    }
    ~VehicleBundle() {
        delete star; delete ins; delete newt; delete kin; delete eul;
        delete forc; delete prop; delete env;
    }
};

int run_star_off() {
    std::printf("\n=== Test 1: Startrack off (mstar=0) ===\n");
    VehicleBundle vb;
    vb.star->mstar = 0;
    vb.wire();
    auto stages = std::vector<std::vector<osk::Block*>>{ vb.stage() };
    double dts[] = { 0.01 };
    osk::Sim sim(dts, 5.0, stages);
    sim.run();
    std::printf("  meas_count=%d (expect 0)  star_update_avail=%d (expect 0)\n",
                vb.star->meas_count, vb.star->star_update_avail);
    bool ok = (vb.star->meas_count == 0) && (vb.star->star_update_avail == 0);
    std::printf("  %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

int run_star_timing() {
    std::printf("\n=== Test 2: Perfect Startrack timing (mstar=1) ===\n");
    VehicleBundle vb;
    vb.star->mstar          = 1;
    vb.star->startrack_step = 1.0;
    vb.star->star_acqtime   = 0.5;
    vb.ins->mins = 0;
    vb.wire();
    auto stages = std::vector<std::vector<osk::Block*>>{ vb.stage() };
    double dts[] = { 0.01 };
    osk::Sim sim(dts, 10.0, stages);
    sim.run();
    int meas = vb.star->meas_count;
    std::printf("  meas_count = %d  (expect 9-11)\n", meas);
    bool ok = (meas >= 9 && meas <= 11);
    std::printf("  %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

int run_star_noise() {
    std::printf("\n=== Test 3: Noisy Startrack noise stats (mstar=2) ===\n");
    VehicleBundle vb;
    vb.star->mstar          = 2;
    vb.star->startrack_step = 0.1;
    vb.star->star_acqtime   = 0.0;
    vb.star->tilt_noise     = 1.0e-4;
    vb.star->noise_seed     = 54321;
    vb.ins->mins = 0;
    vb.wire();
    auto stages = std::vector<std::vector<osk::Block*>>{ vb.stage() };
    double dts[] = { 0.01 };
    osk::Sim sim(dts, 10.0, stages);
    sim.run();

    // Manual sweep for statistics
    int N = 1000;
    double sum_sq = 0.0;
    int seen = 0;
    for (int i = 0; i < N; i++) {
        osk::State::t = vb.star->starfix_epoch + vb.star->startrack_step + 0.001;
        vb.star->update();
        if (vb.star->star_update_avail) {
            sum_sq += vb.star->URIC.dot(vb.star->URIC);
            seen++;
            vb.star->star_update_avail = 0;
        }
    }
    double rms = std::sqrt(sum_sq / seen);
    double expected = std::sqrt(3.0) * vb.star->tilt_noise;
    std::printf("  Manual sweep: %d samples\n", seen);
    std::printf("  RMS |URIC| = %.5e rad  (expect ~%.5e)\n", rms, expected);
    double rel = std::fabs(rms - expected) / expected;
    std::printf("  Relative error = %.2f %%\n", 100.0 * rel);
    osk::State::t = 0.0;
    bool ok = (rel < 0.1);
    std::printf("  %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

// TESTS_4_5_MARKER

int run_ins_bounded_by_star() {
    std::printf("\n=== Test 4: INS attitude bounded by Startrack ===\n");

    auto build = [](bool star_active, double& att_err_final, int& nupdates) {
        VehicleBundle vb;
        vb.ins->mins       = 1;
        vb.ins->bias_accel = osk::Vec(0, 0, 0);
        vb.ins->bias_gyro  = osk::Vec(1.0e-4, 0, 0);
        vb.star->mstar          = star_active ? 1 : 0;
        vb.star->startrack_step = 1.0;
        vb.star->star_acqtime   = 0.0;
        vb.wire();
        auto stages = std::vector<std::vector<osk::Block*>>{ vb.stage() };
        double dts[] = { 0.01 };
        osk::Sim sim(dts, 30.0, stages);
        sim.run();
        att_err_final = vb.ins->ins_att_err;
        nupdates      = vb.ins->startrack_update_count;
    };

    double err_off, err_on;
    int n_off, n_on;
    build(false, err_off, n_off);
    build(true,  err_on,  n_on);

    std::printf("  Setup: gyro bias 1e-4 rad/s along x, 30 s run\n");
    std::printf("  Predicted no-star drift = b*t = 0.003 rad\n");
    std::printf("  Without star: ins_att_err = %.6e rad, star_updates = %d\n",
                err_off, n_off);
    std::printf("  With perfect star (1 Hz): ins_att_err = %.6e rad, star_updates = %d\n",
                err_on, n_on);

    bool off_drifts = (err_off > 0.001);
    bool on_bounded = (err_on  < 0.001);
    bool updates_ok = (n_on >= 25);

    std::printf("  off_drifts: %s  on_bounded: %s  updates_seen: %s (n=%d)\n",
                off_drifts ? "YES" : "NO",
                on_bounded ? "YES" : "NO",
                updates_ok ? "YES" : "NO", n_on);

    bool ok = off_drifts && on_bounded && updates_ok;
    std::printf("  %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

int run_star_full_triad() {
    std::printf("\n=== Test 5: Full triad-based Startrack (mstar=3) ===\n");
    VehicleBundle vb;
    vb.star->mstar          = 3;
    vb.star->startrack_step = 1.0;
    vb.star->star_acqtime   = 2.0;
    vb.star->startrack_alt  = 100000.0;
    vb.star->star_el_min    = 20.0;
    for (int i = 0; i < 3; i++) {
        vb.star->az_bias[i]  = 0.0;
        vb.star->az_noise[i] = 1.0e-4;
        vb.star->el_bias[i]  = 0.0;
        vb.star->el_noise[i] = 1.0e-4;
    }
    vb.star->noise_seed     = 12345;
    vb.ins->mins = 1;
    vb.ins->bias_accel = osk::Vec(0, 0, 0);
    vb.ins->bias_gyro  = osk::Vec(0, 0, 0);
    vb.wire();

    auto stages = std::vector<std::vector<osk::Block*>>{ vb.stage() };
    double dts[] = { 0.01 };
    osk::Sim sim(dts, 10.0, stages);
    sim.run();

    std::printf("  After 10 s with mstar=3:\n");
    std::printf("    meas_count        = %d (expect >= 5)\n", vb.star->meas_count);
    std::printf("    triad slots       = {%d, %d, %d}\n",
                vb.star->triad[0], vb.star->triad[1], vb.star->triad[2]);
    std::printf("    parallelepiped volume = %.4f (expect > 0.1)\n",
                vb.star->star_volume);
    std::printf("    last |URIC|       = %.2e rad\n", vb.star->last_tilt_mag);
    std::printf("    ins_att_err       = %.2e rad\n", vb.ins->ins_att_err);

    bool meas_ok    = vb.star->meas_count >= 5;
    bool slots_ok   =  vb.star->triad[0] >= 1 && vb.star->triad[0] <= 25
                    && vb.star->triad[1] >= 1 && vb.star->triad[1] <= 25
                    && vb.star->triad[2] >= 1 && vb.star->triad[2] <= 25
                    && vb.star->triad[0] != vb.star->triad[1]
                    && vb.star->triad[0] != vb.star->triad[2]
                    && vb.star->triad[1] != vb.star->triad[2];
    bool volume_ok  = (vb.star->star_volume > 0.1);
    bool finite_ok  = std::isfinite(vb.star->last_tilt_mag)
                   && std::isfinite(vb.ins->ins_att_err);
    bool tilt_small = vb.star->last_tilt_mag < 0.01;

    std::printf("    meas_ok:%s slots_ok:%s volume_ok:%s finite_ok:%s tilt_small:%s\n",
                meas_ok    ? "Y":"N", slots_ok ? "Y":"N",
                volume_ok  ? "Y":"N", finite_ok? "Y":"N",
                tilt_small ? "Y":"N");

    bool ok = meas_ok && slots_ok && volume_ok && finite_ok && tilt_small;
    std::printf("  %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

} // anon

// MAIN_MARKER

int main() {
    int fails = 0;
    fails += run_star_off();
    fails += run_star_timing();
    fails += run_star_noise();
    fails += run_ins_bounded_by_star();
    fails += run_star_full_triad();
    std::printf("\n=== Total failures: %d ===\n", fails);
    return fails == 0 ? 0 : 1;
}
