//  ins_test.cpp  --  Verify INS outputs and drift behavior
//
//  Three tests:
//    (1) Ideal INS (mins=0): outputs identically match truth.  No drift.
//    (2) Biased INS (mins=1): with a constant accelerometer bias, INS
//        position drifts at the predicted rate ~ 0.5 * bias * t^2 over
//        a no-motion trajectory.
//    (3) Control reads INS: Control's autopilot output matches when
//        wired through ideal INS vs read directly from Newton/Euler.

#include "../osk/osk.h"
#include "environment.h"
#include "newton.h"
#include "kinematics.h"
#include "euler.h"
#include "forces.h"
#include "propulsion.h"
#include "aerodynamics.h"
#include "tvc.h"
#include "control.h"
#include "ins.h"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace rocket6dof;

namespace {

// ---- Test 1: ideal INS tracks truth ----
int run_ideal_ins() {
    std::printf("\n=== Test 1: Ideal INS (mins=0) tracks truth ===\n");

    Environment* env   = new Environment();
    Propulsion*  prop  = new Propulsion();
    Forces*      forc  = new Forces();
    Euler*       eul   = new Euler();
    Kinematics*  kin   = new Kinematics();
    Newton*      newt  = new Newton();
    INS*         ins   = new INS();

    // Inert vehicle dropping from 1000 m
    prop->mprop  = 0;
    prop->vmass0 = 100.0;
    forc->FAPB_ext = osk::Vec(0, 0, 0);
    forc->FMB_ext  = osk::Vec(0, 0, 0);
    eul->ppx0 = eul->qqx0 = eul->rrx0 = 0.0;

    newt->lonx0  = 0.0; newt->latx0 = 0.0; newt->alt0 = 1000.0;
    newt->dvbe0  = 0.0;
    kin->psibdx0 = 0.0; kin->thtbdx0 = 0.0; kin->phibdx0 = 0.0;

    // INS in ideal mode
    ins->mins       = 0;
    ins->bias_accel = osk::Vec(0, 0, 0);
    ins->bias_gyro  = osk::Vec(0, 0, 0);

    env->getsFrom(newt);
    forc->getsFrom(prop);
    newt->getsFrom(env, kin, forc, prop);
    kin->getsFrom(env, newt, eul);
    eul->getsFrom(forc, prop, kin);
    ins->getsFrom(newt, eul, kin);

    // Stage: physics first, then INS (which reads physics)
    std::vector<osk::Block*> stage0 = { kin, env, prop, forc,
                                         newt, eul, ins };
    std::vector<std::vector<osk::Block*>> stages = { stage0 };

    double dts[] = { 0.01 };
    osk::Sim sim(dts, 5.0, stages);
    sim.run();

    // Verify INS outputs equal truth
    osk::Vec pos_diff = ins->SBIIC - newt->SBII;
    osk::Vec vel_diff = ins->VBIIC - newt->VBII;
    osk::Vec rate_diff = ins->WBICB - eul->WBIB;
    osk::Vec spec_diff = ins->FSPCB - newt->FSPB;

    std::printf("  At t=5 s (ideal INS, should be zero drift):\n");
    std::printf("    |SBIIC - SBII| = %.3e m\n", pos_diff.mag());
    std::printf("    |VBIIC - VBII| = %.3e m/s\n", vel_diff.mag());
    std::printf("    |WBICB - WBIB| = %.3e rad/s\n", rate_diff.mag());
    std::printf("    |FSPCB - FSPB| = %.3e m/s^2\n", spec_diff.mag());
    std::printf("    ins_pos_err    = %.3e m  (diagnostic)\n", ins->ins_pos_err);
    std::printf("    ins_vel_err    = %.3e m/s\n", ins->ins_vel_err);

    bool ok =    pos_diff.mag()  < 1e-9
              && vel_diff.mag()  < 1e-9
              && rate_diff.mag() < 1e-12
              && spec_diff.mag() < 1e-12;

    std::printf("  %s\n", ok ? "PASS" : "FAIL");

    delete ins; delete newt; delete kin; delete eul; delete forc; delete prop; delete env;
    return ok ? 0 : 1;
}

// ---- Test 2: biased INS drifts at predicted rate ----
int run_biased_ins() {
    std::printf("\n=== Test 2: Biased INS (mins=1) drift behavior ===\n");

    Environment* env   = new Environment();
    Propulsion*  prop  = new Propulsion();
    Forces*      forc  = new Forces();
    Euler*       eul   = new Euler();
    Kinematics*  kin   = new Kinematics();
    Newton*      newt  = new Newton();
    INS*         ins   = new INS();

    // Stationary vehicle on the pad -- no motion, no forces, just bias drift
    prop->mprop  = 0;
    prop->vmass0 = 100.0;
    forc->FAPB_ext = osk::Vec(0, 0, 0);
    forc->FMB_ext  = osk::Vec(0, 0, 0);
    eul->ppx0 = eul->qqx0 = eul->rrx0 = 0.0;

    newt->lonx0  = 0.0; newt->latx0 = 0.0; newt->alt0 = 1000.0;
    newt->dvbe0  = 0.0;
    kin->psibdx0 = 0.0; kin->thtbdx0 = 0.0; kin->phibdx0 = 0.0;

    // 1 m/s^2 accel bias along body-x (= +up in inertial at this point)
    double accel_bias = 0.01;    // 0.01 m/s^2 = 1 milli-g
    ins->mins       = 1;
    ins->bias_accel = osk::Vec(accel_bias, 0, 0);
    ins->bias_gyro  = osk::Vec(0, 0, 0);

    env->getsFrom(newt);
    forc->getsFrom(prop);
    newt->getsFrom(env, kin, forc, prop);
    kin->getsFrom(env, newt, eul);
    eul->getsFrom(forc, prop, kin);
    ins->getsFrom(newt, eul, kin);

    std::vector<osk::Block*> stage0 = { kin, env, prop, forc,
                                         newt, eul, ins };
    std::vector<std::vector<osk::Block*>> stages = { stage0 };

    double t_end = 10.0;
    double dts[] = { 0.01 };
    osk::Sim sim(dts, t_end, stages);
    sim.run();

    // Predicted drift: with constant accel bias b over time t,
    // velocity error grows as b*t, position error as 0.5*b*t^2
    double v_err_predicted = accel_bias * t_end;
    double p_err_predicted = 0.5 * accel_bias * t_end * t_end;

    std::printf("  Setup: vehicle stationary (no truth motion), bias=%.3f m/s^2 on body-x\n",
                accel_bias);
    std::printf("  After %.0f s:\n", t_end);
    std::printf("    Predicted vel drift = b*t   = %.4f m/s\n", v_err_predicted);
    std::printf("    INS    vel drift           = %.4f m/s  (ins_vel_err)\n",
                ins->ins_vel_err);
    std::printf("    Predicted pos drift = 0.5*b*t^2 = %.4f m\n", p_err_predicted);
    std::printf("    INS    pos drift           = %.4f m  (ins_pos_err)\n",
                ins->ins_pos_err);

    // The predicted formula assumes constant bias direction in body frame,
    // which equals constant in inertial only if the body isn't rotating.
    // Vehicle is stationary -> no rotation -> formula exact.
    double v_err_rel = std::fabs(ins->ins_vel_err - v_err_predicted) / v_err_predicted;
    double p_err_rel = std::fabs(ins->ins_pos_err - p_err_predicted) / p_err_predicted;
    std::printf("    vel drift rel error = %.3f %%\n", 100.0 * v_err_rel);
    std::printf("    pos drift rel error = %.3f %%\n", 100.0 * p_err_rel);

    bool ok =    v_err_rel < 0.02     // < 2% (RK4 integration of constant)
              && p_err_rel < 0.02;

    std::printf("  %s\n", ok ? "PASS" : "FAIL");

    delete ins; delete newt; delete kin; delete eul; delete forc; delete prop; delete env;
    return ok ? 0 : 1;
}

// ---- Test 3: gyro bias produces attitude drift at predicted rate ----
// A constant gyro bias b_g (rad/s) on body axis n causes the INS's
// integrated attitude to rotate about that axis at rate b_g.  Over time
// t, the angular displacement is b_g * t.  We extract the angle between
// truth TBI and computed TBIC and verify it matches.
int run_gyro_drift() {
    std::printf("\n=== Test 3: Gyro bias attitude drift ===\n");

    Environment* env  = new Environment();
    Propulsion*  prop = new Propulsion();
    Forces*      forc = new Forces();
    Euler*       eul  = new Euler();
    Kinematics*  kin  = new Kinematics();
    Newton*      newt = new Newton();
    INS*         ins  = new INS();

    prop->mprop  = 0;
    prop->vmass0 = 100.0;
    forc->FAPB_ext = osk::Vec(0, 0, 0);
    forc->FMB_ext  = osk::Vec(0, 0, 0);
    eul->ppx0 = eul->qqx0 = eul->rrx0 = 0.0;

    newt->lonx0  = 0.0; newt->latx0 = 0.0; newt->alt0 = 1000.0;
    newt->dvbe0  = 0.0;
    kin->psibdx0 = 0.0; kin->thtbdx0 = 0.0; kin->phibdx0 = 0.0;

    // 0.001 rad/s gyro bias about body-z (yaw axis)
    double gyro_bias = 0.001;
    ins->mins       = 1;
    ins->bias_accel = osk::Vec(0, 0, 0);
    ins->bias_gyro  = osk::Vec(0, 0, gyro_bias);

    env->getsFrom(newt);
    forc->getsFrom(prop);
    newt->getsFrom(env, kin, forc, prop);
    kin->getsFrom(env, newt, eul);
    eul->getsFrom(forc, prop, kin);
    ins->getsFrom(newt, eul, kin);

    std::vector<osk::Block*> stage0 = { kin, env, prop, forc,
                                         newt, eul, ins };
    std::vector<std::vector<osk::Block*>> stages = { stage0 };

    double t_end = 10.0;
    double dts[] = { 0.01 };
    osk::Sim sim(dts, t_end, stages);
    sim.run();

    // Truth TBI hasn't rotated (no torque); INS TBIC should have rotated
    // by gyro_bias * t_end radians about body-z.
    osk::Mat T_truth = kin->TBI;
    osk::Mat T_ins   = ins->TBIC;

    // Relative rotation: T_truth^T * T_ins.  For small angle, this is
    // approximately I + skew(angle_vec), so trace approximately
    // 3 - angle_mag^2 / 2 isn't quite right -- the correct relation
    // for rotation matrices is trace = 1 + 2*cos(angle), so
    // angle = acos((trace-1)/2).
    osk::Mat R = T_truth.transpose() * T_ins;
    double tr = R[0][0] + R[1][1] + R[2][2];
    if (tr > 3.0)  tr = 3.0;
    if (tr < -1.0) tr = -1.0;
    double angle_actual = std::acos((tr - 1.0) / 2.0);

    double angle_predicted = gyro_bias * t_end;

    std::printf("  Setup: stationary, gyro bias = %.4f rad/s about body-z\n", gyro_bias);
    std::printf("  After %.0f s:\n", t_end);
    std::printf("    Predicted attitude drift = b*t = %.4f rad (%.2f deg)\n",
                angle_predicted, angle_predicted * 180.0 / osk::PI);
    std::printf("    INS computed drift      = %.4f rad (%.2f deg)\n",
                angle_actual, angle_actual * 180.0 / osk::PI);

    double rel_err = std::fabs(angle_actual - angle_predicted) / angle_predicted;
    std::printf("    rel error = %.3f %%\n", 100.0 * rel_err);

    bool ok = rel_err < 0.02;

    std::printf("  %s\n", ok ? "PASS" : "FAIL");

    delete ins; delete newt; delete kin; delete eul; delete forc; delete prop; delete env;
    return ok ? 0 : 1;
}

// ---- Test 4: INS tracks a moving vehicle exactly in ideal mode ----
// With mins=0 the INS outputs equal truth.  But we also want to verify
// the INTEGRATOR STATES (sxc, vxc, t00c, etc.) track truth -- because
// when a future test flips mins to 1, those states are the starting
// point and need to be consistent.  We do a vehicle drop, then check
// that all 15 integrator states equal the corresponding truth values
// after several seconds.
int run_moving_vehicle_tracking() {
    std::printf("\n=== Test 4: INS integrator states track truth (mins=0, moving) ===\n");

    Environment* env  = new Environment();
    Propulsion*  prop = new Propulsion();
    Forces*      forc = new Forces();
    Euler*       eul  = new Euler();
    Kinematics*  kin  = new Kinematics();
    Newton*      newt = new Newton();
    INS*         ins  = new INS();

    prop->mprop  = 0;
    prop->vmass0 = 100.0;
    forc->FAPB_ext = osk::Vec(0, 0, 0);
    forc->FMB_ext  = osk::Vec(0, 0, 0);
    eul->ppx0 = 0.0;  eul->qqx0 = 5.0;  eul->rrx0 = 0.0;   // spinning

    newt->lonx0  = 0.0; newt->latx0 = 0.0; newt->alt0 = 5000.0;
    newt->dvbe0  = 0.0;
    kin->psibdx0 = 0.0; kin->thtbdx0 = 0.0; kin->phibdx0 = 0.0;

    ins->mins = 0;     // ideal

    env->getsFrom(newt);
    forc->getsFrom(prop);
    newt->getsFrom(env, kin, forc, prop);
    kin->getsFrom(env, newt, eul);
    eul->getsFrom(forc, prop, kin);
    ins->getsFrom(newt, eul, kin);

    std::vector<osk::Block*> stage0 = { kin, env, prop, forc,
                                         newt, eul, ins };
    std::vector<std::vector<osk::Block*>> stages = { stage0 };

    double dts[] = { 0.01 };
    osk::Sim sim(dts, 5.0, stages);
    sim.run();

    // Check that INS integrator scalars match truth exactly
    double dpos = std::fabs(ins->sxc - newt->SBII.x)
                + std::fabs(ins->syc - newt->SBII.y)
                + std::fabs(ins->szc - newt->SBII.z);
    double dvel = std::fabs(ins->vxc - newt->VBII.x)
                + std::fabs(ins->vyc - newt->VBII.y)
                + std::fabs(ins->vzc - newt->VBII.z);
    double dtbi = std::fabs(ins->t00c - kin->TBI[0][0])
                + std::fabs(ins->t01c - kin->TBI[0][1])
                + std::fabs(ins->t02c - kin->TBI[0][2])
                + std::fabs(ins->t11c - kin->TBI[1][1])
                + std::fabs(ins->t22c - kin->TBI[2][2]);

    std::printf("  Vehicle dropped from 5 km, spinning at 5 deg/s pitch.\n");
    std::printf("  After 5 s (ideal INS):\n");
    std::printf("    L1(SBIIC_state - SBII_truth) = %.3e m\n", dpos);
    std::printf("    L1(VBIIC_state - VBII_truth) = %.3e m/s\n", dvel);
    std::printf("    L1(TBIC_state  - TBI_truth)  = %.3e (5 components)\n", dtbi);
    std::printf("    altitude truth = %.1f m, INS  = %.1f m\n",
                newt->alt, ins->dbic - 6378137.0);

    bool ok =    dpos < 1e-9
              && dvel < 1e-9
              && dtbi < 1e-9;

    std::printf("  %s\n", ok ? "PASS" : "FAIL");

    delete ins; delete newt; delete kin; delete eul; delete forc; delete prop; delete env;
    return ok ? 0 : 1;
}

// ---- Test 5: Control wired through ideal INS gives same result as direct ----
int run_control_via_ins() {
    std::printf("\n=== Test 5: Control sees same sensors through ideal INS ===\n");

    auto build = [](bool use_ins) -> double {
        Environment*  env  = new Environment();
        Propulsion*   prop = new Propulsion();
        Aerodynamics* aero = new Aerodynamics();
        TVC*          tvc  = new TVC();
        Control*      ctrl = new Control();
        Forces*       forc = new Forces();
        Euler*        eul  = new Euler();
        Kinematics*   kin  = new Kinematics();
        Newton*       newt = new Newton();
        INS*          ins  = new INS();

        prop->mprop          = 3;
        prop->vmass0         = 1000.0;
        prop->fmass0         = 600.0;
        prop->spi            = 300.0;
        prop->fuel_flow_rate = 50.0;
        prop->moi_roll_0     = 10.0;  prop->moi_roll_1  = 5.0;
        prop->moi_trans_0    = 500.0; prop->moi_trans_1 = 250.0;
        prop->xcg_0          = 3.0;   prop->xcg_1       = 2.5;

        aero->maero    = 13;
        aero->refa     = 1.0;
        aero->refd     = 1.0;
        aero->xcg_ref  = 3.0;
        aero->aero_file = "aero.txt";
        aero->tag_ca0   = "ca0slv3_vs_mach";
        aero->tag_caa   = "caaslv3_vs_mach";
        aero->tag_ca0b  = "ca0bslv3_vs_mach";
        aero->tag_cn0   = "cn0slv3_vs_mach_alpha";
        aero->tag_clm0  = "clm0slv3_vs_mach_alpha";
        aero->tag_clmq  = "clmqslv3_vs_mach";

        tvc->mtvc    = 1;
        tvc->gtvc    = 1.0;
        tvc->parm    = 5.0;
        tvc->del_max = 15.0;

        ctrl->maut    = 53;
        ctrl->delimx  = 15.0;
        ctrl->drlimx  = 15.0;
        ctrl->gnmax   = 5.0;
        ctrl->gymax   = 5.0;
        ctrl->ancomx  = 0.1;
        ctrl->alcomx  = 0.0;

        ins->mins       = 0;          // ideal: outputs == truth
        ins->bias_accel = osk::Vec(0, 0, 0);
        ins->bias_gyro  = osk::Vec(0, 0, 0);

        // Horizontal launch to avoid gimbal lock
        newt->lonx0  = 0.0; newt->latx0 = 0.0; newt->alt0 = 1000.0;
        newt->dvbe0  = 200.0; newt->psivdx0 = 0.0; newt->thtvdx0 = 0.0;
        kin->psibdx0 = 0.0; kin->thtbdx0 = 0.0; kin->phibdx0 = 0.0;
        eul->ppx0 = eul->qqx0 = eul->rrx0 = 0.0;
        forc->FAPB_ext = osk::Vec(0, 0, 0);
        forc->FMB_ext  = osk::Vec(0, 0, 0);

        env->getsFrom(newt);
        aero->getsFrom(env, kin, prop, tvc);
        tvc->getsFrom(prop);
        forc->getsFrom(prop, aero, env, tvc);
        newt->getsFrom(env, kin, forc, prop);
        kin->getsFrom(env, newt, eul);
        eul->getsFrom(forc, prop, kin);
        ins->getsFrom(newt, eul, kin);
        if (use_ins)
            ctrl->getsFrom(env, newt, eul, aero, prop, ins);
        else
            ctrl->getsFrom(env, newt, eul, aero, prop);

        // Order matters: INS must update AFTER newt/eul/kin to read fresh
        // truth, BEFORE ctrl so ctrl reads fresh INS outputs.
        std::vector<osk::Block*> stage0 = { kin, env, prop, aero, tvc, forc,
                                             newt, eul, ins, ctrl };
        std::vector<std::vector<osk::Block*>> stages = { stage0 };

        double dts[] = { 0.01 };
        osk::Sim sim(dts, 2.0, stages);
        sim.run();

        double result = ctrl->delecx;

        delete ins; delete newt; delete kin; delete eul; delete forc;
        delete ctrl; delete tvc; delete aero; delete prop; delete env;
        return result;
    };

    double delecx_direct  = build(false);  // Control reads Newton/Euler directly
    double delecx_via_ins = build(true);   // Control reads INS (ideal -> same)

    std::printf("  Control direct (no INS):  delecx = %+.6f deg\n", delecx_direct);
    std::printf("  Control via ideal INS:    delecx = %+.6f deg\n", delecx_via_ins);
    std::printf("  difference: %.3e deg\n", std::fabs(delecx_direct - delecx_via_ins));

    // With ideal INS, the two paths should give identical results (modulo
    // numerical precision and potentially a one-step lag if INS updates
    // after Control would read sensors).
    bool ok = std::fabs(delecx_direct - delecx_via_ins) < 0.001;

    std::printf("  %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

} // anon

int main() {
    int fails = 0;
    fails += run_ideal_ins();
    fails += run_biased_ins();
    fails += run_gyro_drift();
    fails += run_moving_vehicle_tracking();
    fails += run_control_via_ins();
    std::printf("\n=== Total failures: %d ===\n", fails);
    return fails == 0 ? 0 : 1;
}
