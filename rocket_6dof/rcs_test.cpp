//  rcs_test.cpp -- Verify RCS thruster outputs
#include "../osk/osk.h"
#include "environment.h"
#include "newton.h"
#include "kinematics.h"
#include "euler.h"
#include "forces.h"
#include "propulsion.h"
#include "ins.h"
#include "guidance.h"
#include "rcs.h"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace rocket6dof;

namespace {

// TEST1_MARKER

int run_rcs_off() {
    std::printf("\n=== Test 1: RCS off ===\n");
    RCS* rcs = new RCS();
    rcs->mrcs_moment = 0;
    rcs->mrcs_force  = 0;
    rcs->update();
    bool ok =    rcs->FMRCS.mag() == 0.0
              && rcs->FARCS.mag() == 0.0;
    std::printf("  FMRCS=(%g,%g,%g)  FARCS=(%g,%g,%g)\n",
                rcs->FMRCS.x, rcs->FMRCS.y, rcs->FMRCS.z,
                rcs->FARCS.x, rcs->FARCS.y, rcs->FARCS.z);
    std::printf("  %s\n", ok ? "PASS" : "FAIL");
    delete rcs;
    return ok ? 0 : 1;
}

// TEST2_MARKER

int run_rcs_prop_saturation() {
    std::printf("\n=== Test 2: Proportional saturation ===\n");
    Propulsion* prop = new Propulsion();
    Kinematics* kin  = new Kinematics();
    Newton*     newt = new Newton();
    Euler*      eul  = new Euler();
    INS*        ins  = new INS();
    RCS*        rcs  = new RCS();

    prop->IBBB = osk::Mat(5,0,0, 0,50,0, 0,0,50);
    kin->phibdx = 100.0;  // 100 deg roll error
    kin->thtbdx = 0.0; kin->psibdx = 0.0;

    rcs->mrcs_moment   = 11;    // proportional, Euler mode
    rcs->rcs_freq      = 1.0;
    rcs->rcs_zeta      = 0.707;
    rcs->roll_mom_max  = 10.0;
    rcs->pitch_mom_max = 10.0;
    rcs->yaw_mom_max   = 10.0;
    rcs->phibdcomx     = 0.0;
    rcs->getsFrom(prop, newt, kin, ins, nullptr);

    rcs->update();
    std::printf("  e_roll=%.3f  FMRCS.x=%.3f Nm  (expect -10.0)\n",
                rcs->e_roll, rcs->FMRCS.x);
    bool ok = std::fabs(rcs->FMRCS.x + 10.0) < 1e-9;
    std::printf("  %s\n", ok ? "PASS" : "FAIL");

    delete rcs; delete ins; delete eul; delete newt; delete kin; delete prop;
    return ok ? 0 : 1;
}

// TEST3_MARKER

int run_rcs_schmitt() {
    std::printf("\n=== Test 3: Schmitt trigger on/off ===\n");
    Propulsion* prop = new Propulsion();
    Kinematics* kin  = new Kinematics();
    Newton*     newt = new Newton();
    Euler*      eul  = new Euler();
    INS*        ins  = new INS();
    RCS*        rcs  = new RCS();

    prop->IBBB = osk::Mat(1,0,0, 0,1,0, 0,0,1);
    kin->phibdx = 0.0; kin->thtbdx = 0.0; kin->psibdx = 0.0;

    rcs->mrcs_moment   = 21;    // on-off, Euler mode
    rcs->dead_zone     = 1.0;
    rcs->hysteresis    = 0.0;
    rcs->rcs_tau       = 0.0;
    rcs->roll_mom_max  = 5.0;
    rcs->pitch_mom_max = 5.0;
    rcs->yaw_mom_max   = 5.0;
    rcs->getsFrom(prop, newt, kin, ins, nullptr);

    struct Case { double cmd; double expected_FMRCS_x; const char* note; };
    Case cases[] = {
        {  0.1,   0.0, "well inside dead zone" },
        { 10.0,   5.0, "large +" },
        { -10.0, -5.0, "large -" },
    };

    bool all_ok = true;
    for (const auto& c : cases) {
        rcs->phibdcomx = c.cmd;
        rcs->update();
        // Schmitt depends on previous state; one call may not settle.
        // Trigger again so trend stabilizes:
        rcs->update();
        bool ok = std::fabs(rcs->FMRCS.x - c.expected_FMRCS_x) < 1e-9;
        std::printf("  cmd=%6.2f -> FMRCS.x=%+.2f Nm  (expect %+.2f, %s)  %s\n",
                    c.cmd, rcs->FMRCS.x, c.expected_FMRCS_x, c.note,
                    ok ? "OK" : "FAIL");
        if (!ok) all_ok = false;
    }
    std::printf("  %s\n", all_ok ? "PASS" : "FAIL");
    delete rcs; delete ins; delete eul; delete newt; delete kin; delete prop;
    return all_ok ? 0 : 1;
}

// TEST4_MARKER

int run_rcs_prop_closed_loop() {
    std::printf("\n=== Test 4: Proportional closed-loop roll convergence ===\n");

    Environment* env  = new Environment();
    Propulsion*  prop = new Propulsion();
    Forces*      forc = new Forces();
    Euler*       eul  = new Euler();
    Kinematics*  kin  = new Kinematics();
    Newton*      newt = new Newton();
    INS*         ins  = new INS();
    RCS*         rcs  = new RCS();

    prop->mprop  = 0;  prop->vmass0 = 100.0;
    prop->moi_roll_0  = 5.0;  prop->moi_roll_1  = 5.0;
    prop->moi_trans_0 = 50.0; prop->moi_trans_1 = 50.0;

    eul->ppx0 = eul->qqx0 = eul->rrx0 = 0.0;
    newt->lonx0 = 0; newt->latx0 = 0; newt->alt0 = 500000.0;
    newt->dvbe0 = 0.0;
    kin->psibdx0 = 0; kin->thtbdx0 = 0; kin->phibdx0 = 30.0;

    rcs->mrcs_moment   = 11;
    rcs->rcs_freq      = 1.0;
    rcs->rcs_zeta      = 1.0;
    rcs->roll_mom_max  = 100.0;
    rcs->pitch_mom_max = 100.0;
    rcs->yaw_mom_max   = 100.0;
    rcs->phibdcomx     = 0.0;

    forc->FAPB_ext = osk::Vec(0,0,0); forc->FMB_ext = osk::Vec(0,0,0);
    ins->mins = 0;

    env->getsFrom(newt);
    forc->getsFrom(prop, nullptr, env, nullptr, rcs);
    newt->getsFrom(env, kin, forc, prop);
    kin->getsFrom(env, newt, eul);
    eul->getsFrom(forc, prop, kin);
    ins->getsFrom(newt, eul, kin);
    rcs->getsFrom(prop, newt, kin, ins, nullptr);

    std::vector<osk::Block*> stage0 = { kin, env, prop, ins, rcs, forc,
                                         newt, eul };
    std::vector<std::vector<osk::Block*>> stages = { stage0 };
    double dts[] = { 0.01 };
    osk::Sim sim(dts, 15.0, stages);
    sim.run();

    double final_roll = kin->phibdx;
    std::printf("  Initial roll = 30 deg\n");
    std::printf("  Final roll   = %.4f deg  (expect close to 0)\n", final_roll);
    std::printf("  Final FMRCS.x = %+.3f Nm\n", rcs->FMRCS.x);

    bool ok = std::fabs(final_roll) < 0.5;   // 0.5 deg residual
    std::printf("  %s\n", ok ? "PASS" : "FAIL");

    delete rcs; delete ins; delete newt; delete kin; delete eul;
    delete forc; delete prop; delete env;
    return ok ? 0 : 1;
}

// TEST5_MARKER

int run_rcs_force_aggregation() {
    std::printf("\n=== Test 5: RCS torque appears in Forces.FMB ===\n");

    Environment* env  = new Environment();
    Propulsion*  prop = new Propulsion();
    Forces*      forc = new Forces();
    Kinematics*  kin  = new Kinematics();
    Newton*      newt = new Newton();
    INS*         ins  = new INS();
    RCS*         rcs  = new RCS();

    prop->IBBB = osk::Mat(5,0,0, 0,50,0, 0,0,50);
    kin->phibdx = 30.0;       // roll error
    kin->thtbdx = 0.0; kin->psibdx = 0.0;

    rcs->mrcs_moment   = 11;
    rcs->rcs_freq      = 1.0;
    rcs->rcs_zeta      = 0.707;
    rcs->roll_mom_max  = 100.0;
    rcs->pitch_mom_max = 100.0;
    rcs->yaw_mom_max   = 100.0;
    rcs->phibdcomx     = 0.0;
    rcs->getsFrom(prop, newt, kin, ins, nullptr);

    forc->FAPB_ext = osk::Vec(0,0,0);
    forc->FMB_ext  = osk::Vec(0,0,0);
    forc->getsFrom(prop, nullptr, env, nullptr, rcs);

    rcs->update();
    forc->update();

    std::printf("  rcs->FMRCS = (%+.3f, %+.3f, %+.3f)\n",
                rcs->FMRCS.x, rcs->FMRCS.y, rcs->FMRCS.z);
    std::printf("  forc->FMB  = (%+.3f, %+.3f, %+.3f)\n",
                forc->FMB.x,  forc->FMB.y,  forc->FMB.z);

    bool ok =    std::fabs(forc->FMB.x - rcs->FMRCS.x) < 1e-9
              && std::fabs(forc->FMB.y - rcs->FMRCS.y) < 1e-9
              && std::fabs(forc->FMB.z - rcs->FMRCS.z) < 1e-9;
    std::printf("  %s\n", ok ? "PASS" : "FAIL");

    delete rcs; delete ins; delete newt; delete kin;
    delete forc; delete prop; delete env;
    return ok ? 0 : 1;
}

} // anon

int main() {
    int fails = 0;
    fails += run_rcs_off();
    fails += run_rcs_prop_saturation();
    fails += run_rcs_schmitt();
    fails += run_rcs_prop_closed_loop();
    fails += run_rcs_force_aggregation();
    std::printf("\n=== Total failures: %d ===\n", fails);
    return fails == 0 ? 0 : 1;
}
