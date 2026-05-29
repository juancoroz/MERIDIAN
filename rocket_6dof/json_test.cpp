//  json_test.cpp  --  Verify the minimal JSON parser

#include "json.h"
#include "mission_config.h"
#include "propulsion.h"
#include "aerodynamics.h"
#include "tvc.h"
#include "rcs.h"
#include "control.h"
#include "guidance.h"
#include "ins.h"
#include "gps.h"
#include "startrack.h"
#include "intercept.h"
#include "newton.h"
#include "kinematics.h"
#include "euler.h"
#include "forces.h"

#include <cmath>
#include <cstdio>
#include <string>

using namespace rocket6dof::json;
using namespace rocket6dof;

namespace {

int fails = 0;

bool check(bool cond, const char* desc) {
    std::printf("  %s %s\n", cond ? "OK  " : "FAIL", desc);
    if (!cond) fails++;
    return cond;
}

int test_primitives() {
    std::printf("\n=== Test 1: Primitive types ===\n");

    Value v = parse_string("42");
    check(v.isNumber(),                  "parse '42' -> NUMBER");
    check(v.asNumber() == 42.0,          "  value == 42.0");

    v = parse_string("-3.14e2");
    check(std::fabs(v.asNumber() + 314.0) < 1e-9, "parse '-3.14e2' -> -314.0");

    v = parse_string("\"hello world\"");
    check(v.isString(),                  "parse string -> STRING");
    check(v.asString() == "hello world", "  value == 'hello world'");

    v = parse_string("true");
    check(v.isBool() && v.asBool(),      "parse 'true' -> BOOL true");

    v = parse_string("false");
    check(v.isBool() && !v.asBool(),     "parse 'false' -> BOOL false");

    v = parse_string("null");
    check(v.isNull(),                    "parse 'null' -> NUL");

    return 0;
}

int test_objects_arrays() {
    std::printf("\n=== Test 2: Objects and arrays ===\n");

    Value v = parse_string("{ \"a\": 1, \"b\": 2.5, \"c\": \"x\" }");
    check(v.isObject(),                   "parse object");
    check(v.size() == 3,                  "  size == 3");
    check(v["a"].asNumber() == 1.0,       "  v['a'] == 1");
    check(v["b"].asNumber() == 2.5,       "  v['b'] == 2.5");
    check(v["c"].asString() == "x",       "  v['c'] == 'x'");
    check(!v["missing"].exists(),         "  v['missing'].exists() == false");
    check(v["missing"].asNumber(7.0) == 7.0, "  v['missing'].asNumber(7) returns default");

    v = parse_string("[10, 20, 30]");
    check(v.isArray(),                    "parse array");
    check(v.size() == 3,                  "  size == 3");
    check(v[0].asNumber() == 10.0,        "  v[0] == 10");
    check(v[2].asNumber() == 30.0,        "  v[2] == 30");

    // Nested
    v = parse_string("{ \"outer\": { \"inner\": [1, 2, [3, 4]] } }");
    check(v["outer"]["inner"][2][1].asNumber() == 4.0,
          "  nested object/array access");

    return 0;
}

int test_comments_and_trailing() {
    std::printf("\n=== Test 3: Comments and trailing commas ===\n");

    Value v = parse_string(
        "// This is a comment\n"
        "{\n"
        "  \"x\": 1,  // trailing comment\n"
        "  \"y\": 2,\n"   // trailing comma below
        "}\n"
    );
    check(v.isObject(),                 "comments and trailing commas accepted");
    check(v["x"].asNumber() == 1.0,     "  v['x'] == 1");
    check(v["y"].asNumber() == 2.0,     "  v['y'] == 2");

    return 0;
}

int test_escapes() {
    std::printf("\n=== Test 4: String escapes ===\n");

    Value v = parse_string("\"line1\\nline2\\ttab\\\\\\\"\"");
    std::string s = v.asString();
    check(s == "line1\nline2\ttab\\\"",  "  escapes parsed correctly");

    return 0;
}

int test_error_lines() {
    std::printf("\n=== Test 5: Error line/column reporting ===\n");

    bool caught = false;
    try {
        parse_string("{\n  \"x\": 1,\n  \"y\": @broken\n}");
    } catch (const ParseError& e) {
        caught = true;
        std::printf("  caught ParseError at line %d col %d: %s\n",
                    e.line, e.col, e.what());
        check(e.line == 3, "  error reported on correct line");
    }
    check(caught, "  malformed input throws ParseError");

    return 0;
}

int test_realistic_config() {
    std::printf("\n=== Test 6: Realistic mission config snippet ===\n");

    Value v = parse_string(R"(
{
  "propulsion": {
    "mprop": 3,
    "vmass0": 1000.0,
    "fmass0": 500.0,
    "spi": 300.0,
    "fuel_flow_rate": 20.0
  },
  "control": {
    "maut": 60,
    "tau_att": 1.0,
    "q_max": 15.0,
    "theta_com_inertial": 89.0,
    "psi_com_inertial": 90.0
  },
  "guidance": {
    "mguide": 2,
    "t_att_start": 3.0,
    "t_att_end": 15.0,
    "theta_com_start": 89.0,
    "theta_com_end": 70.0
  },
  "launch": {
    "lonx0": -80.0,
    "latx0": 28.5,
    "alt0": 0.0,
    "dvbe0": 0.01,
    "psivdx0": 90.0,
    "thtvdx0": 89.0
  }
}
)");

    check(v["propulsion"]["vmass0"].asNumber() == 1000.0,    "propulsion.vmass0 = 1000");
    check(v["propulsion"]["mprop"].asInt() == 3,             "propulsion.mprop = 3");
    check(v["control"]["theta_com_inertial"].asNumber() == 89.0,
          "control.theta_com_inertial = 89");
    check(v["guidance"]["mguide"].asInt() == 2,              "guidance.mguide = 2");
    check(v["launch"]["lonx0"].asNumber() == -80.0,          "launch.lonx0 = -80");

    // Missing keys return defaults
    check(v["rcs"]["mrcs_moment"].asInt(99) == 99,           "missing keys take defaults");
    check(v["control"]["nonexistent"].asNumber(42.0) == 42.0,
          "missing sub-keys take defaults");

    return 0;
}

// ---- Test 7: apply_config end-to-end against real blocks ----
// Verify that parsing a config and applying it to actual block
// instances produces the expected member values, including:
//   - scalar overrides (numbers, bools)
//   - vector overrides ([x, y, z] arrays into osk::Vec)
//   - string overrides (aero file tags)
//   - missing keys leaving block defaults intact
int test_apply_config_integration() {
    std::printf("\n=== Test 7: apply_config integration ===\n");

    Propulsion*   prop = new Propulsion();
    Aerodynamics* aero = new Aerodynamics();
    TVC*          tvc  = new TVC();
    RCS*          rcs  = new RCS();
    Control*      ctrl = new Control();
    Guidance*     guid = new Guidance();
    INS*          ins  = new INS();
    GPS*          gps  = new GPS();
    Startrack*    star = new Startrack();
    Intercept*    icpt = new Intercept();
    Newton*       newt = new Newton();
    Kinematics*   kin  = new Kinematics();
    Euler*        eul  = new Euler();
    Forces*       forc = new Forces();

    // Capture a few defaults before applying config
    double default_tvc_parm   = tvc->parm;       // un-set in config below
    double default_ctrl_gnmax = ctrl->gnmax;     // un-set in config below

    Value cfg = parse_string(R"(
{
  "propulsion": {
    "vmass0": 2500.0,
    "fmass0": 1200.0,
    "fuel_flow_rate": 40.0
  },
  "aerodynamics": {
    "maero": 13,
    "refa": 1.5,
    "aero_file": "custom_aero.txt",
    "tag_ca0": "ca0_custom_vs_mach"
  },
  "tvc": { "mtvc": 1, "del_max": 12.0 },
  "control": {
    "maut": 60,
    "tau_att": 0.5,
    "theta_com_inertial": 75.0
  },
  "guidance": {
    "mguide": 2,
    "t_att_start": 5.0,
    "t_att_end": 20.0,
    "theta_com_start": 89.0,
    "theta_com_end": 60.0
  },
  "ins": {
    "mins": 1,
    "bias_accel": [0.002, -0.001, 0.0005],
    "bias_gyro":  [1.0e-5, 2.0e-5, 3.0e-5]
  },
  "gps": { "mgps": 2, "gps_step": 0.5, "rpos": 8.0, "noise_seed": 99 },
  "intercept": {
    "check_alt_max": true,
    "check_time_max": false,
    "alt_max": 50000.0
  },
  "launch": {
    "lonx0": -90.0, "latx0": 30.0, "alt0": 100.0,
    "psivdx0": 0.0, "thtvdx0": 88.0,
    "psibdx0": 0.0, "thtbdx0": 88.0
  }
}
)");

    int rc = apply_config(cfg, prop, aero, tvc, rcs, ctrl, guid,
                          ins, gps, star, icpt, newt, kin, eul, forc);
    check(rc == 0, "apply_config returns 0");

    // Scalar overrides
    check(prop->vmass0          == 2500.0, "prop.vmass0 = 2500");
    check(prop->fmass0          == 1200.0, "prop.fmass0 = 1200");
    check(prop->fuel_flow_rate  == 40.0,   "prop.fuel_flow_rate = 40");
    check(aero->maero           == 13,     "aero.maero = 13");
    check(aero->refa            == 1.5,    "aero.refa = 1.5");
    check(tvc->mtvc             == 1,      "tvc.mtvc = 1");
    check(tvc->del_max          == 12.0,   "tvc.del_max = 12");
    check(ctrl->maut            == 60,     "ctrl.maut = 60");
    check(ctrl->tau_att         == 0.5,    "ctrl.tau_att = 0.5");
    check(ctrl->theta_com_inertial == 75.0, "ctrl.theta_com_inertial = 75");
    check(guid->mguide          == 2,      "guid.mguide = 2");
    check(guid->theta_com_end   == 60.0,   "guid.theta_com_end = 60");
    check(gps->mgps             == 2,      "gps.mgps = 2");
    check(gps->gps_step         == 0.5,    "gps.gps_step = 0.5");
    check(gps->rpos             == 8.0,    "gps.rpos = 8");
    check(gps->noise_seed       == 99,     "gps.noise_seed = 99");

    // Bool overrides
    check(icpt->check_alt_max   == true,   "icpt.check_alt_max = true");
    check(icpt->check_time_max  == false,  "icpt.check_time_max = false");
    check(icpt->alt_max         == 50000.0,"icpt.alt_max = 50000");

    // String overrides
    check(aero->aero_file == "custom_aero.txt",     "aero.aero_file overridden");
    check(aero->tag_ca0   == "ca0_custom_vs_mach",  "aero.tag_ca0 overridden");

    // Vector overrides (osk::Vec from [x, y, z])
    check(std::fabs(ins->bias_accel.x - 0.002)   < 1e-9, "ins.bias_accel.x = 0.002");
    check(std::fabs(ins->bias_accel.y - (-0.001))< 1e-9, "ins.bias_accel.y = -0.001");
    check(std::fabs(ins->bias_accel.z - 0.0005)  < 1e-9, "ins.bias_accel.z = 0.0005");
    check(std::fabs(ins->bias_gyro.x - 1.0e-5)   < 1e-12,"ins.bias_gyro.x = 1e-5");
    check(std::fabs(ins->bias_gyro.z - 3.0e-5)   < 1e-12,"ins.bias_gyro.z = 3e-5");

    // Launch ICs into Newton + Kinematics
    check(newt->lonx0  == -90.0, "newt.lonx0 = -90");
    check(newt->latx0  ==  30.0, "newt.latx0 = 30");
    check(newt->alt0   == 100.0, "newt.alt0 = 100");
    check(newt->thtvdx0 == 88.0, "newt.thtvdx0 = 88");
    check(kin->thtbdx0  == 88.0, "kin.thtbdx0 = 88");

    // Defaults preserved when config doesn't mention a key
    check(tvc->parm   == default_tvc_parm,
          "tvc.parm preserved (not in config)");
    check(ctrl->gnmax == default_ctrl_gnmax,
          "ctrl.gnmax preserved (not in config)");

    // Sim config getters
    check(config_sim_dt(cfg, 0.01) == 0.01,    "sim.dt missing -> default 0.01");
    check(config_sim_tmax(cfg, 25.0) == 25.0,  "sim.t_max missing -> default 25");

    Value cfg2 = parse_string(R"({ "sim": { "dt": 0.005, "t_max": 60.0 } })");
    check(config_sim_dt(cfg2, 0.01)   == 0.005, "sim.dt loaded from config");
    check(config_sim_tmax(cfg2, 25.0) == 60.0,  "sim.t_max loaded from config");

    delete forc; delete eul; delete kin; delete newt;
    delete icpt; delete star; delete gps; delete ins;
    delete guid; delete ctrl; delete rcs; delete tvc;
    delete aero; delete prop;
    return 0;
}

} // anon

int main() {
    test_primitives();
    test_objects_arrays();
    test_comments_and_trailing();
    test_escapes();
    test_error_lines();
    test_realistic_config();
    test_apply_config_integration();
    std::printf("\n=== Total failures: %d ===\n", fails);
    return fails == 0 ? 0 : 1;
}
