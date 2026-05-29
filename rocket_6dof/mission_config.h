//  mission_config.h  --  Apply JSON config to rocket_6dof blocks
//
//  Loads a JSON config file and applies it to the simulation blocks.
//  Missing keys retain whatever default the block's constructor set,
//  so a partial config file works -- you only override what you need.
//
//  Top-level config schema (all sections optional):
//
//    {
//      "propulsion":   { ... },     // see Propulsion class for keys
//      "aerodynamics": { ... },
//      "tvc":          { ... },
//      "rcs":          { ... },
//      "control":      { ... },
//      "guidance":     { ... },
//      "ins":          { ... },
//      "gps":          { ... },
//      "startrack":    { ... },
//      "intercept":    { ... },
//      "launch":       { lonx0, latx0, alt0, dvbe0, psivdx0, thtvdx0,
//                        psibdx0, thtbdx0, phibdx0 },
//      "logger":       { sample_period, filename },
//      "sim":          { dt, t_max }
//    }
//
//  Each block's section accepts the public scalar members of that
//  block.  Vector members (e.g. ins.bias_accel) accept a 3-element
//  array: [x, y, z].
//
//  Usage:
//
//      auto cfg = rocket6dof::json::parse_file("mission.json");
//      apply_config(cfg, prop, aero, tvc, rcs, ctrl, guid,
//                   ins, gps, star, icpt, newt, kin, eul, forc);
//
//  The function applies the config in place; nothing is returned.
//  Unknown keys are silently ignored (so you can have extra fields
//  for human comments, version stamps, etc).

#ifndef ROCKET6DOF_MISSION_CONFIG_H
#define ROCKET6DOF_MISSION_CONFIG_H

#include "json.h"

namespace rocket6dof {

class Propulsion;
class Aerodynamics;
class TVC;
class RCS;
class Control;
class Guidance;
class INS;
class GPS;
class Startrack;
class Intercept;
class Newton;
class Kinematics;
class Euler;
class Forces;

// Apply a parsed JSON config to all the simulation blocks.
// Returns 0 on success, nonzero if a fatal config error was detected.
// (Currently always returns 0; this is a placeholder for future schema
// validation.)
int apply_config(const json::Value& cfg,
                 Propulsion*   prop,
                 Aerodynamics* aero,
                 TVC*          tvc,
                 RCS*          rcs,
                 Control*      ctrl,
                 Guidance*     guid,
                 INS*          ins,
                 GPS*          gps,
                 Startrack*    star,
                 Intercept*    icpt,
                 Newton*       newt,
                 Kinematics*   kin,
                 Euler*        eul,
                 Forces*       forc);

// Returns the configured simulation time-step (default 0.01 s) and
// maximum simulation duration (default = the Intercept block's t_max).
// These are read from the "sim" section of the config.
double config_sim_dt(const json::Value& cfg, double def = 0.01);
double config_sim_tmax(const json::Value& cfg, double def = 26.0);

} // namespace rocket6dof

#endif
