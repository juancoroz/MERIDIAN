//  propulsion.h  --  Multi-stage rocket propulsion + mass properties
//
//  Cross-reference to Zipfel Section 10.3.x (hyper[10-49]).
//
//  Multi-stage support: set num_stages > 1 and fill in the per-stage
//  arrays (vmass0_stage[], fmass0_stage[], spi_stage[], etc.).  The
//  block automatically:
//
//    1. Burns fuel of the current stage at mdot_stage[current_stage]
//       until fmasse >= fmass0_stage[current_stage]  (stage burnout)
//    2. Optionally coasts for coast_to_next[current_stage] seconds (no
//       thrust, mass held constant)
//    3. At end of coast: discards dry_mass_dropped[current_stage] of
//       structural mass, switches all parameters to stage N+1, resets
//       fmasse to 0 and resumes burning
//    4. Continues until all stages are exhausted
//
//  Backward compatibility: when num_stages == 1 (default), only the
//  legacy scalar fields are used (vmass0, fmass0, spi, fuel_flow_rate,
//  moi_*, xcg_*).  Existing test code and mission scenarios with a
//  single-stage rocket continue to work unchanged.
//
//  Mode 'mprop':
//     0 = no thrust (commanded off, or all stages exhausted)
//     3 = constant-thrust rocket; engine fires per current stage params
//
//  Operational state (read-only outputs for diagnostics):
//     current_stage:  index 0..num_stages-1 of the active stage
//     phase:          0 = burning
//                     1 = coasting (between stages)
//                     2 = spent (all stages exhausted)

#ifndef ROCKET6DOF_PROPULSION_H
#define ROCKET6DOF_PROPULSION_H

#include "../osk/osk.h"

namespace rocket6dof {

class Environment;

class Propulsion : public osk::Block {
public:
    // Maximum number of propulsion stages this build supports.  Can be
    // overridden at compile time with -DROCKET6DOF_MAX_STAGES=N (e.g.
    // for a 6-stage research vehicle).  Default 4 covers everything
    // currently flying (Saturn V was 3; Falcon Heavy is effectively 2-3).
    // Raising this only changes the per-stage array sizes here; the
    // guidance LTG path is still 3-stage-only (see mission_config.cpp).
#ifdef ROCKET6DOF_MAX_STAGES
    static constexpr int MAX_STAGES = ROCKET6DOF_MAX_STAGES;
#else
    static constexpr int MAX_STAGES = 4;
#endif
    static_assert(MAX_STAGES >= 1, "MAX_STAGES must be >= 1");

    // ---- Inputs ----
    Environment* env;
    void getsFrom(Environment* e) { env = e; }

    // ---- Mode selector ----
    int mprop;

    // ---- Single-stage parameters (backward-compatible interface) ----
    // Used directly when num_stages == 1.  When num_stages > 1, these
    // are overwritten at init() from the *_stage[] arrays at stage 0.
    double vmass0;
    double fmass0;
    double spi;
    double fuel_flow_rate;
    double xcg_0;
    double xcg_1;
    double moi_roll_0;
    double moi_roll_1;
    double moi_trans_0;
    double moi_trans_1;

    // ---- Multi-stage configuration ----
    int    num_stages;
    double vmass0_stage  [MAX_STAGES];
    double fmass0_stage  [MAX_STAGES];
    double spi_stage     [MAX_STAGES];
    double mdot_stage    [MAX_STAGES];
    double xcg_0_stage   [MAX_STAGES];
    double xcg_1_stage   [MAX_STAGES];
    double moi_roll_0_stage [MAX_STAGES];
    double moi_roll_1_stage [MAX_STAGES];
    double moi_trans_0_stage[MAX_STAGES];
    double moi_trans_1_stage[MAX_STAGES];
    double dry_mass_dropped [MAX_STAGES];
    double coast_to_next    [MAX_STAGES];

    // ---- Integrator state ----
    double fmasse;
    double fmassd;

    // ---- Operational state ----
    int    current_stage;
    int    phase;
    double t_stage_start;
    double t_coast_end;

    // ---- Outputs ----
    double   thrust;
    double   vmass;
    double   fmassr;
    double   xcg;
    osk::Mat IBBB;

    Propulsion();
    void init()   override;
    void update() override;
    void rpt()    override;

    ACCESS_FN(double,   thrust)
    ACCESS_FN(double,   vmass)
    ACCESS_FN(double,   fmassr)
    ACCESS_FN(double,   xcg)
    ACCESS_FN(osk::Mat, IBBB)
};

} // namespace rocket6dof

#endif
