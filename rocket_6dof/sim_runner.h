//  sim_runner.h  --  Single-mission build/run/read helper
//
//  Factors the per-run lifecycle out of mc_runner so it can be shared
//  with sobol_runner.  Each call:
//
//    1. Allocates fresh block instances (no cross-run state).
//    2. Applies the JSON config tree via apply_config().
//    3. Wires blocks via getsFrom().
//    4. Runs the sim with stdout/stderr silenced.
//    5. Reads the requested output names from the block instances.
//    6. Frees all blocks.
//
//  Supported output names (case-sensitive): see read_output() in the
//  .cpp file.  Returns SingleRunResult with values aligned to the
//  output_names vector; NaN for unrecognized names.

#ifndef ROCKET6DOF_SIM_RUNNER_H
#define ROCKET6DOF_SIM_RUNNER_H

#include "json.h"

#include <string>
#include <vector>

namespace rocket6dof {

struct SingleRunResult {
    std::vector<double> outputs;
    double t_end;
    bool   nan_seen;
};

// Run a single mission with the provided config tree, returning the
// requested outputs.  stdout/stderr from the kernel are silenced
// (redirected to /dev/null) during the sim's run() call.
SingleRunResult run_single_mission(const json::Value& cfg,
                                   const std::vector<std::string>& output_names);

} // namespace rocket6dof

#endif
