//  mc_runner.h  --  Monte Carlo simulation driver
//
//  Reads a mission config containing a "monte_carlo" section, runs N
//  simulations with each input variable drawn from its specified
//  probability distribution, and writes per-run results to CSV files.
//
//  monte_carlo schema:
//
//    "monte_carlo": {
//      "seed":   12345,                // top-level seed
//      "n_runs": 100,                  // number of samples
//      "variables": [
//        {
//          "path":     "propulsion.fmass0",
//          "distribution": "normal",
//          "mean": 500.0,
//          "stddev": 10.0
//        },
//        {
//          "path":     "launch.thtvdx0",
//          "distribution": "uniform",
//          "low":  88.5,
//          "high": 89.5
//        }
//      ],
//      "outputs": ["alt", "dvbi", "fpa", "theta", "psi"]
//    }
//
//  Supported "path" syntax: any path that json::Value::set_path can
//  navigate, e.g. "propulsion.fmass0" or "ins.bias_accel[0]".
//
//  Supported output names (case-sensitive):
//    Trajectory: alt, dvbi, dvbe, fpa, heading, lon, lat
//    Vehicle:    vmass, fmassr, fmasse
//    Attitude:   phi, theta, psi, alpha, beta
//    Body rates: p_deg, q_deg, r_deg
//    Control:    delecx, delrcx
//    INS:        ins_pos_err, ins_vel_err, ins_att_err
//    Sensors:    gps_updates, startrack_updates
//    Time:       t_end                     (sim termination time)
//
//  Reproducibility:
//    Per-variable RNG uses seed (top_seed XOR hash(path)).  The
//    per-run sample is drawn from one RNG per variable.  This means
//    adding/removing a variable does not perturb the samples of
//    unrelated variables -- a property useful when iterating on the
//    MC setup.
//
//  Output files (paths relative to working directory):
//    monte_carlo_inputs.csv   one row per run, columns = variable paths
//    monte_carlo_outputs.csv  one row per run, columns = output names
//                             plus a 'run_index' column at the start
//    monte_carlo_summary.txt  human-readable summary statistics

#ifndef ROCKET6DOF_MC_RUNNER_H
#define ROCKET6DOF_MC_RUNNER_H

#include "json.h"

#include <string>

namespace rocket6dof {

// Run the Monte Carlo sweep described in the config file.  Returns 0
// on success, nonzero on any error.  Prints a progress line per run
// and a summary at the end.
int run_monte_carlo(const std::string& config_path,
                    const std::string& outputs_csv = "monte_carlo_outputs.csv",
                    const std::string& inputs_csv  = "monte_carlo_inputs.csv",
                    const std::string& summary_txt = "monte_carlo_summary.txt");

} // namespace rocket6dof

#endif
