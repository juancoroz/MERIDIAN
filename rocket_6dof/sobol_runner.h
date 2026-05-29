//  sobol_runner.h  --  Sobol sensitivity-index sweep
//
//  Computes first-order (S_i) and total-effect (ST_i) Sobol sensitivity
//  indices for each (input, output) pair using the pick-freeze
//  estimator due to Saltelli (2010).
//
//  Reuses the existing monte_carlo schema, with an optional override
//  in a sobol sub-section:
//
//    "sobol": {
//      "n_base":  1024,    // base sample count N
//      "seed":    42       // optional, defaults to monte_carlo.seed
//    }
//
//  Total simulations run = N * (k + 2) where k = number of MC variables.
//  Per-output Sobol indices are written to a CSV file.  At N=1024 with
//  k=9, that's 11264 sims = ~5.5 min at 30 ms/run.
//
//  Outputs:
//    sobol_indices.csv    rows = (output, input) pairs, cols = S1, ST
//    sobol_summary.txt    formatted human-readable table
//
//  Method reference: Saltelli et al., "Variance based sensitivity
//  analysis of model output. Design and estimator for the total
//  sensitivity index", Comp. Phys. Comm. 181:259-270 (2010).

#ifndef ROCKET6DOF_SOBOL_RUNNER_H
#define ROCKET6DOF_SOBOL_RUNNER_H

#include <string>

namespace rocket6dof {

int run_sobol(const std::string& config_path,
              const std::string& indices_csv = "sobol_indices.csv",
              const std::string& summary_txt = "sobol_summary.txt");

// Pilot-N advisory.  Runs Sobol at a small N (default 64) with
// bootstrap enabled, measures the worst-case ST CI half-width
// across all (output, input) pairs, and prints a recommended
// production N to hit a target half-width (sobol.target_ci_width
// in the config, default 0.05).  Does not write CSV/summary files.
int suggest_n(const std::string& config_path);

} // namespace rocket6dof

#endif
