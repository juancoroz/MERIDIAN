//  parallel_runner.h  --  Subprocess-based parallel execution
//
//  Given a sample matrix (rows = input vectors), a config path, a list
//  of output names, and a worker count, fan out execution across
//  mission_worker subprocesses and aggregate the results back into a
//  per-row Y matrix.
//
//  Each worker is a separate OS process with its own copy of the OSK
//  kernel's static globals (osk::State::t etc.), so simulations are
//  fully independent.  Communication is via temporary CSV files in
//  the working directory (or /tmp, configurable).
//
//  Use when n_workers >= 2.  For n_workers == 1, callers should use
//  run_single_mission() directly to avoid the IPC overhead.
//
//  ParallelResult.row_outputs[r][o] = j-th output of run row r
//  ParallelResult.nan_flags[r]      = true if any output was NaN

#ifndef ROCKET6DOF_PARALLEL_RUNNER_H
#define ROCKET6DOF_PARALLEL_RUNNER_H

#include <string>
#include <vector>

namespace rocket6dof {

struct ParallelResult {
    // row_outputs[row][output] -- N rows, n_out outputs each
    std::vector<std::vector<double>> row_outputs;
    std::vector<bool>                nan_flags;
    bool   ok;        // false if any worker failed
    double elapsed_s; // wallclock for the parallel phase
};

// Run a chunked-parallel sweep.
//
// inputs:
//   config_path     -- the JSON config file each worker reads
//   input_paths     -- variable paths corresponding to sample columns
//   sample_matrix   -- N rows, each with input_paths.size() values
//   output_names    -- outputs to extract per run
//   n_workers       -- number of subprocess workers (>= 1)
//   worker_bin      -- path to the mission_worker binary
//   tmp_dir         -- where to place temp CSV files (e.g. "/tmp")
//   progress_label  -- short string for progress lines (e.g. "MC", "Sobol-A")
//
// returns ParallelResult.  On failure, ok=false and stderr is printed.
ParallelResult run_chunked_parallel(
    const std::string& config_path,
    const std::vector<std::string>& input_paths,
    const std::vector<std::vector<double>>& sample_matrix,
    const std::vector<std::string>& output_names,
    int n_workers,
    const std::string& worker_bin = "./mission_worker",
    const std::string& tmp_dir    = "/tmp",
    const std::string& progress_label = "parallel");

// Resolve the user-specified n_workers value to a concrete count.
//
//   requested == 0  -> auto-detect: std::thread::hardware_concurrency(),
//                       capped at AUTO_CAP (default 8) to avoid spawning
//                       a huge number of workers on big servers.
//   requested >= 1  -> use as-is (no cap; user knows their environment).
//   requested  < 0  -> clamp to 1.
//
// If hardware_concurrency() reports 0 (rare; happens in some containers),
// the auto path falls back to 1.
//
// Returns the resolved count plus a short string explaining the choice,
// suitable for inclusion in startup logs.
struct ResolvedWorkers {
    int n;
    std::string source;  // e.g. "auto-detected: 8 cores, capped at 8"
};
ResolvedWorkers resolve_n_workers(int requested, int auto_cap = 8);

} // namespace rocket6dof

#endif
