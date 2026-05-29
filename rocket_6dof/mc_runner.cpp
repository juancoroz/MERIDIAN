//  mc_runner.cpp  --  Monte Carlo driver implementation
//
//  Per-run flow:
//    1. Copy the baseline JSON config tree.
//    2. For each MC variable, sample its distribution and set the
//       corresponding path in the per-run tree.
//    3. Construct fresh block instances.
//    4. Apply the per-run config (apply_config from mission_config).
//    5. Wire blocks, build stages, run the sim with logging suppressed.
//    6. Read the requested outputs from the block instances.
//    7. Free everything.
//
//  Output rows are buffered in memory and flushed to disk at the end
//  (or on Ctrl-C).  Per-run progress is printed to stdout.

#include "mc_runner.h"
#include "json.h"
#include "distributions.h"
#include "mission_config.h"
#include "propulsion.h"
#include "sim_runner.h"
#include "parallel_runner.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

namespace rocket6dof {

namespace {

//  Per-variable RNG seeding via path-hash
unsigned long path_hash(const std::string& s) {
    // FNV-1a 64-bit, truncated to 32 bits
    unsigned long h = 14695981039346656037UL;
    for (char c : s) {
        h ^= static_cast<unsigned char>(c);
        h *= 1099511628211UL;
    }
    return static_cast<unsigned long>(h & 0xFFFFFFFFu);
}

//  Summary statistics
struct ColStats {
    double mean, stddev, min, max;
    int    nan_count;
};

ColStats compute_stats(const std::vector<double>& col) {
    ColStats s = { 0, 0, 0, 0, 0 };
    int n = 0;
    double sum = 0, sum2 = 0;
    bool first = true;
    for (double v : col) {
        if (!std::isfinite(v)) { s.nan_count++; continue; }
        if (first) { s.min = s.max = v; first = false; }
        else { if (v < s.min) s.min = v; if (v > s.max) s.max = v; }
        sum  += v;
        sum2 += v * v;
        ++n;
    }
    if (n > 0) {
        s.mean   = sum / n;
        double var = sum2 / n - s.mean * s.mean;
        s.stddev = (var > 0) ? std::sqrt(var) : 0.0;
    } else {
        // All values were NaN -- report NaN rather than 0 so the
        // pathology is visible (zero would silently look like a
        // valid result).
        double nan_v = std::nan("");
        s.mean   = nan_v;
        s.stddev = nan_v;
        s.min    = nan_v;
        s.max    = nan_v;
    }
    return s;
}

} // anon

//  Top-level driver
int run_monte_carlo(const std::string& config_path,
                    const std::string& outputs_csv,
                    const std::string& inputs_csv,
                    const std::string& summary_txt)
{
    // Parse the config
    json::Value cfg;
    try {
        cfg = json::parse_file(config_path);
    } catch (const json::ParseError& e) {
        std::fprintf(stderr,
                     "ERROR parsing %s at line %d col %d: %s\n",
                     config_path.c_str(), e.line, e.col, e.what());
        return 1;
    }

    const json::Value& mc = cfg["monte_carlo"];
    if (!mc.isObject()) {
        std::fprintf(stderr,
                     "ERROR: config %s has no 'monte_carlo' section\n",
                     config_path.c_str());
        return 1;
    }

    unsigned long seed = static_cast<unsigned long>(mc["seed"].asNumber(12345));
    int    n_runs      = mc["n_runs"].asInt(10);
    // Default 0 means auto-detect; resolver picks min(hw_concurrency, 8).
    // Explicit user values >= 1 are honored without cap.
    int    raw_workers = mc["n_workers"].asInt(0);
    if (n_runs <= 0) {
        std::fprintf(stderr, "ERROR: n_runs must be > 0 (got %d)\n", n_runs);
        return 1;
    }
    ResolvedWorkers rw = resolve_n_workers(raw_workers);
    int    n_workers   = rw.n;

    // Parse variable list
    const json::Value& vars = mc["variables"];
    if (!vars.isArray() || vars.size() == 0) {
        std::fprintf(stderr, "ERROR: monte_carlo.variables must be a non-empty array\n");
        return 1;
    }

    struct VarSpec {
        std::string path;
        std::unique_ptr<Distribution> dist;
        std::mt19937 rng;
    };
    std::vector<VarSpec> spec_list;
    spec_list.reserve(vars.size());
    for (size_t i = 0; i < vars.size(); ++i) {
        VarSpec vs;
        vs.path = vars[i]["path"].asString();
        if (vs.path.empty()) {
            std::fprintf(stderr, "ERROR: variable %zu has empty path\n", i);
            return 1;
        }
        try {
            vs.dist = make_distribution(vars[i]);
        } catch (const std::exception& e) {
            std::fprintf(stderr, "ERROR for variable '%s': %s\n",
                         vs.path.c_str(), e.what());
            return 1;
        }
        // Per-variable RNG seeded as (top_seed XOR hash(path))
        unsigned long vseed = seed ^ path_hash(vs.path);
        vs.rng.seed(vseed);
        spec_list.push_back(std::move(vs));
    }

    // Output list
    const json::Value& outs = mc["outputs"];
    std::vector<std::string> output_names;
    if (outs.isArray() && outs.size() > 0) {
        for (size_t i = 0; i < outs.size(); ++i) {
            output_names.push_back(outs[i].asString());
        }
    } else {
        // Default outputs
        output_names = { "alt", "dvbi", "fpa", "theta", "psi",
                         "ins_pos_err", "t_end" };
    }
    // Output name validation happens implicitly: unknown names return
    // NaN from the run, which we record in the nan_flag column.

    // ---- Print banner ----
    std::printf("================================================================\n");
    std::printf("    rocket_6dof  --  Monte Carlo sweep\n");
    std::printf("================================================================\n");
    std::printf("Config:        %s\n", config_path.c_str());
    std::printf("Build:         MAX_STAGES=%d\n", Propulsion::MAX_STAGES);
    std::printf("Seed:          %lu\n", seed);
    std::printf("Runs:          %d\n", n_runs);
    std::printf("Workers:       %d  [%s]%s\n", n_workers, rw.source.c_str(),
                n_workers > 1 ? "  (subprocess parallel)" : "  (sequential)");
    std::printf("Variables:     %zu\n", spec_list.size());
    for (const auto& vs : spec_list) {
        std::printf("    %-30s  %s\n", vs.path.c_str(), vs.dist->describe().c_str());
    }
    std::printf("Outputs (%zu): ", output_names.size());
    for (size_t i = 0; i < output_names.size(); ++i) {
        std::printf("%s%s", output_names[i].c_str(),
                    i + 1 < output_names.size() ? ", " : "\n");
    }
    std::printf("================================================================\n\n");

    // ---- Open output files and write headers ----
    std::FILE* fout = std::fopen(outputs_csv.c_str(), "w");
    if (!fout) {
        std::fprintf(stderr, "ERROR: cannot open %s for writing\n", outputs_csv.c_str());
        return 1;
    }
    std::fprintf(fout, "run_index");
    for (const auto& n : output_names) std::fprintf(fout, ",%s", n.c_str());
    std::fprintf(fout, ",nan_flag\n");

    std::FILE* fin = std::fopen(inputs_csv.c_str(), "w");
    if (!fin) {
        std::fprintf(stderr, "ERROR: cannot open %s for writing\n", inputs_csv.c_str());
        std::fclose(fout);
        return 1;
    }
    std::fprintf(fin, "run_index");
    for (const auto& vs : spec_list) std::fprintf(fin, ",%s", vs.path.c_str());
    std::fprintf(fin, "\n");

    // ---- Run sweep ----
    // Phase 1: draw all sample rows up front (deterministic given the
    // per-variable RNG streams seeded above).  Doing this in one pass
    // means the sample matrix is identical regardless of n_workers, so
    // the inputs CSV is bit-identical to a single-process run with the
    // same seed.
    std::vector<std::vector<double>> sample_matrix(n_runs,
        std::vector<double>(spec_list.size()));
    std::vector<std::string> input_paths;
    input_paths.reserve(spec_list.size());
    for (const auto& vs : spec_list) input_paths.push_back(vs.path);

    for (int run = 0; run < n_runs; ++run) {
        for (size_t i = 0; i < spec_list.size(); ++i) {
            sample_matrix[run][i] = spec_list[i].dist->sample(spec_list[i].rng);
        }
    }

    // Phase 2: execute.  When n_workers >= 2, use the subprocess
    // orchestrator; otherwise run sequentially in this process.
    std::vector<std::vector<double>> all_inputs (spec_list.size());
    std::vector<std::vector<double>> all_outputs(output_names.size());
    int nan_runs = 0;
    auto t0 = std::chrono::steady_clock::now();

    // Per-row outputs and nan flags (filled by either branch)
    std::vector<std::vector<double>> row_outputs(n_runs,
        std::vector<double>(output_names.size(), std::nan("")));
    std::vector<bool> row_nan(n_runs, false);

    if (n_workers >= 2) {
        // Parallel: dispatch all sample rows in one chunked call
        ParallelResult pr = run_chunked_parallel(
            config_path, input_paths, sample_matrix, output_names,
            n_workers, "./mission_worker", "/tmp", "MC");
        if (!pr.ok) {
            std::fprintf(stderr, "ERROR: parallel MC sweep failed\n");
            std::fclose(fout);
            std::fclose(fin);
            return 1;
        }
        row_outputs = pr.row_outputs;
        row_nan     = pr.nan_flags;
    } else {
        // Sequential: same logic as the original loop
        for (int run = 0; run < n_runs; ++run) {
            json::Value run_cfg = cfg;
            const auto& samples = sample_matrix[run];
            for (size_t i = 0; i < spec_list.size(); ++i) {
                if (!run_cfg.set_path(spec_list[i].path, samples[i])) {
                    std::fprintf(stderr,
                                 "WARNING run %d: failed to set path '%s' in config\n",
                                 run, spec_list[i].path.c_str());
                }
            }
            SingleRunResult res = run_single_mission(run_cfg, output_names);
            row_outputs[run] = res.outputs;
            row_nan[run]     = res.nan_seen;

            // Sequential progress every ~5%
            if ((run + 1) % std::max(1, n_runs / 20) == 0 || run == n_runs - 1) {
                auto t1 = std::chrono::steady_clock::now();
                double elapsed = std::chrono::duration<double>(t1 - t0).count();
                double per_run = elapsed / (run + 1);
                double eta     = per_run * (n_runs - run - 1);
                std::printf("  run %4d / %d  (%.2f s/run, ETA %.1f s)\n",
                            run + 1, n_runs, per_run, eta);
            }
        }
    }

    // Write CSVs and aggregate stats
    for (int run = 0; run < n_runs; ++run) {
        const auto& samples = sample_matrix[run];
        const auto& outputs = row_outputs[run];
        bool nan_seen = row_nan[run];

        for (size_t i = 0; i < samples.size(); ++i)
            all_inputs[i].push_back(samples[i]);
        for (size_t i = 0; i < outputs.size(); ++i)
            all_outputs[i].push_back(outputs[i]);
        if (nan_seen) nan_runs++;

        std::fprintf(fin, "%d", run);
        for (double s : samples) std::fprintf(fin, ",%.9g", s);
        std::fprintf(fin, "\n");

        std::fprintf(fout, "%d", run);
        for (double v : outputs) std::fprintf(fout, ",%.9g", v);
        std::fprintf(fout, ",%d\n", nan_seen ? 1 : 0);
    }
    std::fflush(fout);
    std::fflush(fin);

    std::fclose(fout);
    std::fclose(fin);

    // ---- Summary statistics ----
    std::FILE* fsum = std::fopen(summary_txt.c_str(), "w");
    if (fsum) {
        std::fprintf(fsum, "rocket_6dof Monte Carlo summary\n");
        std::fprintf(fsum, "================================\n");
        std::fprintf(fsum, "Config: %s\n", config_path.c_str());
        std::fprintf(fsum, "Seed:   %lu\n", seed);
        std::fprintf(fsum, "Runs:   %d (%d with NaN outputs)\n\n", n_runs, nan_runs);

        std::fprintf(fsum, "INPUT VARIABLES\n");
        std::fprintf(fsum, "  %-30s  %12s  %12s  %12s  %12s\n",
                     "path", "mean", "stddev", "min", "max");
        for (size_t i = 0; i < spec_list.size(); ++i) {
            auto s = compute_stats(all_inputs[i]);
            std::fprintf(fsum, "  %-30s  %12.4g  %12.4g  %12.4g  %12.4g\n",
                         spec_list[i].path.c_str(),
                         s.mean, s.stddev, s.min, s.max);
        }

        std::fprintf(fsum, "\nOUTPUT VARIABLES\n");
        std::fprintf(fsum, "  %-30s  %12s  %12s  %12s  %12s  %5s\n",
                     "name", "mean", "stddev", "min", "max", "n_nan");
        for (size_t i = 0; i < output_names.size(); ++i) {
            auto s = compute_stats(all_outputs[i]);
            std::fprintf(fsum, "  %-30s  %12.4g  %12.4g  %12.4g  %12.4g  %5d\n",
                         output_names[i].c_str(),
                         s.mean, s.stddev, s.min, s.max, s.nan_count);
        }
        std::fclose(fsum);
    }

    // Print summary to stdout too
    std::printf("\n================================================================\n");
    std::printf("MONTE CARLO COMPLETE\n");
    std::printf("================================================================\n");
    std::printf("Runs:           %d\n", n_runs);
    std::printf("NaN outputs:    %d\n", nan_runs);
    std::printf("Outputs CSV:    %s\n", outputs_csv.c_str());
    std::printf("Inputs CSV:     %s\n", inputs_csv.c_str());
    std::printf("Summary:        %s\n", summary_txt.c_str());

    std::printf("\nOUTPUT STATISTICS\n");
    std::printf("  %-25s  %12s  %12s  %12s  %12s  %5s\n",
                "name", "mean", "stddev", "min", "max", "n_nan");
    for (size_t i = 0; i < output_names.size(); ++i) {
        auto s = compute_stats(all_outputs[i]);
        std::printf("  %-25s  %12.4g  %12.4g  %12.4g  %12.4g  %5d\n",
                    output_names[i].c_str(),
                    s.mean, s.stddev, s.min, s.max, s.nan_count);
    }
    std::printf("================================================================\n");

    return 0;
}

} // namespace rocket6dof
