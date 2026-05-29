//  mission_worker.cpp -- Headless sim worker for parallel sweeps
//
//  Reads a config JSON, an output-names list, and a CSV of input
//  overrides (one row per sim).  Runs each row through
//  run_single_mission() and writes the output values to stdout (or
//  the file path passed as the third argument).
//
//  This is the building block for parallel MC/Sobol sweeps when
//  the orchestrator (typically a Python script) wants to fan out
//  across processes.  Each worker has its own copy of all kernel
//  globals (osk::State::t etc.), avoiding the kernel's static-member
//  constraint that prevents threading.
//
//  Usage:
//    ./mission_worker <config.json> <overrides.csv> [<output.csv>]
//
//  overrides.csv format (header row required):
//    row_id,path1,path2,...,pathK
//    0,500.3,300.1,...,12345.0
//    1,498.7,301.5,...,67890.0
//    ...
//  Each row sets the named JSON paths to those values via set_path()
//  before running the sim.
//
//  Output CSV columns:
//    row_id,<output_names>,nan_flag

#include "json.h"
#include "sim_runner.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : line) {
        if (c == ',') { out.push_back(cur); cur.clear(); }
        else if (c == '\r') {}  // strip CR
        else cur.push_back(c);
    }
    out.push_back(cur);
    return out;
}

}  // anon

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr,
            "Usage: %s <config.json> <overrides.csv> [<output.csv>]\n", argv[0]);
        return 2;
    }
    std::string config_path    = argv[1];
    std::string overrides_path = argv[2];
    std::string output_path    = (argc >= 4) ? argv[3] : "";

    // Parse config
    rocket6dof::json::Value cfg;
    try {
        cfg = rocket6dof::json::parse_file(config_path);
    } catch (const rocket6dof::json::ParseError& e) {
        std::fprintf(stderr, "ERROR parsing %s at line %d col %d: %s\n",
                     config_path.c_str(), e.line, e.col, e.what());
        return 1;
    }

    // Determine output names from monte_carlo.outputs (fall back to defaults)
    std::vector<std::string> output_names;
    const auto& mc = cfg["monte_carlo"];
    if (mc.isObject() && mc["outputs"].isArray() && mc["outputs"].size() > 0) {
        for (size_t i = 0; i < mc["outputs"].size(); ++i) {
            output_names.push_back(mc["outputs"][i].asString());
        }
    } else {
        output_names = { "alt", "dvbi", "fpa", "theta", "psi",
                         "ins_pos_err", "t_end" };
    }

    // Read overrides CSV
    std::ifstream fin(overrides_path);
    if (!fin) {
        std::fprintf(stderr, "ERROR: cannot open %s\n", overrides_path.c_str());
        return 1;
    }
    std::string header_line;
    if (!std::getline(fin, header_line)) {
        std::fprintf(stderr, "ERROR: empty overrides file\n");
        return 1;
    }
    auto header = split_csv_line(header_line);
    if (header.size() < 2 || header[0] != "row_id") {
        std::fprintf(stderr, "ERROR: overrides header must start with row_id\n");
        return 1;
    }
    std::vector<std::string> input_paths(header.begin() + 1, header.end());

    // Output stream
    std::FILE* fout = stdout;
    if (!output_path.empty()) {
        fout = std::fopen(output_path.c_str(), "w");
        if (!fout) {
            std::fprintf(stderr, "ERROR: cannot open %s for writing\n", output_path.c_str());
            return 1;
        }
    }

    // Write output header
    std::fprintf(fout, "row_id");
    for (const auto& n : output_names) std::fprintf(fout, ",%s", n.c_str());
    std::fprintf(fout, ",nan_flag\n");

    // Process each row
    std::string line;
    int rows_processed = 0;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        auto cells = split_csv_line(line);
        if (cells.size() != header.size()) {
            std::fprintf(stderr, "WARN: row has %zu cells, expected %zu (skipping)\n",
                         cells.size(), header.size());
            continue;
        }
        // Apply overrides to a fresh copy of the config
        rocket6dof::json::Value run_cfg = cfg;
        for (size_t i = 0; i < input_paths.size(); ++i) {
            double v = std::strtod(cells[i + 1].c_str(), nullptr);
            run_cfg.set_path(input_paths[i], v);
        }

        // Run the sim
        auto res = rocket6dof::run_single_mission(run_cfg, output_names);
        std::fprintf(fout, "%s", cells[0].c_str());
        for (double v : res.outputs) std::fprintf(fout, ",%.17g", v);
        std::fprintf(fout, ",%d\n", res.nan_seen ? 1 : 0);
        std::fflush(fout);
        rows_processed++;

        // Heartbeat on stderr so a parent process (the GUI's
        // subprocess wrapper) can track progress.  Workers' stdout
        // is redirected to /dev/null by parallel_runner.cpp to
        // prevent interleaving; stderr stays attached.  The exact
        // format is a stable contract with gui_sweeps.cpp:
        //   PROGRESS row=<N>\n
        // The parser counts these lines across all workers and
        // displays N_completed / N_total.
        std::fprintf(stderr, "PROGRESS row=%d\n", rows_processed);
        std::fflush(stderr);
    }

    if (!output_path.empty()) std::fclose(fout);
    std::fprintf(stderr, "worker: %d rows processed\n", rows_processed);
    return 0;
}
