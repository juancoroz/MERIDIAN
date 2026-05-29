//  mission_sobol.cpp  --  Sobol sensitivity-analysis binary
//
//  Usage:
//    ./mission_sobol                            # uses monte_carlo.json
//    ./mission_sobol my_mc.json
//    ./mission_sobol my_mc.json indices.csv summary.txt
//    ./mission_sobol --suggest-n my_mc.json     # pilot mode (advisory)
//
//  In --suggest-n mode:
//    Runs the configured Sobol analysis at a small pilot N (default 64)
//    with bootstrap enabled, observes the worst-case ST CI half-width,
//    and prints a recommended production N to hit a target half-width.
//    The pilot N is forced via env var SOBOL_PILOT_N (default 64) and
//    the target via sobol.target_ci_width in the config (default 0.05).
//    No production run is performed.

#include "sobol_runner.h"

#include <cstdio>
#include <cstring>
#include <string>

int main(int argc, char** argv) {
    // Line-buffer stdout so a piped consumer (the GUI's subprocess
    // wrapper) sees progress lines in real time.  See mission_mc.cpp
    // for the rationale.
    std::setvbuf(stdout, nullptr, _IOLBF, 0);

    bool suggest_n = false;

    // Strip --suggest-n if present and shift remaining args
    int argi = 1;
    if (argi < argc && std::strcmp(argv[argi], "--suggest-n") == 0) {
        suggest_n = true;
        argi++;
    }

    std::string config_path  = (argi < argc) ? argv[argi++] : "monte_carlo.json";
    std::string indices_csv  = (argi < argc) ? argv[argi++] : "sobol_indices.csv";
    std::string summary_txt  = (argi < argc) ? argv[argi++] : "sobol_summary.txt";

    if (suggest_n) {
        return rocket6dof::suggest_n(config_path);
    }
    return rocket6dof::run_sobol(config_path, indices_csv, summary_txt);
}
