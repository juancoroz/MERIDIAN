//  mission_mc.cpp -- Monte Carlo driver binary
//
//  Usage:
//    ./mission_mc                     # uses monte_carlo.json
//    ./mission_mc my_mc.json          # custom config
//    ./mission_mc my_mc.json out.csv inputs.csv summary.txt   # custom paths

#include "mc_runner.h"

#include <cstdio>
#include <string>

int main(int argc, char** argv) {
    // Line-buffer stdout so a piped consumer (the GUI's subprocess
    // wrapper) sees progress lines in real time.  Without this,
    // block-buffered stdout would withhold "[MC] N rows split" until
    // the program exits, leaving the GUI's progress parser unable
    // to determine the total.
    std::setvbuf(stdout, nullptr, _IOLBF, 0);

    std::string config_path  = (argc >= 2) ? argv[1] : "monte_carlo.json";
    std::string outputs_csv  = (argc >= 3) ? argv[2] : "monte_carlo_outputs.csv";
    std::string inputs_csv   = (argc >= 4) ? argv[3] : "monte_carlo_inputs.csv";
    std::string summary_txt  = (argc >= 5) ? argv[4] : "monte_carlo_summary.txt";
    return rocket6dof::run_monte_carlo(config_path, outputs_csv,
                                       inputs_csv, summary_txt);
}
