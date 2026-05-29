//  gui_sweeps.h -- Sweeps tab: MC / Sobol / viewer integration
//
//  Encapsulates the long-lived state for the Sweeps tab so the main
//  loop only needs:
//      static gui::SweepsTabState sweeps;
//      if (ImGui::BeginTabItem("Sweeps")) {
//          gui::render_sweeps_tab(sweeps);
//          ImGui::EndTabItem();
//      }
//
//  The state struct owns three Subprocess instances (one per workflow:
//  MC, Sobol, viewer) and three FilePickers (for choosing MC config,
//  Sobol config, and viewer-output path).  Pickers are shown
//  on-demand.
#ifndef ROCKET6DOF_GUI_SWEEPS_H
#define ROCKET6DOF_GUI_SWEEPS_H

#include "gui_subprocess.h"
#include "gui_file_picker.h"
#include <string>
#include <vector>

namespace gui {

// Parsed Sobol row: one input's S1/ST for a particular output.
struct SobolRow {
    std::string input;
    double      S1     = 0.0;
    double      ST     = 0.0;
    double      S1_lo  = 0.0;
    double      S1_hi  = 0.0;
    double      ST_lo  = 0.0;
    double      ST_hi  = 0.0;
};

struct SweepsTabState {
    // ---- Monte Carlo ----
    std::string  mc_config_path  = "monte_carlo_reduced.json";
    Subprocess   mc_proc;
    std::string  mc_last_status;    // "Done in 14.4s" / "Failed (rc=2)" etc.
    FilePicker   mc_picker;
    // Summary parsed from monte_carlo_summary.txt after success
    std::string  mc_summary_text;

    // ---- Sobol ----
    std::string  sobol_config_path = "monte_carlo_sobol.json";
    Subprocess   sobol_proc;
    std::string  sobol_last_status;
    FilePicker   sobol_picker;
    std::string  sobol_summary_text;

    // ---- Viewer ----
    Subprocess   viewer_proc;
    std::string  viewer_last_status;
    std::string  viewer_html_path;   // "mc_viewer.html" once built
    bool         viewer_built = false;

    // ---- In-GUI results viewer ----
    // Parsed monte_carlo_outputs.csv: column names + column data.
    // Empty until Load is clicked.
    std::vector<std::string>          mc_output_names;
    std::vector<std::vector<double>>  mc_output_data;
    bool                              mc_results_loaded = false;
    int                               mc_selected_output = 0;
    int                               mc_n_bins          = 30;

    // Parsed sobol_indices.csv grouped by output: map output -> rows.
    // For simplicity stored as parallel vectors (output names + their rows).
    std::vector<std::string>            sobol_output_names;
    std::vector<std::vector<SobolRow>>  sobol_by_output;
    bool                                sobol_results_loaded = false;
    int                                 sobol_selected_output = 0;
};

// Render the Sweeps tab.  Call from inside an ImGui::BeginTabItem block.
void render_sweeps_tab(SweepsTabState& s);

}  // namespace gui

#endif
