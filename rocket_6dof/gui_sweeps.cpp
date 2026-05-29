//  gui_sweeps.cpp -- Sweeps tab implementation
//
//  Three workflows in one tab:
//    1. Monte Carlo  -- launches ./mission_mc with a JSON config
//    2. Sobol        -- launches ./mission_sobol with a JSON config
//    3. Build viewer -- launches ./build_mc_viewer to bundle the
//                       latest MC + Sobol CSVs into mc_viewer.html
//
//  Each workflow has its own Subprocess and FilePicker.  The UI for
//  each is a section with a config path field, Browse, Run, status
//  display, and (on success) a "Open results" button.
//
//  Why not run MC/Sobol in-process?  Both binaries use fork()ed
//  workers and were architected as external processes.  Linking them
//  into the GUI would require non-trivial refactoring.  Subprocess
//  is the right boundary.
#include "gui_sweeps.h"
#include "imgui.h"
#include "implot.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>

namespace gui {

namespace {

// Small helpers

// Read a small text file fully into a std::string.  Used to load the
// MC/Sobol summary files after a successful run.  Returns empty on
// failure.
std::string read_file(const std::string& path) {
    std::ifstream f(path);
    if (!f) return "";
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

// Try platform-appropriate launchers in order to open a local file
// (typically an HTML viewer) in the user's default browser.
//
// Order:
//   explorer.exe                -- WSL primary, when PATH includes
//                                  the Windows directories (default
//                                  in WSL's bundled bash).
//   /mnt/c/Windows/explorer.exe -- WSL fallback when PATH is stripped
//                                  (e.g. XRDP-launched terminals,
//                                  some remote SSH sessions).
//   xdg-open                    -- native Linux desktop environments.
//   wslview                     -- last-resort fallback (wslu package).
//
// Returns the name of the launcher that was invoked, or an empty
// string if none was available.  Launchers are run fire-and-forget;
// success here only means dispatch succeeded.
std::string open_in_browser(const std::string& path) {
    auto have = [](const char* cmd) {
        // `command -v` succeeds for both PATH lookup and absolute
        // paths to existing executables.
        std::string c = std::string("command -v ") + cmd + " >/dev/null 2>&1";
        return std::system(c.c_str()) == 0;
    };
    const char* candidates[] = {
        "explorer.exe",
        "/mnt/c/Windows/explorer.exe",
        "xdg-open",
        "wslview",
        nullptr
    };
    for (const char** c = candidates; *c; ++c) {
        if (have(*c)) {
            std::string cmd = std::string("'") + *c + "' '" + path + "' &";
            std::system(cmd.c_str());
            return *c;
        }
    }
    return "";
}

// Render a labeled config-file row: text input + Browse button.
// Returns true if the text field was edited this frame.  The picker
// is drawn separately by the caller.
bool config_path_row(const char* id_suffix,
                     std::string& path,
                     FilePicker& picker,
                     const char* picker_title,
                     const char* tooltip = nullptr) {
    bool edited = false;
    char buf[512];
    std::snprintf(buf, sizeof(buf), "%s", path.c_str());
    ImGui::PushID(id_suffix);

    // Reserve room for the Browse button on the right.
    float avail = ImGui::GetContentRegionAvail().x;
    ImGui::SetNextItemWidth(avail - 90.0f);
    if (ImGui::InputText("##cfg", buf, sizeof(buf))) {
        path = buf;
        edited = true;
    }
    if (tooltip && ImGui::IsItemHovered()) {
        ImGui::SetTooltip("%s", tooltip);
    }
    ImGui::SameLine();
    if (ImGui::Button("Browse...")) {
        picker.open(picker_title, false, path, ".json");
    }
    ImGui::PopID();
    return edited;
}

// Style the Run button green when idle, grey + "Running..." when busy.
// Returns true if clicked while idle.
bool run_button(const char* id, bool busy, const char* idle_label) {
    bool clicked = false;
    ImGui::PushID(id);
    if (busy) {
        ImGui::PushStyleColor(ImGuiCol_Button,        IM_COL32(80,  80,  80, 255));
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, IM_COL32(80,  80,  80, 255));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive,  IM_COL32(80,  80,  80, 255));
        ImGui::BeginDisabled();
    } else {
        ImGui::PushStyleColor(ImGuiCol_Button,        IM_COL32(40, 110, 60, 255));
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, IM_COL32(60, 150, 80, 255));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive,  IM_COL32(30,  90, 50, 255));
    }
    const char* label = busy ? "Running..." : idle_label;
    if (ImGui::Button(label, ImVec2(160, 30))) {
        clicked = true;
    }
    if (busy) ImGui::EndDisabled();
    ImGui::PopStyleColor(3);
    ImGui::PopID();
    return clicked && !busy;
}

// "Cancel" button enabled only while busy; sends SIGTERM.
void cancel_button(const char* id, Subprocess& sp) {
    ImGui::PushID(id);
    bool busy = sp.is_running();
    if (!busy) ImGui::BeginDisabled();
    if (ImGui::Button("Cancel", ImVec2(80, 30))) {
        sp.signal(15);  // SIGTERM
    }
    if (!busy) ImGui::EndDisabled();
    ImGui::PopID();
}

// Parse the subprocess log for progress info emitted by mission_worker.cpp
// and parallel_runner.cpp.  Returns {completed, total}.  If total is
// 0, no "N rows split" line has appeared yet.
//
// Format from workers (one per row completion):
//     PROGRESS row=<N>\n
// Format from parallel_runner (one line near start):
//     [MC] <N> rows split across <W> workers (~<R> rows/worker)
//     [Sobol] <N> rows split across <W> workers (~<R> rows/worker)
//
// The completed count is the total number of PROGRESS rows seen
// across all workers (the GUI is running ONE subprocess that wraps
// ALL workers because parallel_runner.cpp is invoked directly; all
// worker stderrs are merged into our combined stdout+stderr capture).
struct ProgressInfo {
    int completed = 0;
    int total     = 0;
};

ProgressInfo parse_progress(const std::vector<std::string>& lines) {
    ProgressInfo p;
    for (const auto& line : lines) {
        // Worker heartbeat: count occurrences
        if (line.compare(0, 13, "PROGRESS row=") == 0) {
            ++p.completed;
            continue;
        }
        // Parent total line.  Both "[MC]" and "[Sobol]" use the same
        // template; find " rows split across".
        size_t pos = line.find(" rows split across");
        if (pos != std::string::npos) {
            // Walk backwards to find the integer before " rows".
            size_t end = pos;
            size_t start = end;
            while (start > 0 && line[start - 1] >= '0' && line[start - 1] <= '9') {
                --start;
            }
            if (start < end) {
                p.total = std::atoi(line.substr(start, end - start).c_str());
            }
        }
    }
    return p;
}

// Render a status line for a subprocess.  Shows the last stdout line
// while running, and a summary after completion.  Caller passes a
// status string that reflects "last terminal state we observed".
//
// If progress markers are present in the log, also renders a
// progress bar above the status text.
void status_line(const Subprocess& sp, const std::string& post_status) {
    if (sp.is_running()) {
        // Pull the log to look for progress markers.  This is cheap
        // (log_lines() returns a copy of a moderate vector); we do
        // it once per frame.
        auto lines = sp.log_lines();
        ProgressInfo pg = parse_progress(lines);
        if (pg.total > 0) {
            float frac = static_cast<float>(pg.completed) / static_cast<float>(pg.total);
            if (frac > 1.0f) frac = 1.0f;
            if (frac < 0.0f) frac = 0.0f;
            char pbuf[64];
            std::snprintf(pbuf, sizeof(pbuf),
                          "%d / %d  (%.1fs)",
                          pg.completed, pg.total, sp.elapsed_seconds());
            ImGui::ProgressBar(frac, ImVec2(-1.0f, 0.0f), pbuf);
        }

        char buf[64];
        std::snprintf(buf, sizeof(buf), "Running... (%.1fs)  ",
                      sp.elapsed_seconds());
        ImGui::TextColored(ImVec4(1.0f, 0.8f, 0.3f, 1.0f), "%s", buf);
        ImGui::SameLine();
        // Show the latest output line (truncated)
        std::string line = sp.last_line();
        if (line.size() > 120) line = line.substr(0, 117) + "...";
        ImGui::TextDisabled("%s", line.c_str());
    } else if (!post_status.empty()) {
        ImGui::TextDisabled("%s", post_status.c_str());
    }
}

// Workflow sections

void render_mc_section(SweepsTabState& s) {
    ImGui::Text("Monte Carlo");
    ImGui::Spacing();

    config_path_row("mc", s.mc_config_path, s.mc_picker,
                    "Choose MC config",
                    "JSON config defining variables, distributions, "
                    "n_runs, n_workers");
    std::string picked;
    if (s.mc_picker.draw(picked)) {
        s.mc_config_path = picked;
    }

    bool busy = s.mc_proc.is_running();
    if (run_button("run_mc", busy, "Run MC")) {
        // Reset summary; launch
        s.mc_summary_text.clear();
        s.mc_last_status = "Launching...";
        bool ok = s.mc_proc.start({"./mission_mc", s.mc_config_path});
        if (!ok) {
            s.mc_last_status = "Could not start subprocess";
        }
    }
    ImGui::SameLine();
    cancel_button("cancel_mc", s.mc_proc);

    status_line(s.mc_proc, s.mc_last_status);

    // Poll for completion and load summary on success.
    if (s.mc_proc.done()) {
        // done() stays true until the next start().  Use a sticky
        // post-status string so we know whether we've already
        // processed the completion.
        if (s.mc_last_status == "Launching..." ||
            s.mc_last_status.rfind("Running", 0) == 0 ||
            s.mc_last_status == "Could not start subprocess") {
            int rc = s.mc_proc.exit_code();
            char buf[128];
            std::snprintf(buf, sizeof(buf),
                          "%s in %.1fs (rc=%d)",
                          rc == 0 ? "Done" : "Failed",
                          s.mc_proc.elapsed_seconds(), rc);
            s.mc_last_status = buf;
            if (rc == 0) {
                s.mc_summary_text = read_file("monte_carlo_summary.txt");
            }
        }
    }

    // Show the summary in a scrollable child after a successful run.
    if (!s.mc_summary_text.empty()) {
        ImGui::Spacing();
        ImGui::Text("Last MC summary:");
        ImGui::PushStyleColor(ImGuiCol_ChildBg, IM_COL32(20, 20, 25, 255));
        ImGui::BeginChild("mc_summary",
                          ImVec2(0, 140), true,
                          ImGuiWindowFlags_HorizontalScrollbar);
        ImGui::TextUnformatted(s.mc_summary_text.c_str());
        ImGui::EndChild();
        ImGui::PopStyleColor();
    }
}

void render_sobol_section(SweepsTabState& s) {
    ImGui::Text("Sobol sensitivity analysis");
    ImGui::Spacing();

    config_path_row("sobol", s.sobol_config_path, s.sobol_picker,
                    "Choose Sobol config",
                    "JSON config defining input variables, "
                    "n_samples, output expressions");
    std::string picked;
    if (s.sobol_picker.draw(picked)) {
        s.sobol_config_path = picked;
    }

    bool busy = s.sobol_proc.is_running();
    if (run_button("run_sobol", busy, "Run Sobol")) {
        s.sobol_summary_text.clear();
        s.sobol_last_status = "Launching...";
        bool ok = s.sobol_proc.start({"./mission_sobol", s.sobol_config_path});
        if (!ok) s.sobol_last_status = "Could not start subprocess";
    }
    ImGui::SameLine();
    cancel_button("cancel_sobol", s.sobol_proc);

    status_line(s.sobol_proc, s.sobol_last_status);

    if (s.sobol_proc.done()) {
        if (s.sobol_last_status == "Launching..." ||
            s.sobol_last_status.rfind("Running", 0) == 0 ||
            s.sobol_last_status == "Could not start subprocess") {
            int rc = s.sobol_proc.exit_code();
            char buf[128];
            std::snprintf(buf, sizeof(buf),
                          "%s in %.1fs (rc=%d)",
                          rc == 0 ? "Done" : "Failed",
                          s.sobol_proc.elapsed_seconds(), rc);
            s.sobol_last_status = buf;
            if (rc == 0) {
                s.sobol_summary_text = read_file("sobol_summary.txt");
            }
        }
    }

    if (!s.sobol_summary_text.empty()) {
        ImGui::Spacing();
        ImGui::Text("Last Sobol summary:");
        ImGui::PushStyleColor(ImGuiCol_ChildBg, IM_COL32(20, 20, 25, 255));
        ImGui::BeginChild("sobol_summary",
                          ImVec2(0, 140), true,
                          ImGuiWindowFlags_HorizontalScrollbar);
        ImGui::TextUnformatted(s.sobol_summary_text.c_str());
        ImGui::EndChild();
        ImGui::PopStyleColor();
    }
}

void render_viewer_section(SweepsTabState& s) {
    ImGui::Text("Build interactive viewer");
    ImGui::Spacing();

    ImGui::TextDisabled(
        "Bundles the latest monte_carlo_*.csv and sobol_indices*.csv "
        "into a single self-contained HTML file (mc_viewer.html) with "
        "interactive histograms, scatter plots, and sensitivity bars.");
    ImGui::Spacing();

    bool busy = s.viewer_proc.is_running();
    if (run_button("run_viewer", busy, "Build viewer")) {
        s.viewer_last_status = "Launching...";
        s.viewer_built = false;
        bool ok = s.viewer_proc.start({"./build_mc_viewer"});
        if (!ok) s.viewer_last_status = "Could not start subprocess";
    }
    ImGui::SameLine();

    // Open in browser button: only enabled if viewer was built.
    if (!s.viewer_built) ImGui::BeginDisabled();
    if (ImGui::Button("Open in browser", ImVec2(160, 30))) {
        // Fire-and-forget dispatch via the first available launcher.
        // On WSL: explorer.exe wins (handled by WSL kernel interop).
        // On native Linux: xdg-open is used.
        std::string used = open_in_browser(s.viewer_html_path);
        if (!used.empty()) {
            s.viewer_last_status = "Opened " + s.viewer_html_path
                                 + " via " + used;
        } else {
            s.viewer_last_status =
                "No launcher found (need xdg-open on Linux, "
                "or run from WSL with explorer.exe on PATH)";
        }
    }
    if (!s.viewer_built) ImGui::EndDisabled();

    status_line(s.viewer_proc, s.viewer_last_status);

    if (s.viewer_proc.done()) {
        if (s.viewer_last_status == "Launching..." ||
            s.viewer_last_status.rfind("Running", 0) == 0 ||
            s.viewer_last_status == "Could not start subprocess") {
            int rc = s.viewer_proc.exit_code();
            char buf[128];
            std::snprintf(buf, sizeof(buf),
                          "%s in %.1fs (rc=%d)",
                          rc == 0 ? "Done" : "Failed",
                          s.viewer_proc.elapsed_seconds(), rc);
            s.viewer_last_status = buf;
            if (rc == 0) {
                s.viewer_html_path = "mc_viewer.html";
                s.viewer_built = true;
            }
        }
    }
}

// In-GUI results viewer: MC histograms + Sobol bar charts (ImPlot)

// Tokenize one CSV line, honoring commas inside double-quoted fields.
// Used by both parsers.  Trims surrounding whitespace from fields.
std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    bool in_quotes = false;
    for (char c : line) {
        if (c == '"') { in_quotes = !in_quotes; continue; }
        if (c == ',' && !in_quotes) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    out.push_back(cur);
    return out;
}

// Load monte_carlo_outputs.csv into column-major form.  Skips the
// "run_index" and "nan_flag" bookkeeping columns (the user wants to
// see histograms of the actual outputs).  Returns false on parse
// failure.
bool load_mc_outputs(const std::string& path,
                     std::vector<std::string>& names,
                     std::vector<std::vector<double>>& cols)
{
    names.clear();
    cols.clear();
    std::ifstream f(path);
    if (!f) return false;
    std::string header_line;
    if (!std::getline(f, header_line)) return false;
    auto all_headers = split_csv_line(header_line);
    if (all_headers.empty()) return false;

    // Decide which columns to keep: drop run_index, nan_flag.
    std::vector<int> keep_idx;
    for (size_t i = 0; i < all_headers.size(); ++i) {
        if (all_headers[i] == "run_index") continue;
        if (all_headers[i] == "nan_flag")  continue;
        keep_idx.push_back(static_cast<int>(i));
        names.push_back(all_headers[i]);
    }
    cols.resize(names.size());

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        auto cells = split_csv_line(line);
        if (cells.size() != all_headers.size()) continue;
        for (size_t k = 0; k < keep_idx.size(); ++k) {
            cols[k].push_back(std::strtod(cells[keep_idx[k]].c_str(),
                                          nullptr));
        }
    }
    return !cols.empty() && !cols[0].empty();
}

// Load sobol_indices.csv grouped by output name.  Each row is
// (output, input, S1, ST, S1_lo, S1_hi, ST_lo, ST_hi, mean, variance).
// We discard the mean/variance for plotting purposes.
bool load_sobol_indices(const std::string& path,
                        std::vector<std::string>& output_names,
                        std::vector<std::vector<SobolRow>>& by_output)
{
    output_names.clear();
    by_output.clear();
    std::ifstream f(path);
    if (!f) return false;

    std::string header_line;
    if (!std::getline(f, header_line)) return false;
    auto headers = split_csv_line(header_line);
    // Column lookup: figure out where each field is.
    auto col_idx = [&](const char* name) -> int {
        for (size_t i = 0; i < headers.size(); ++i) {
            if (headers[i] == name) return static_cast<int>(i);
        }
        return -1;
    };
    int i_out   = col_idx("output");
    int i_in    = col_idx("input");
    int i_S1    = col_idx("S1");
    int i_ST    = col_idx("ST");
    int i_S1lo  = col_idx("S1_lo");
    int i_S1hi  = col_idx("S1_hi");
    int i_STlo  = col_idx("ST_lo");
    int i_SThi  = col_idx("ST_hi");
    if (i_out < 0 || i_in < 0 || i_S1 < 0 || i_ST < 0) return false;

    // Track output -> index for grouping
    std::unordered_map<std::string, size_t> out_index;

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        auto cells = split_csv_line(line);
        if (static_cast<int>(cells.size()) <= i_ST) continue;
        const std::string& out_name = cells[i_out];

        auto it = out_index.find(out_name);
        if (it == out_index.end()) {
            out_index[out_name] = output_names.size();
            output_names.push_back(out_name);
            by_output.push_back({});
            it = out_index.find(out_name);
        }

        SobolRow r;
        r.input = cells[i_in];
        r.S1    = std::strtod(cells[i_S1].c_str(),   nullptr);
        r.ST    = std::strtod(cells[i_ST].c_str(),   nullptr);
        r.S1_lo = (i_S1lo >= 0) ? std::strtod(cells[i_S1lo].c_str(), nullptr) : r.S1;
        r.S1_hi = (i_S1hi >= 0) ? std::strtod(cells[i_S1hi].c_str(), nullptr) : r.S1;
        r.ST_lo = (i_STlo >= 0) ? std::strtod(cells[i_STlo].c_str(), nullptr) : r.ST;
        r.ST_hi = (i_SThi >= 0) ? std::strtod(cells[i_SThi].c_str(), nullptr) : r.ST;
        by_output[it->second].push_back(r);
    }

    return !output_names.empty();
}

// Compute and render a histogram using ImPlot.  We pre-bin the data
// rather than using PlotHistogram's built-in binning so we have control
// over the bin range and count (ImPlot's defaults can be quirky for
// constant or near-constant outputs).
void plot_histogram(const std::string& label,
                    const std::vector<double>& data,
                    int n_bins)
{
    if (data.empty()) {
        ImGui::TextDisabled("No data for %s", label.c_str());
        return;
    }
    if (n_bins < 4) n_bins = 4;
    if (n_bins > 200) n_bins = 200;

    double mn =  std::numeric_limits<double>::infinity();
    double mx = -std::numeric_limits<double>::infinity();
    int n_valid = 0;
    for (double v : data) {
        if (!std::isfinite(v)) continue;
        if (v < mn) mn = v;
        if (v > mx) mx = v;
        ++n_valid;
    }
    if (n_valid == 0 || mn == mx) {
        ImGui::TextDisabled(
            "%s: %s",
            label.c_str(),
            n_valid == 0 ? "no finite values" :
                           "all values identical (no histogram)");
        return;
    }

    // Bin the data
    std::vector<int> counts(n_bins, 0);
    const double w = (mx - mn) / n_bins;
    for (double v : data) {
        if (!std::isfinite(v)) continue;
        int b = static_cast<int>((v - mn) / w);
        if (b < 0) b = 0;
        if (b >= n_bins) b = n_bins - 1;
        ++counts[b];
    }

    std::vector<double> centers(n_bins), heights(n_bins);
    for (int i = 0; i < n_bins; ++i) {
        centers[i] = mn + (i + 0.5) * w;
        heights[i] = static_cast<double>(counts[i]);
    }

    if (ImPlot::BeginPlot(("Histogram of " + label).c_str(),
                          ImVec2(-1, 240.0f))) {
        ImPlot::SetupAxes(label.c_str(), "count");
        ImPlot::PlotBars("count", centers.data(), heights.data(),
                         n_bins, w * 0.9);
        ImPlot::EndPlot();
    }

    // Compute and show basic stats
    double sum = 0, sumsq = 0;
    for (double v : data) if (std::isfinite(v)) { sum += v; sumsq += v*v; }
    double mean = sum / n_valid;
    double var  = (sumsq / n_valid) - mean * mean;
    if (var < 0) var = 0;
    double stddev = std::sqrt(var);
    ImGui::Text("n=%d  mean=%.4g  stddev=%.4g  range=[%.4g, %.4g]",
                n_valid, mean, stddev, mn, mx);
}

// Render the Sobol S1 + ST bar chart for one output's rows.
void plot_sobol_for_output(const std::string& output_name,
                           const std::vector<SobolRow>& rows)
{
    if (rows.empty()) {
        ImGui::TextDisabled("No Sobol data for %s", output_name.c_str());
        return;
    }

    const int N = static_cast<int>(rows.size());
    std::vector<double> xs(N), s1(N), st(N), s1_lo(N), s1_hi(N);
    std::vector<const char*> labels(N);
    // ImPlot wants C-string pointers that live for the duration of
    // BeginPlot/EndPlot.  We can borrow rows[i].input.c_str() directly.
    for (int i = 0; i < N; ++i) {
        xs[i]    = i;
        s1[i]    = rows[i].S1;
        st[i]    = rows[i].ST;
        // Convert (S1_lo, S1_hi) to symmetric-ish error bars centered
        // on S1.  Some estimators can give asymmetric CIs; ImPlot's
        // PlotErrorBars takes neg/pos magnitudes (always positive).
        s1_lo[i] = rows[i].S1 - rows[i].S1_lo;  // neg magnitude
        s1_hi[i] = rows[i].S1_hi - rows[i].S1;  // pos magnitude
        // Clamp to non-negative (rare for floating-point noise to give
        // S1_lo > S1 etc.).
        if (s1_lo[i] < 0) s1_lo[i] = 0;
        if (s1_hi[i] < 0) s1_hi[i] = 0;
        labels[i] = rows[i].input.c_str();
    }

    if (ImPlot::BeginPlot(("Sobol indices for " + output_name).c_str(),
                          ImVec2(-1, 320.0f))) {
        ImPlot::SetupAxes("input", "index value");
        ImPlot::SetupAxisTicks(ImAxis_X1, xs.data(), N, labels.data());
        // Tick labels use ImPlot's default wrapping for long names;
        // no manual rotation is applied.

        const double bar_w = 0.35;
        // S1 bars centered slightly left, ST slightly right, so they
        // sit side by side rather than overlapping.
        std::vector<double> xs_s1(N), xs_st(N);
        for (int i = 0; i < N; ++i) {
            xs_s1[i] = xs[i] - 0.2;
            xs_st[i] = xs[i] + 0.2;
        }
        ImPlot::PlotBars("S1", xs_s1.data(), s1.data(), N, bar_w);
        ImPlot::PlotBars("ST", xs_st.data(), st.data(), N, bar_w);
        // Error bars on S1 (the most-used index; ST CIs are visually
        // similar and add clutter)
        ImPlot::PlotErrorBars("S1", xs_s1.data(), s1.data(),
                              s1_lo.data(), s1_hi.data(), N);
        ImPlot::EndPlot();
    }
}

void render_results_viewer_section(SweepsTabState& s) {
    ImGui::Text("Results viewer (in-GUI)");
    ImGui::Spacing();
    ImGui::TextDisabled(
        "Read MC and Sobol CSV outputs directly and plot them inline.  "
        "Click \"Reload results\" after a fresh sweep completes.");
    ImGui::Spacing();

    // Reload buttons
    if (ImGui::Button("Reload MC results", ImVec2(170, 30))) {
        s.mc_results_loaded = load_mc_outputs(
            "monte_carlo_outputs.csv",
            s.mc_output_names, s.mc_output_data);
        if (!s.mc_results_loaded) {
            s.mc_output_names.clear();
            s.mc_output_data.clear();
        }
        // Pick the first output that has variance (not all identical
        // values).  This makes the post-load default a USEFUL plot
        // rather than the boring t_end-is-constant case.
        s.mc_selected_output = 0;
        for (size_t i = 0; i < s.mc_output_data.size(); ++i) {
            const auto& col = s.mc_output_data[i];
            if (col.size() < 2) continue;
            double a = col[0], b = col[0];
            for (double v : col) { if (v < a) a = v; if (v > b) b = v; }
            if (b > a) { s.mc_selected_output = static_cast<int>(i); break; }
        }
    }
    ImGui::SameLine();
    if (ImGui::Button("Reload Sobol results", ImVec2(190, 30))) {
        // Try sobol_indices.csv first; fall back to other commonly-
        // produced filenames if the user named them differently.
        const char* candidates[] = {
            "sobol_indices.csv",
            "sobol_indices_reduced.csv",
        };
        s.sobol_results_loaded = false;
        for (const char* path : candidates) {
            if (load_sobol_indices(path,
                                   s.sobol_output_names,
                                   s.sobol_by_output)) {
                s.sobol_results_loaded = true;
                break;
            }
        }
        if (!s.sobol_results_loaded) {
            s.sobol_output_names.clear();
            s.sobol_by_output.clear();
        }
        // Pick the first output whose Sobol rows aren't all zero.
        // Same rationale as MC -- avoid a useless first-glance plot.
        s.sobol_selected_output = 0;
        for (size_t i = 0; i < s.sobol_by_output.size(); ++i) {
            const auto& rows = s.sobol_by_output[i];
            bool nonzero = false;
            for (const auto& r : rows) {
                if (std::fabs(r.S1) > 1e-9 || std::fabs(r.ST) > 1e-9) {
                    nonzero = true; break;
                }
            }
            if (nonzero) { s.sobol_selected_output = static_cast<int>(i); break; }
        }
    }
    ImGui::SameLine();
    ImGui::TextDisabled(
        "MC: %s   Sobol: %s",
        s.mc_results_loaded ? "loaded" : "not loaded",
        s.sobol_results_loaded ? "loaded" : "not loaded");

    ImGui::Spacing();

    // ---- MC histogram subsection ----
    if (s.mc_results_loaded && !s.mc_output_names.empty()) {
        ImGui::PushID("mc_hist");
        ImGui::Text("MC histogram");
        ImGui::SetNextItemWidth(220.0f);
        // Output dropdown
        if (ImGui::BeginCombo("Output",
                              s.mc_output_names[s.mc_selected_output].c_str())) {
            for (int i = 0; i < static_cast<int>(s.mc_output_names.size()); ++i) {
                bool is_sel = (i == s.mc_selected_output);
                if (ImGui::Selectable(s.mc_output_names[i].c_str(), is_sel)) {
                    s.mc_selected_output = i;
                }
            }
            ImGui::EndCombo();
        }
        ImGui::SameLine();
        ImGui::SetNextItemWidth(180.0f);
        ImGui::SliderInt("Bins", &s.mc_n_bins, 4, 100);

        plot_histogram(s.mc_output_names[s.mc_selected_output],
                       s.mc_output_data[s.mc_selected_output],
                       s.mc_n_bins);
        ImGui::PopID();
    }

    ImGui::Spacing();

    // ---- Sobol bars subsection ----
    if (s.sobol_results_loaded && !s.sobol_output_names.empty()) {
        ImGui::PushID("sobol_bars");
        ImGui::Text("Sobol sensitivity");
        ImGui::SetNextItemWidth(220.0f);
        if (ImGui::BeginCombo("Output",
                              s.sobol_output_names[s.sobol_selected_output].c_str())) {
            for (int i = 0; i < static_cast<int>(s.sobol_output_names.size()); ++i) {
                bool is_sel = (i == s.sobol_selected_output);
                if (ImGui::Selectable(s.sobol_output_names[i].c_str(), is_sel)) {
                    s.sobol_selected_output = i;
                }
            }
            ImGui::EndCombo();
        }
        plot_sobol_for_output(
            s.sobol_output_names[s.sobol_selected_output],
            s.sobol_by_output[s.sobol_selected_output]);
        ImGui::PopID();
    }
}

}  // anon

void render_sweeps_tab(SweepsTabState& s) {
    ImGui::TextDisabled(
        "Run Monte Carlo / Sobol sweeps as external subprocesses.  "
        "Each section is independent; you can run MC and Sobol "
        "back-to-back or interleave.  When both have completed, "
        "build the interactive viewer to see the results.");
    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    render_mc_section(s);

    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    render_sobol_section(s);

    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    render_viewer_section(s);

    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    render_results_viewer_section(s);
}

}  // namespace gui
