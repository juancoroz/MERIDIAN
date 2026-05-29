//  gui_config.cpp -- Config editor and CSV loader for mission_gui
//
//  Renders editor widgets over a JSON config tree (rocket6dof::json::Value)
//  using Dear ImGui.  Two tabs:
//    Basic    -- curated subset of ~30 parameters
//    Advanced -- full structured tree, every numeric leaf editable
//
//  Edits are pushed back into the live tree via Value::set_path() so
//  changes flow straight to the run_mission() path.
//
//  Limits this implementation knows about:
//    - The Value tree is in-place editable for numeric paths only.
//      Boolean and string leaves work through set() helpers.
//    - Array indices are not user-extendable in the GUI; if the
//      config has 3 stages, the GUI shows 3 stages.
//    - This file does NOT validate ranges -- the underlying simulator
//      already validates at load time, so bad values just produce
//      visible run failures (rc != 0).

#include "imgui.h"
#include "implot.h"
#include "json.h"
#include "mission_progress.h"
#include "gui_file_picker.h"
#include "gui_plot_panel.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace gui {

using rocket6dof::json::Value;

// File I/O

bool load_config_file(const std::string& path, Value& out) {
    try {
        out = rocket6dof::json::parse_file(path);
        return true;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "load_config_file('%s') failed: %s\n",
                     path.c_str(), e.what());
        return false;
    }
}

bool save_config_file(const std::string& path, const Value& cfg) {
    std::ofstream f(path);
    if (!f) {
        std::fprintf(stderr, "save_config_file('%s'): cannot open\n",
                     path.c_str());
        return false;
    }
    f << cfg.dump(0);
    f << "\n";
    return f.good();
}

// Small editor primitives
//
// These all take a path string and a reference to the root config.
// They display a labeled widget and, if edited, write the new value
// back via cfg.set_path() (for numbers) or via direct map mutation
// (for bools/strings -- but the simulator's mission_config.cpp uses
// scalars heavily, so set_path covers the common case).

// Read a numeric leaf by dotted path; returns def if missing.  Mirrors
// the json::Value chained-access pattern, but only supports
// dot-and-bracket paths that set_path also supports.
static double read_number_path(const Value& cfg, const std::string& path,
                               double def = 0.0) {
    // Walk the path manually -- json::Value doesn't expose a "get_path".
    // Tokenize on '.', then handle bracket-indexed segments.
    const Value* cur = &cfg;
    std::string seg;
    auto descend = [&](const std::string& key) -> bool {
        if (!cur->isObject()) return false;
        const Value& next = (*cur)[key];
        if (!next.exists()) return false;
        cur = &next;
        return true;
    };
    auto descend_idx = [&](size_t idx) -> bool {
        if (!cur->isArray() || idx >= cur->size()) return false;
        cur = &(*cur)[idx];
        return true;
    };

    size_t i = 0, n = path.size();
    while (i < n) {
        seg.clear();
        // Read identifier up to '.', '[' or end
        while (i < n && path[i] != '.' && path[i] != '[') {
            seg.push_back(path[i++]);
        }
        if (!descend(seg)) return def;
        // Bracket indices: [0][1]...
        while (i < n && path[i] == '[') {
            ++i;
            size_t idx = 0;
            bool digits = false;
            while (i < n && path[i] >= '0' && path[i] <= '9') {
                idx = idx * 10 + static_cast<size_t>(path[i] - '0');
                digits = true;
                ++i;
            }
            if (!digits || i >= n || path[i] != ']') return def;
            ++i;  // skip ']'
            if (!descend_idx(idx)) return def;
        }
        if (i < n && path[i] == '.') ++i;  // skip dot
    }
    return cur->isNumber() ? cur->asNumber() : def;
}

// Same walk, string leaf.  Returns def if the path is missing or the
// leaf isn't a string.  Used by Basic-tab widgets that read string
// fields like aerodynamics.aero_file.
static std::string read_string_path(const Value& cfg,
                                    const std::string& path,
                                    const std::string& def = "") {
    const Value* cur = &cfg;
    std::string seg;
    auto descend = [&](const std::string& key) -> bool {
        if (!cur->isObject()) return false;
        const Value& next = (*cur)[key];
        if (!next.exists()) return false;
        cur = &next;
        return true;
    };
    auto descend_idx = [&](size_t idx) -> bool {
        if (!cur->isArray() || idx >= cur->size()) return false;
        cur = &(*cur)[idx];
        return true;
    };

    size_t i = 0, n = path.size();
    while (i < n) {
        seg.clear();
        while (i < n && path[i] != '.' && path[i] != '[') {
            seg.push_back(path[i++]);
        }
        if (!descend(seg)) return def;
        while (i < n && path[i] == '[') {
            ++i;
            size_t idx = 0;
            bool digits = false;
            while (i < n && path[i] >= '0' && path[i] <= '9') {
                idx = idx * 10 + static_cast<size_t>(path[i] - '0');
                digits = true;
                ++i;
            }
            if (!digits || i >= n || path[i] != ']') return def;
            ++i;
            if (!descend_idx(idx)) return def;
        }
        if (i < n && path[i] == '.') ++i;
    }
    return cur->isString() ? cur->asString() : def;
}

// Float-edit widget bound to a config path.  Uses InputDouble.
static void edit_double(Value& cfg, const char* label,
                        const std::string& path,
                        const char* fmt = "%.6g",
                        const char* tooltip = nullptr) {
    double v = read_number_path(cfg, path);
    ImGui::PushID(path.c_str());
    if (ImGui::InputDouble(label, &v, 0.0, 0.0, fmt)) {
        cfg.set_path(path, v);
    }
    if (tooltip && ImGui::IsItemHovered()) {
        ImGui::SetTooltip("%s", tooltip);
    }
    ImGui::PopID();
}

// Slider edit, useful for bounded params like FPA in [0, 90]
static void edit_slider(Value& cfg, const char* label,
                        const std::string& path,
                        double lo, double hi,
                        const char* fmt = "%.3f",
                        const char* tooltip = nullptr) {
    double v = read_number_path(cfg, path);
    float vf = static_cast<float>(v);
    ImGui::PushID(path.c_str());
    if (ImGui::SliderFloat(label, &vf,
                           static_cast<float>(lo),
                           static_cast<float>(hi), fmt)) {
        cfg.set_path(path, static_cast<double>(vf));
    }
    if (tooltip && ImGui::IsItemHovered()) {
        ImGui::SetTooltip("%s", tooltip);
    }
    ImGui::PopID();
}

// Read-only display of a string field (e.g. aero_file path).  Currently
// unused -- kept here for future Basic-tab fields like an aero deck
// path display.  Marked __attribute__((unused)) so -Wall doesn't fire.
__attribute__((unused))
static void show_string(const Value& cfg, const char* label,
                        const std::string& path) {
    // For strings we need to walk manually -- read_number_path returns
    // 0 for non-numbers, but here we want the actual string.
    // For an MVP, just show that the field exists and let Advanced
    // tab handle the (rare) case of editing.
    const Value* cur = &cfg;
    std::string seg;
    size_t i = 0, n = path.size();
    while (i < n) {
        seg.clear();
        while (i < n && path[i] != '.' && path[i] != '[') seg.push_back(path[i++]);
        if (!cur->isObject() || !((*cur)[seg]).exists()) {
            ImGui::TextDisabled("%s: <missing>", label);
            return;
        }
        cur = &(*cur)[seg];
        if (i < n && path[i] == '.') ++i;
    }
    ImGui::Text("%s: %s", label, cur->asString().c_str());
}

// Basic tab: curated parameter subset
//
// The fields here are the ones most engineers reach for first when
// setting up a new mission scenario.  They map directly onto the
// JSON paths that mission_config.cpp consumes.

bool render_basic_tab(Value& cfg) {
    bool run_clicked = false;

    // Brighter header bars (see render_advanced_tab for rationale).
    ImGui::PushStyleColor(ImGuiCol_Header,        IM_COL32( 70,  90, 130, 255));
    ImGui::PushStyleColor(ImGuiCol_HeaderHovered, IM_COL32( 90, 120, 170, 255));
    ImGui::PushStyleColor(ImGuiCol_HeaderActive,  IM_COL32(110, 140, 200, 255));

    if (ImGui::CollapsingHeader("Launch state", ImGuiTreeNodeFlags_DefaultOpen)) {
        edit_slider(cfg, "Latitude [deg]",   "launch.lat",     -90.0,  90.0, "%.2f",
                    "Launch latitude in geographic frame");
        edit_slider(cfg, "Longitude [deg]",  "launch.lon",    -180.0, 180.0, "%.2f",
                    "Launch longitude in geographic frame");
        edit_slider(cfg, "FPA [deg]",        "launch.thtbdx0", 0.0,   90.0,  "%.2f",
                    "Initial flight-path angle from local horizontal");
        edit_slider(cfg, "Heading [deg]",    "launch.psibdx0", 0.0,  360.0,  "%.2f",
                    "Initial heading clockwise from north");
        edit_double(cfg, "Alt [m]",          "launch.alt0",   "%.2f");
    }

    if (ImGui::CollapsingHeader("Propulsion (stage 1)", ImGuiTreeNodeFlags_DefaultOpen)) {
        edit_double(cfg, "Prop mass [kg]",       "propulsion.stages[0].fmass0",         "%.1f",
                    "Stage-1 propellant mass at ignition");
        edit_double(cfg, "Isp [s]",              "propulsion.stages[0].spi",            "%.1f",
                    "Specific impulse, defines exhaust velocity = Isp * g0");
        edit_double(cfg, "Fuel flow [kg/s]",     "propulsion.stages[0].fuel_flow_rate", "%.3f",
                    "Mass flow rate; thrust = flow * Isp * g0");
        edit_double(cfg, "Burn time max [s]",    "propulsion.stages[0].t_burn_max",     "%.1f");
        edit_double(cfg, "Vehicle dry mass [kg]","propulsion.vmass0",                   "%.1f",
                    "Dry mass of the integrated vehicle (excl. propellant)");
    }

    if (ImGui::CollapsingHeader("Aerodynamics")) {
        // Aero deck file.  The simulator auto-detects OSK or Missile
        // DATCOM format by extension (.asc = DATCOM, anything else =
        // OSK).  Edit in place via the text field, or click Browse
        // to open a file picker.  The picker filter accepts both
        // .asc (DATCOM) and .txt (OSK) so users can switch formats
        // without reconfiguring anything else.
        //
        // The static FilePicker outlives the function call; it
        // remembers the last directory the user browsed to, which is
        // exactly the right behavior for a tool that lets you flip
        // between several aero decks during a tuning session.
        static gui::FilePicker aero_picker;

        // Editable text field for the path
        char path_buf[256];
        std::string current = read_string_path(cfg, "aerodynamics.aero_file");
        std::snprintf(path_buf, sizeof(path_buf), "%s", current.c_str());
        ImGui::PushID("aero_file");
        // Use the full available width minus the Browse button (~70 px)
        float avail_w = ImGui::GetContentRegionAvail().x;
        ImGui::SetNextItemWidth(avail_w - 75.0f);
        if (ImGui::InputText("##aero_file", path_buf, sizeof(path_buf))) {
            cfg.set_path("aerodynamics.aero_file", std::string(path_buf));
        }
        ImGui::SameLine();
        if (ImGui::Button("Browse...")) {
            aero_picker.open("Choose aero deck", false, current,
                             ".asc,.txt");
        }
        ImGui::SameLine();
        ImGui::Text("Aero deck");
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip(
                "Path to aero deck file.  Auto-detects format:\n"
                "  .asc  -- Missile DATCOM format\n"
                "  other -- OSK Table1/Table2 format (e.g. aero.txt)");
        }
        ImGui::PopID();

        // Draw the picker (cheap when closed); on confirm, write the
        // chosen path back into the config.
        std::string picked;
        if (aero_picker.draw(picked)) {
            cfg.set_path("aerodynamics.aero_file", picked);
        }

        // The dimensional reference fields next to the path -- these
        // are intimately tied to whatever DATCOM run produced the
        // deck, so it makes sense to expose them next to the file
        // picker.
        edit_double(cfg, "Ref area [m^2]", "aerodynamics.refa",    "%.4f",
                    "Reference area (matches the area used in DATCOM run)");
        edit_double(cfg, "Ref diam [m]",   "aerodynamics.refd",    "%.4f",
                    "Reference diameter");
        edit_double(cfg, "xCG ref [m]",    "aerodynamics.xcg_ref", "%.4f",
                    "Reference CG location (DATCOM moments computed about this)");
    }

    if (ImGui::CollapsingHeader("Control / Autopilot")) {
        edit_slider(cfg, "Pitch cmd [g]",     "control.ancomx_program", -3.0, 3.0, "%.3f",
                    "Time-programmed pitch acceleration command");
        edit_slider(cfg, "Yaw cmd [g]",       "control.alcomx_program", -3.0, 3.0, "%.3f",
                    "Time-programmed yaw acceleration command");
        edit_double(cfg, "Pitch start [s]",   "control.t_pitch_start", "%.2f");
        edit_double(cfg, "Pitch end [s]",     "control.t_pitch_end",   "%.2f");
        edit_double(cfg, "Thrust loc [m]",    "control.thrust_loc",    "%.2f",
                    "Thrust attachment point from nose (used for moment arm)");

        ImGui::Spacing();
        ImGui::TextDisabled("Load relief (set gain=0 to disable)");
        edit_slider(cfg, "Accel relief gain",  "control.accel_relief_gain",        0.0, 2.0, "%.2f");
        edit_double(cfg, "Q threshold [Pa]",   "control.accel_relief_q_threshold", "%.0f");
        edit_double(cfg, "Q width [Pa]",       "control.accel_relief_q_width",     "%.0f");
        edit_double(cfg, "Washout tau [s]",    "control.accel_relief_tau",         "%.2f");
        edit_slider(cfg, "Autopilot detune",   "control.load_relief_tau_factor",   1.0, 10.0, "%.2f",
                    "Multiplies autopilot tau during load relief engagement");
    }

    if (ImGui::CollapsingHeader("Termination")) {
        edit_double(cfg, "End time [s]",     "sim.t_end",          "%.1f");
        edit_double(cfg, "Step size [s]",    "sim.dt",             "%.4f");
        edit_double(cfg, "Max altitude [m]", "intercept.max_alt",  "%.0f");
    }

    ImGui::PopStyleColor(3);  // header colors

    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    // Big run button.  Greyed out and styled neutral while a mission
    // is in flight on the worker thread; the GUI's own dispatch logic
    // also ignores a click in that state, so this is purely a UX hint.
    const bool busy = rocket6dof::progress::is_running.load();
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
    const char* label = busy ? "Running..." : "Run mission";
    if (ImGui::Button(label, ImVec2(200, 36))) {
        run_clicked = true;
    }
    if (busy) ImGui::EndDisabled();
    ImGui::PopStyleColor(3);

    return run_clicked;
}

// Advanced tab: full tree walk
//
// Recursively renders the JSON tree as nested collapsible sections.
// Every numeric leaf becomes an InputDouble; every string is shown
// read-only; every bool gets a checkbox.  Booleans and strings can't
// be edited via set_path -- the GUI shows them but routes the user
// back to the file Save/Load for those.

static void render_subtree(Value& cfg, const std::string& base_path,
                           const std::string& display_name,
                           int depth) {
    // Pull the subtree for read access
    // We re-walk from cfg root each time to find leaves (simple,
    // correct, and the configs are small so cost doesn't matter).
    const Value* node = &cfg;
    if (!base_path.empty()) {
        std::string seg;
        size_t i = 0;
        size_t n = base_path.size();
        while (i < n) {
            seg.clear();
            while (i < n && base_path[i] != '.' && base_path[i] != '[') {
                seg.push_back(base_path[i++]);
            }
            if (!seg.empty()) {
                if (!node->isObject()) return;
                node = &(*node)[seg];
            }
            while (i < n && base_path[i] == '[') {
                ++i;
                size_t idx = 0;
                while (i < n && base_path[i] >= '0' && base_path[i] <= '9') {
                    idx = idx * 10 + static_cast<size_t>(base_path[i++] - '0');
                }
                if (i < n && base_path[i] == ']') ++i;
                if (!node->isArray() || idx >= node->size()) return;
                node = &(*node)[idx];
            }
            if (i < n && base_path[i] == '.') ++i;
        }
    }

    if (node->isObject()) {
        bool open = (depth == 0 ? ImGui::CollapsingHeader(display_name.c_str())
                                : ImGui::TreeNode(display_name.c_str()));
        if (open) {
            ImGui::PushID(display_name.c_str());
            for (const auto& kv : node->obj()) {
                const std::string& key = kv.first;
                std::string child_path = base_path.empty()
                    ? key : (base_path + "." + key);
                render_subtree(cfg, child_path, key, depth + 1);
            }
            ImGui::PopID();
            if (depth != 0) ImGui::TreePop();
        }
    } else if (node->isArray()) {
        if (depth == 0 ? ImGui::CollapsingHeader(display_name.c_str())
                       : ImGui::TreeNode(display_name.c_str())) {
            ImGui::PushID(display_name.c_str());
            for (size_t i = 0; i < node->size(); ++i) {
                char idx_name[32];
                std::snprintf(idx_name, sizeof(idx_name), "[%zu]", i);
                std::string child_path = base_path + idx_name;
                render_subtree(cfg, child_path, idx_name, depth + 1);
            }
            ImGui::PopID();
            if (depth != 0) ImGui::TreePop();
        }
    } else if (node->isNumber()) {
        double v = node->asNumber();
        ImGui::PushID(base_path.c_str());
        if (ImGui::InputDouble(display_name.c_str(), &v, 0.0, 0.0, "%.6g")) {
            cfg.set_path(base_path, v);
        }
        ImGui::PopID();
    } else if (node->isString()) {
        // Editable text input bound to the JSON path.  Common cases
        // here: aerodynamics.aero_file (path to aero deck),
        // aerodynamics.tag_* (named tables inside the deck), and
        // per-input distribution names in monte_carlo configs.
        //
        // 256 chars is enough for any reasonable filesystem path.
        // The Value's stored string is the source of truth; we copy
        // into a temporary buffer for ImGui to mutate, then write
        // back via set_path on change.
        char buf[256];
        std::snprintf(buf, sizeof(buf), "%s", node->asString().c_str());
        ImGui::PushID(base_path.c_str());
        if (ImGui::InputText(display_name.c_str(), buf, sizeof(buf))) {
            cfg.set_path(base_path, std::string(buf));
        }
        ImGui::PopID();
    } else if (node->isBool()) {
        // Bool leaves are displayed read-only.  The configs in this
        // project don't expose bool leaves the user needs to toggle
        // from the GUI, so a Checkbox + set_path<bool> overload is
        // not provided.
        ImGui::Text("%s = %s", display_name.c_str(),
                    node->asBool() ? "true" : "false");
    } else if (node->isNull()) {
        ImGui::TextDisabled("%s = null", display_name.c_str());
    }
}

bool render_advanced_tab(Value& cfg) {
    bool run_clicked = false;

    ImGui::TextDisabled(
        "Full config tree.  Numeric and string leaves are editable; "
        "booleans are read-only here -- edit the JSON directly and "
        "reload.");
    ImGui::Spacing();

    // Push a brighter header color so the top-level sections stand
    // out against the dark window background.  ImGui's default header
    // color is very close to the panel color, which makes a list of
    // collapsed sections look like an empty panel.
    ImGui::PushStyleColor(ImGuiCol_Header,        IM_COL32( 70,  90, 130, 255));
    ImGui::PushStyleColor(ImGuiCol_HeaderHovered, IM_COL32( 90, 120, 170, 255));
    ImGui::PushStyleColor(ImGuiCol_HeaderActive,  IM_COL32(110, 140, 200, 255));

    // Walk top-level keys
    if (cfg.isObject()) {
        for (const auto& kv : cfg.obj()) {
            render_subtree(cfg, kv.first, kv.first, 0);
        }
    }

    ImGui::PopStyleColor(3);

    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    // Run button mirrors the Basic tab's: greyed out while a mission
    // is in flight on the worker thread.
    const bool busy = rocket6dof::progress::is_running.load();
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
    const char* label = busy ? "Running..." : "Run mission";
    if (ImGui::Button(label, ImVec2(200, 36))) {
        run_clicked = true;
    }
    if (busy) ImGui::EndDisabled();
    ImGui::PopStyleColor(3);

    return run_clicked;
}

// CSV reader for the post-run mission_log.csv
//
// Minimal parser, comma-separated, header row, all numeric.  Failures
// (file missing, non-numeric data) return false and the caller leaves
// the plots empty.

bool load_csv_columns(const std::string& path,
                      std::vector<std::string>& names,
                      std::vector<std::vector<double>>& cols) {
    std::ifstream f(path);
    if (!f) return false;
    names.clear();
    cols.clear();

    std::string line;
    if (!std::getline(f, line)) return false;
    // Header
    {
        std::stringstream ss(line);
        std::string field;
        while (std::getline(ss, field, ',')) {
            // Strip whitespace
            while (!field.empty() && (field.front() == ' ' || field.front() == '\t')) field.erase(field.begin());
            while (!field.empty() && (field.back()  == ' ' || field.back()  == '\t' || field.back() == '\r')) field.pop_back();
            names.push_back(field);
            cols.emplace_back();
        }
    }

    while (std::getline(f, line)) {
        std::stringstream ss(line);
        std::string field;
        size_t col = 0;
        while (std::getline(ss, field, ',')) {
            if (col >= cols.size()) break;
            try {
                cols[col].push_back(std::stod(field));
            } catch (...) {
                cols[col].push_back(std::nan(""));
            }
            ++col;
        }
        // Pad short rows
        while (col < cols.size()) {
            cols[col].push_back(std::nan(""));
            ++col;
        }
    }
    return !names.empty() && !cols[0].empty();
}

// Plot panel: altitude, velocity, q-alpha vs time

namespace {

// Find a named column.  Returns nullptr if missing.
const std::vector<double>* find_col(
    const std::vector<std::string>& names,
    const std::vector<std::vector<double>>& cols,
    const char* name) {
    for (size_t i = 0; i < names.size(); ++i) {
        if (names[i] == name) return &cols[i];
    }
    return nullptr;
}

}  // anon

PlotPanelActions render_plot_panel(
    const std::vector<std::string>& names,
    const std::vector<std::vector<double>>& cols,
    const std::vector<std::string>& cmp_names = {},
    const std::vector<std::vector<double>>& cmp_cols = {},
    const std::string& cmp_label = std::string())
{
    PlotPanelActions actions;

    if (names.empty() || cols.empty()) {
        ImGui::TextDisabled(
            "No mission log loaded.  Run a mission from the Basic or "
            "Advanced tab; results will appear here.");
        return actions;
    }

    const std::vector<double>* t = find_col(names, cols, "t");
    if (!t) {
        ImGui::TextDisabled("mission_log.csv has no 't' column?");
        return actions;
    }

    const std::vector<double>* alt    = find_col(names, cols, "alt");
    const std::vector<double>* dvbi   = find_col(names, cols, "dvbi");
    const std::vector<double>* q_dyn  = find_col(names, cols, "q_dyn");
    const std::vector<double>* alpha  = find_col(names, cols, "alpha");
    const std::vector<double>* mach   = find_col(names, cols, "mach");

    // Comparison-series shortcuts.  Only used when caller provided
    // a second CSV.
    const bool have_cmp = !cmp_names.empty() && !cmp_cols.empty();
    const std::vector<double>* tc      = have_cmp ? find_col(cmp_names, cmp_cols, "t")      : nullptr;
    const std::vector<double>* alt_c   = have_cmp ? find_col(cmp_names, cmp_cols, "alt")    : nullptr;
    const std::vector<double>* dvbi_c  = have_cmp ? find_col(cmp_names, cmp_cols, "dvbi")   : nullptr;
    const std::vector<double>* q_dyn_c = have_cmp ? find_col(cmp_names, cmp_cols, "q_dyn")  : nullptr;
    const std::vector<double>* alpha_c = have_cmp ? find_col(cmp_names, cmp_cols, "alpha")  : nullptr;
    const std::vector<double>* mach_c  = have_cmp ? find_col(cmp_names, cmp_cols, "mach")   : nullptr;

    // ---- Header row: status + buttons ----
    ImGui::Text("Loaded %zu rows from mission_log.csv (%zu columns)",
                t->size(), names.size());
    if (have_cmp && tc) {
        ImGui::SameLine();
        ImGui::TextDisabled("| comparison: %zu rows (%s)",
                            tc->size(),
                            cmp_label.empty() ? "compare" : cmp_label.c_str());
    }

    // Right-aligned buttons.  Layout right-to-left: Export, Clear cmp,
    // Load cmp.  Skip Clear if no comparison is loaded.
    //
    // Place all the buttons on the same line as the status text.
    // GetContentRegionAvail().x AFTER SameLine() is the remaining
    // space from the current cursor to the right edge, which is
    // exactly what we need to compute the right-edge offset.
    float export_w     = 150.0f;
    float load_cmp_w   = 170.0f;
    float clear_cmp_w  = 140.0f;
    float spacing      = 6.0f;
    float total_w = export_w + spacing + load_cmp_w;
    if (have_cmp) total_w += spacing + clear_cmp_w;

    ImGui::SameLine();
    float avail = ImGui::GetContentRegionAvail().x;
    if (avail > total_w) {
        ImGui::SetCursorPosX(ImGui::GetCursorPosX() + (avail - total_w));
    }
    if (ImGui::Button("Load comparison...", ImVec2(load_cmp_w, 0))) {
        actions.load_compare_clicked = true;
    }
    if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip(
            "Load a second mission_log CSV and overlay its curves on "
            "the current plots.  Useful for comparing two trajectories "
            "or configurations side by side.");
    }
    if (have_cmp) {
        ImGui::SameLine();
        if (ImGui::Button("Clear comparison", ImVec2(clear_cmp_w, 0))) {
            actions.clear_compare_clicked = true;
        }
    }
    ImGui::SameLine();
    if (ImGui::Button("Export as PNG", ImVec2(export_w, 0))) {
        actions.export_clicked = true;
    }
    if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip(
            "Save the current Plots view as plots_export.png in the "
            "working directory.  Captures the full window contents.");
    }
    ImGui::Spacing();

    // ---- Plots ----
    // Each plot draws the primary series, plus the comparison series
    // (against its own t vector) if available.  The labels in the
    // legend make it clear which line is which.
    auto plot_one = [&](const char* title, const char* yaxis,
                        const char* pri_label,
                        const std::vector<double>* pri_series,
                        const std::vector<double>* cmp_series)
    {
        std::vector<std::pair<const char*, const std::vector<double>*>> series;
        std::vector<std::pair<const char*,
                              std::pair<const std::vector<double>*,
                                        const std::vector<double>*>>> cmp_extra;
        series.push_back({pri_label, pri_series});
        std::string cmp_label_full;
        if (cmp_series && tc) {
            // Build a label like "alt (cmp)" or "alt (heavy)" using
            // the caller-supplied cmp_label.
            cmp_label_full = std::string(pri_label) + " (" +
                (cmp_label.empty() ? "cmp" : cmp_label) + ")";
        }
        // Custom draw because we need TWO different time vectors.
        if (ImPlot::BeginPlot(title, ImVec2(-1, 220.0f))) {
            ImPlot::SetupAxes("t [s]", yaxis);
            if (pri_series && !pri_series->empty()) {
                ImPlot::PlotLine(pri_label, t->data(), pri_series->data(),
                                 static_cast<int>(std::min(t->size(),
                                                           pri_series->size())));
            }
            if (cmp_series && !cmp_series->empty() && tc) {
                ImPlot::PlotLine(cmp_label_full.c_str(), tc->data(), cmp_series->data(),
                                 static_cast<int>(std::min(tc->size(),
                                                           cmp_series->size())));
            }
            ImPlot::EndPlot();
        }
    };

    plot_one("Altitude vs time",            "alt [m]",     "alt",   alt,    alt_c);
    plot_one("Inertial velocity vs time",   "dvbi [m/s]",  "dvbi",  dvbi,   dvbi_c);
    plot_one("Dynamic pressure vs time",    "q_dyn [Pa]",  "q_dyn", q_dyn,  q_dyn_c);
    plot_one("Angle of attack vs time",     "alpha [deg]", "alpha", alpha,  alpha_c);
    plot_one("Mach vs time",                "Mach",        "mach",  mach,   mach_c);

    return actions;
}

}  // namespace gui
