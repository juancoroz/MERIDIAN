//  mission_gui.cpp -- Dear ImGui frontend for the rocket_6dof simulator
//
//  Two-tab parameter editor + run button + post-run plots.
//
//  Architecture:
//    1. Load default config (mission.json or last-used) into an
//       in-memory copy of the same JSON structure mission_config.cpp
//       uses.
//    2. ImGui renders an editor over that JSON tree.  Two top-level
//       tabs: "Basic" (curated subset of ~30 params) and "Advanced"
//       (the full structured tree, organized by config section).
//    3. On "Run mission", we serialize the live config to a temp file
//       and invoke run_mission(temp_path).  Re-uses all existing
//       JSON loading and validation -- the GUI can't drift from the
//       CLI's view of what's valid.
//    4. After run completes, the mission_log.csv is read back and
//       altitude/velocity/q-alpha are plotted with ImPlot.
//
//  Build: see Makefile target mission_gui.
//
//  Run modes:
//    ./mission             -- CLI: loads mission.json
//    ./mission config.json -- CLI: loads named config
//    ./mission_gui         -- GUI: opens with mission.json loaded
//    ./mission_gui foo.json-- GUI: opens with foo.json loaded
//
//  Headless safety: this binary requires a display; if GLFW fails to
//  initialize (e.g. server with no X), it prints an error and exits.
//  The CLI ./mission binary remains usable on headless boxes for CI
//  and Sobol chunk runs.

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "implot.h"

#include <GLFW/glfw3.h>

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <thread>
#include <atomic>
#include <chrono>
#include <memory>

// JSON loading is shared with the CLI.
#include "json.h"

// File picker dialog (pure ImGui, no system dep).
#include "gui_file_picker.h"

// Progress atomics published by mission_runner.cpp during a run.
#include "mission_progress.h"

// Sweeps tab (MC/Sobol/viewer subprocess runners).
#include "gui_sweeps.h"

// Framebuffer screenshot helper.
#include "gui_screenshot.h"

// Plot panel action-flags struct.
#include "gui_plot_panel.h"

// Queue tab (sequential multi-config runner).
#include "gui_queue.h"

// Forward declarations
// Implemented in gui_config.cpp.  Both treat the loaded config as a
// mutable Value tree; edits flow back into the same tree.
namespace gui {
    bool load_config_file(const std::string& path, rocket6dof::json::Value& out);
    bool save_config_file(const std::string& path, const rocket6dof::json::Value& cfg);

    // Render the Basic-tab editor on the current cfg.  Returns true
    // if the user clicked "Run mission" this frame.
    bool render_basic_tab(rocket6dof::json::Value& cfg);

    // Render the Advanced-tab editor.  Same return convention.
    bool render_advanced_tab(rocket6dof::json::Value& cfg);

    // Read the post-run mission_log.csv into vectors keyed by column
    // name.  Returns false on failure.
    bool load_csv_columns(const std::string& path,
                          std::vector<std::string>& names,
                          std::vector<std::vector<double>>& cols);

    // Action flags returned by render_plot_panel each frame.  Defined
    // in gui_plot_panel.h.

    // Render the post-run plot panel given the loaded CSV columns and
    // (optionally) a second comparison CSV to overlay.  Returns a set
    // of action flags indicating what the user clicked this frame.
    PlotPanelActions render_plot_panel(
        const std::vector<std::string>& names,
        const std::vector<std::vector<double>>& cols,
        const std::vector<std::string>& cmp_names = {},
        const std::vector<std::vector<double>>& cmp_cols = {},
        const std::string& cmp_label = std::string());
}

// run_mission_ext is defined in mission_runner.cpp (extracted there
// for reuse).  The `_ext` suffix is just a symbol-name reminder that
// this is the external-linkage trampoline; the actual pipeline body
// is in the anonymous namespace inside mission_runner.cpp.
extern "C" int run_mission_ext(const char* config_path);

// GUI state
namespace {

struct AppState {
    rocket6dof::json::Value                       cfg;
    std::string                      cfg_path = "mission.json";
    std::string                      last_run_status;
    int                              last_run_rc = 0;

    // Post-run CSV columns
    std::vector<std::string>         csv_col_names;
    std::vector<std::vector<double>> csv_col_data;
    bool                             have_csv = false;

    // Optional second CSV loaded via "Load comparison..." on the Plots
    // tab.  When present, render_plot_panel overlays its curves on
    // the primary plots.  Empty by default.
    std::vector<std::string>         cmp_col_names;
    std::vector<std::vector<double>> cmp_col_data;
    std::string                      cmp_label;     // shown in legend
    gui::FilePicker                  cmp_picker;

    // Load/Save file picker dialogs.  Separate instances so they
    // remember their own last-visited directories independently.
    gui::FilePicker                  load_picker;
    gui::FilePicker                  save_picker;

    // Sweeps tab state (MC/Sobol/viewer subprocess runners).
    gui::SweepsTabState              sweeps;

    // Queue tab state (sequential multi-config runner).
    gui::QueueTabState               queue;

    // Worker thread for the Run button.  The mission runs on this
    // thread so the GUI stays responsive (rendering at ~60 FPS) while
    // the simulation integrates.  Progress published via the atomics
    // in <mission_progress.h>; the GUI polls them every frame.
    //
    // Owned by unique_ptr because std::thread isn't trivially default-
    // constructible/movable in a struct member without surrounding
    // boilerplate -- unique_ptr<thread> sidesteps that.
    std::unique_ptr<std::thread>     worker;

    // Set by the worker thread on completion; the GUI thread polls
    // this each frame.  Separate from progress::is_running so we can
    // tell "thread finished" from "haven't started yet".
    std::atomic<bool>                worker_done {false};
    int                              worker_rc = 0;
};

// Print a GLFW error to stderr.  GLFW invokes this on startup failures.
void glfw_error_cb(int code, const char* desc) {
    std::fprintf(stderr, "GLFW error %d: %s\n", code, desc);
}

}  // anon

// Main
int main(int argc, char** argv) {
    AppState app;

    // Default config is mission.json.  Allow override:
    //   ./mission_gui foo.json
    if (argc >= 2) {
        app.cfg_path = argv[1];
    }

    if (!gui::load_config_file(app.cfg_path, app.cfg)) {
        std::fprintf(stderr,
            "warn: failed to load %s; starting with empty config\n",
            app.cfg_path.c_str());
        app.cfg = rocket6dof::json::Value::makeObject();
    }

    // If a mission_log.csv exists in the working directory, load it on
    // startup so the Plots tab is immediately useful.  This is purely
    // convenience -- a run never depends on it.
    app.have_csv = gui::load_csv_columns(
        "mission_log.csv", app.csv_col_names, app.csv_col_data);
    if (app.have_csv) {
        app.last_run_status = "Loaded existing mission_log.csv";
    }

    // ---- GLFW init ----
    glfwSetErrorCallback(glfw_error_cb);
    if (!glfwInit()) {
        std::fprintf(stderr,
            "error: glfwInit() failed.  Is a display available?\n"
            "       The CLI './mission %s' still works on headless boxes.\n",
            app.cfg_path.c_str());
        return 1;
    }

    // OpenGL 3.3 core profile -- widely supported across Linux/Mac/Win
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    GLFWwindow* window = glfwCreateWindow(
        1400, 900, "rocket_6dof mission GUI", nullptr, nullptr);
    if (!window) {
        std::fprintf(stderr,
            "error: glfwCreateWindow() failed.  No display?\n");
        glfwTerminate();
        return 1;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);  // vsync on

    // ---- ImGui init ----
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    ImGui::StyleColorsDark();

    if (!ImGui_ImplGlfw_InitForOpenGL(window, true)) {
        std::fprintf(stderr, "error: ImGui_ImplGlfw_InitForOpenGL failed\n");
        return 1;
    }
    if (!ImGui_ImplOpenGL3_Init("#version 330")) {
        std::fprintf(stderr, "error: ImGui_ImplOpenGL3_Init failed\n");
        return 1;
    }

    // ---- Main loop ----
    bool run_requested = false;
    bool plot_export_requested = false;
    bool queue_run_requested = false;
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // Full-window main panel
        ImGuiViewport* vp = ImGui::GetMainViewport();
        ImGui::SetNextWindowPos(vp->WorkPos);
        ImGui::SetNextWindowSize(vp->WorkSize);
        ImGui::Begin("rocket_6dof", nullptr,
            ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoBringToFrontOnFocus |
            ImGuiWindowFlags_NoCollapse);

        // ---- Header: config file + load/save buttons ----
        ImGui::Text("Config: %s", app.cfg_path.c_str());
        ImGui::SameLine();
        if (ImGui::Button("Load")) {
            // Open a file picker dialog rooted at the directory of the
            // currently-loaded config (helps if the user is working in
            // a directory of related scenarios).
            std::string default_name;
            size_t slash = app.cfg_path.find_last_of('/');
            if (slash != std::string::npos) {
                default_name = app.cfg_path.substr(slash + 1);
            } else {
                default_name = app.cfg_path;
            }
            app.load_picker.open("Load config", false,
                                 default_name, ".json");
        }
        ImGui::SameLine();
        if (ImGui::Button("Save")) {
            // Save picker pre-fills with current filename so users can
            // either overwrite (just confirm) or type a new name.
            std::string default_name;
            size_t slash = app.cfg_path.find_last_of('/');
            if (slash != std::string::npos) {
                default_name = app.cfg_path.substr(slash + 1);
            } else {
                default_name = app.cfg_path;
            }
            app.save_picker.open("Save config", true,
                                 default_name, ".json");
        }
        ImGui::SameLine();
        ImGui::TextDisabled("  |  %s", app.last_run_status.c_str());

        // ---- Progress bar (visible only while a mission is running) ----
        // Read both atomics in one shot; minor race is fine -- the bar
        // just updates a frame later.  total_t can be 0 right at the
        // edge of starting; guard against div-by-zero.
        if (rocket6dof::progress::is_running.load()) {
            double ct = rocket6dof::progress::current_t.load();
            double tt = rocket6dof::progress::total_t.load();
            float frac = (tt > 0.0) ? static_cast<float>(ct / tt) : 0.0f;
            if (frac < 0.0f) frac = 0.0f;
            if (frac > 1.0f) frac = 1.0f;
            char buf[64];
            std::snprintf(buf, sizeof(buf), "t = %.1f / %.1f s", ct, tt);
            ImGui::ProgressBar(frac, ImVec2(-1.0f, 0.0f), buf);
        }

        // ---- File picker dialogs (drawn every frame; cheap when closed) ----
        std::string picked_path;
        if (app.load_picker.draw(picked_path)) {
            if (gui::load_config_file(picked_path, app.cfg)) {
                app.cfg_path = picked_path;
                app.last_run_status = "Loaded " + picked_path;
            } else {
                app.last_run_status = "Load failed: " + picked_path;
            }
        }
        if (app.save_picker.draw(picked_path)) {
            // Auto-append .json if missing (the dialog filters but
            // doesn't enforce on typed names).
            if (picked_path.size() < 5 ||
                picked_path.substr(picked_path.size() - 5) != ".json") {
                picked_path += ".json";
            }
            if (gui::save_config_file(picked_path, app.cfg)) {
                app.cfg_path = picked_path;
                app.last_run_status = "Saved " + picked_path;
            } else {
                app.last_run_status = "Save failed: " + picked_path;
            }
        }

        ImGui::Separator();

        // ---- Tabs ----
        if (ImGui::BeginTabBar("MainTabs")) {
            if (ImGui::BeginTabItem("Basic")) {
                if (gui::render_basic_tab(app.cfg)) run_requested = true;
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Advanced")) {
                if (gui::render_advanced_tab(app.cfg)) run_requested = true;
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Plots")) {
                gui::PlotPanelActions plot_actions = gui::render_plot_panel(
                    app.csv_col_names, app.csv_col_data,
                    app.cmp_col_names, app.cmp_col_data,
                    app.cmp_label);
                if (plot_actions.export_clicked) {
                    plot_export_requested = true;
                }
                if (plot_actions.load_compare_clicked) {
                    // Filter accepts .csv only; default to ./mission_log.csv
                    app.cmp_picker.open("Load comparison CSV", false,
                                        "mission_log.csv", ".csv");
                }
                if (plot_actions.clear_compare_clicked) {
                    app.cmp_col_names.clear();
                    app.cmp_col_data.clear();
                    app.cmp_label.clear();
                    app.last_run_status = "Comparison cleared";
                }

                // Draw the comparison picker (cheap when closed).
                std::string cmp_path;
                if (app.cmp_picker.draw(cmp_path)) {
                    std::vector<std::string> nm;
                    std::vector<std::vector<double>> co;
                    if (gui::load_csv_columns(cmp_path, nm, co)) {
                        app.cmp_col_names = std::move(nm);
                        app.cmp_col_data  = std::move(co);
                        // Derive a short label from the filename
                        size_t sl = cmp_path.find_last_of('/');
                        app.cmp_label = (sl == std::string::npos)
                                        ? cmp_path
                                        : cmp_path.substr(sl + 1);
                        app.last_run_status = "Comparison: " + app.cmp_label;
                    } else {
                        app.last_run_status = "Failed to load comparison: " + cmp_path;
                    }
                }
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Sweeps")) {
                gui::render_sweeps_tab(app.sweeps);
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Queue")) {
                if (gui::render_queue_tab(app.queue)) {
                    queue_run_requested = true;
                }
                ImGui::EndTabItem();
            }
            ImGui::EndTabBar();
        }

        ImGui::End();

        // ---- Handle a run request ----
        // (Outside the window so we can synchronously block.  For an
        //  MVP this is fine -- launcher takes <1s.  Future: worker
        //  thread + progress bar.)
        // ---- Run dispatch ----
        // The run_requested flag is set from the Basic/Advanced tabs
        // when the user clicks "Run mission".  We refuse to start a
        // new run if one is already in flight (the buttons in those
        // tabs are styled disabled when busy, but we double-check
        // here to handle a stale request from a previous frame).
        const bool busy = rocket6dof::progress::is_running.load();
        if (run_requested && !busy) {
            run_requested = false;
            // Write the in-memory config to a temp file
            const std::string tmp = "/tmp/mission_gui_run.json";
            if (!gui::save_config_file(tmp, app.cfg)) {
                app.last_run_status = "Failed to write temp config";
            } else {
                // If a previous worker exists, join it before spawning
                // a new one.  worker_done would be true by now since
                // we only get here when !busy.
                if (app.worker && app.worker->joinable()) {
                    app.worker->join();
                }
                app.worker_done.store(false);
                app.last_run_status = "Running...";
                // Launch the mission on a worker thread.  Capture the
                // config path by value (it's a std::string, lifetime
                // outlives the thread by construction).
                std::string tmp_path = tmp;
                app.worker.reset(new std::thread([tmp_path, &app]() {
                    int rc = run_mission_ext(tmp_path.c_str());
                    app.worker_rc = rc;
                    app.worker_done.store(true);
                }));
            }
        }

        // ---- Poll worker completion ----
        // Single-shot on the frame the worker_done flag flips.
        if (app.worker_done.load() && app.worker) {
            // Don't process again until the next run.
            app.worker_done.store(false);
            // Joining is safe and quick here: the worker is finished
            // (it set worker_done last) so join() returns immediately.
            if (app.worker->joinable()) app.worker->join();
            app.worker.reset();

            int rc = app.worker_rc;
            app.last_run_rc = rc;
            if (rc == 0) {
                app.last_run_status = "Run complete (rc=0)";
                app.have_csv = gui::load_csv_columns(
                    "mission_log.csv", app.csv_col_names, app.csv_col_data);
            } else {
                app.last_run_status =
                    "Run failed (rc=" + std::to_string(rc) + ")";
            }
        }

        // ---- Queue dispatch ----
        // Launches a worker thread that iterates all queued items and
        // runs them sequentially.  Each item is treated like a Basic-
        // tab Run -- read JSON from disk, call run_mission_ext, parse
        // mission_log.csv for final altitude.  Cancellation is checked
        // BETWEEN items only (run_mission_ext itself is synchronous).
        if (queue_run_requested && !busy) {
            queue_run_requested = false;
            if (app.queue.worker && app.queue.worker->joinable()) {
                app.queue.worker->join();
            }
            app.queue.worker_done.store(false);
            app.queue.cancel_requested.store(false);
            app.queue.summary_text.clear();

            app.queue.worker.reset(new std::thread([&app]() {
                int n_ok = 0, n_fail = 0, n_cancel = 0;
                for (auto& item : app.queue.items) {
                    if (item.status != gui::QueueItem::Queued) continue;
                    if (app.queue.cancel_requested.load()) {
                        item.status = gui::QueueItem::Cancelled;
                        ++n_cancel;
                        continue;
                    }
                    item.status = gui::QueueItem::Running;
                    auto t0 = std::chrono::steady_clock::now();
                    int rc = run_mission_ext(item.path.c_str());
                    auto t1 = std::chrono::steady_clock::now();
                    item.elapsed_sec =
                        std::chrono::duration<double>(t1 - t0).count();
                    item.rc = rc;
                    if (rc == 0) {
                        item.status = gui::QueueItem::Done;
                        // Parse the final altitude from mission_log.csv.
                        // We read the last row's "alt" column; this is
                        // a single line of text from the bottom of the
                        // file so cost is O(small).
                        std::vector<std::string> cn;
                        std::vector<std::vector<double>> cd;
                        if (gui::load_csv_columns("mission_log.csv", cn, cd)) {
                            // Find alt column
                            for (size_t k = 0; k < cn.size(); ++k) {
                                if (cn[k] == "alt" && !cd[k].empty()) {
                                    item.final_alt_km = cd[k].back() / 1000.0;
                                    break;
                                }
                            }
                        }
                        ++n_ok;
                    } else {
                        item.status = gui::QueueItem::Failed;
                        ++n_fail;
                    }
                }
                char buf[160];
                std::snprintf(buf, sizeof(buf),
                              "Queue done: %d ok, %d failed, %d cancelled",
                              n_ok, n_fail, n_cancel);
                app.queue.summary_text = buf;
                app.queue.worker_done.store(true);
            }));
        }

        // Reap finished queue worker (join cleanly).
        if (app.queue.worker_done.load() && app.queue.worker) {
            app.queue.worker_done.store(false);
            if (app.queue.worker->joinable()) app.queue.worker->join();
            app.queue.worker.reset();
            // Refresh the main Plots tab CSV view to reflect the LAST
            // queue item's output -- usually what the user wants to
            // glance at after the batch completes.
            app.have_csv = gui::load_csv_columns(
                "mission_log.csv", app.csv_col_names, app.csv_col_data);
        }

        // ---- Render ----
        ImGui::Render();
        int dw, dh;
        glfwGetFramebufferSize(window, &dw, &dh);
        glViewport(0, 0, dw, dh);
        glClearColor(0.10f, 0.12f, 0.14f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        // ---- Plot export, if requested ----
        // Read pixels from the back buffer AFTER ImGui rendered into
        // it but BEFORE SwapBuffers makes it the front buffer (after
        // which the back buffer contents become undefined).
        if (plot_export_requested) {
            plot_export_requested = false;
            const std::string out = "plots_export.png";
            if (gui::save_framebuffer_png(out, 0, 0, dw, dh)) {
                app.last_run_status = "Exported " + out;
            } else {
                app.last_run_status = "Export to " + out + " failed";
            }
        }

        glfwSwapBuffers(window);
    }

    // ---- Cleanup ----
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
