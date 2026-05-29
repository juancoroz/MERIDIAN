.SUFFIXES:
# =====================================================================
# Makefile for examples/rocket_6dof
#
#     make            # build both env_sweep and dyn_test
#     make run        # run both
#     make run-env    # run only env_sweep
#     make run-dyn    # run only dyn_test
#     make clean
# =====================================================================

CXX      ?= g++
CXXFLAGS ?= -std=c++14 -Wall -Wextra -O2
# Extra flags passed to every compile.  Use this to set preprocessor
# defines without overriding CXXFLAGS, e.g.
#   make CPPFLAGS=-DROCKET6DOF_MAX_STAGES=8
CPPFLAGS ?=

OSK_DIR  := ../osk

KERNEL_SRCS := \
    $(OSK_DIR)/state.cpp  \
    $(OSK_DIR)/block.cpp  \
    $(OSK_DIR)/sim.cpp    \
    $(OSK_DIR)/filer.cpp  \
    $(OSK_DIR)/table.cpp  \
    $(OSK_DIR)/vec.cpp    \
    $(OSK_DIR)/mat.cpp    \
    $(OSK_DIR)/quat.cpp

KERNEL_OBJS := $(notdir $(KERNEL_SRCS:.cpp=.o))
APP_SRCS    := environment.cpp newton.cpp kinematics.cpp euler.cpp propulsion.cpp aerodynamics.cpp control.cpp tvc.cpp guidance.cpp ins.cpp gps.cpp startrack.cpp rcs.cpp intercept.cpp forces.cpp peak_tracker.cpp cad_utils.cpp json.cpp mission_config.cpp aero_deck.cpp distributions.cpp sim_runner.cpp parallel_runner.cpp mc_runner.cpp sobol_runner.cpp mission_runner.cpp mission_progress.cpp
APP_OBJS    := $(APP_SRCS:.cpp=.o)

KERNEL_HDRS := $(wildcard $(OSK_DIR)/*.h)
APP_HDRS    := $(wildcard *.h)

ENV_BIN := env_sweep
DYN_BIN := dyn_test
CAD_BIN := cad_test
PROP_BIN := prop_test
AERO_BIN := aero_test
ICPT_BIN := intercept_test
CTRL_BIN := control_test
GUID_BIN := guidance_test
INS_BIN  := ins_test
GPS_BIN  := gps_test
STAR_BIN := startrack_test
RCS_BIN  := rcs_test
JSON_BIN := json_test
ADCK_BIN := aero_deck_test
DIST_BIN := distributions_test
INTEG_BIN := integration_test
MISSION_BIN := mission
MC_BIN   := mission_mc
SOBOL_BIN := mission_sobol
WORKER_BIN := mission_worker
VIEWER_BIN := build_mc_viewer
SVIEWER_BIN := build_viewer
MERGE_BIN := sobol_merge
GUI_BIN   := mission_gui

# ---- GUI build (Dear ImGui + ImPlot, vendored under third_party/) ----
# The GUI build is NOT part of the default `make all` target because it
# requires GLFW + OpenGL development headers, which may not be present
# on headless CI boxes.  Build with `make mission_gui` explicitly.
#
# Vendored sources (committed under third_party/):
IMGUI_DIR   := third_party/imgui
IMPLOT_DIR  := third_party/implot

IMGUI_SRCS := \
    $(IMGUI_DIR)/imgui.cpp \
    $(IMGUI_DIR)/imgui_draw.cpp \
    $(IMGUI_DIR)/imgui_tables.cpp \
    $(IMGUI_DIR)/imgui_widgets.cpp \
    $(IMGUI_DIR)/backends/imgui_impl_glfw.cpp \
    $(IMGUI_DIR)/backends/imgui_impl_opengl3.cpp \
    $(IMPLOT_DIR)/implot.cpp \
    $(IMPLOT_DIR)/implot_items.cpp

IMGUI_OBJS := $(notdir $(IMGUI_SRCS:.cpp=.o))
IMGUI_INC  := -I$(IMGUI_DIR) -I$(IMGUI_DIR)/backends -I$(IMPLOT_DIR)
GUI_LIBS   := -lglfw -lGL -lpthread -ldl

.PHONY: all run run-env run-dyn run-cad run-prop run-aero run-icpt run-ctrl run-guid run-ins run-gps run-star run-rcs run-json run-adck run-dist run-mission run-mc run-sobol run-viewer run-mc-viewer clean

all: $(ENV_BIN) $(DYN_BIN) $(CAD_BIN) $(PROP_BIN) $(AERO_BIN) $(ICPT_BIN) $(CTRL_BIN) $(GUID_BIN) $(INS_BIN) $(GPS_BIN) $(STAR_BIN) $(RCS_BIN) $(JSON_BIN) $(ADCK_BIN) $(DIST_BIN) $(INTEG_BIN) $(MISSION_BIN) $(MC_BIN) $(SOBOL_BIN) $(WORKER_BIN) $(VIEWER_BIN) $(SVIEWER_BIN) $(MERGE_BIN)

$(ENV_BIN): env_sweep.o $(APP_OBJS) $(KERNEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(DYN_BIN): dyn_test.o $(APP_OBJS) $(KERNEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(CAD_BIN): cad_test.o $(APP_OBJS) $(KERNEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(PROP_BIN): prop_test.o $(APP_OBJS) $(KERNEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(AERO_BIN): aero_test.o $(APP_OBJS) $(KERNEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(ICPT_BIN): intercept_test.o $(APP_OBJS) $(KERNEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(CTRL_BIN): control_test.o $(APP_OBJS) $(KERNEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(GUID_BIN): guidance_test.o $(APP_OBJS) $(KERNEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(INS_BIN): ins_test.o $(APP_OBJS) $(KERNEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(GPS_BIN): gps_test.o $(APP_OBJS) $(KERNEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(STAR_BIN): startrack_test.o $(APP_OBJS) $(KERNEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(RCS_BIN): rcs_test.o $(APP_OBJS) $(KERNEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(JSON_BIN): json_test.o $(APP_OBJS) $(KERNEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(ADCK_BIN): aero_deck_test.o $(APP_OBJS) $(KERNEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(DIST_BIN): distributions_test.o $(APP_OBJS) $(KERNEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(INTEG_BIN): integration_test.o $(APP_OBJS) $(KERNEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(MISSION_BIN): mission.o $(APP_OBJS) $(KERNEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(MC_BIN): mission_mc.o $(APP_OBJS) $(KERNEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(SOBOL_BIN): mission_sobol.o $(APP_OBJS) $(KERNEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(WORKER_BIN): mission_worker.o $(APP_OBJS) $(KERNEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Standalone viewer-builder: no kernel or app deps, just embeds CSVs
# into the HTML template.  Kept separate so it builds quickly.
$(VIEWER_BIN): build_mc_viewer.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

# Single-mission viewer builder (mission_log.csv -> mission_viewer.html)
$(SVIEWER_BIN): build_viewer.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

# Sobol chunk-merge tool: combines independent chunks (each at N=N_chunk
# with different RNG seeds) into a single CSV with averaged S1/ST point
# estimates and chunk-to-chunk-variance-based CIs.  Statistically valid
# when chunks are independent, which the chunk launcher guarantees by
# using distinct seeds.  Standalone -- no OSK deps -- so builds fast.
$(MERGE_BIN): sobol_merge.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

# ---- GUI binary: mission_gui ----
# Dear ImGui + ImPlot frontend.  Same simulation pipeline as the CLI,
# but with parameter-editor widgets and post-run plots.
#
# Build:  make mission_gui
# Run:    ./mission_gui [config.json]
#
# Not part of `make all` -- requires GLFW + OpenGL devel headers
# (Ubuntu: apt install libglfw3-dev libgl1-mesa-dev).

# ImGui/ImPlot pattern rules.  Their sources live under third_party/
# but we compile them into top-level .o files for simpler linking.
imgui.o:               $(IMGUI_DIR)/imgui.cpp
	$(CXX) $(CXXFLAGS) $(IMGUI_INC) -c $< -o $@
imgui_draw.o:          $(IMGUI_DIR)/imgui_draw.cpp
	$(CXX) $(CXXFLAGS) $(IMGUI_INC) -c $< -o $@
imgui_tables.o:        $(IMGUI_DIR)/imgui_tables.cpp
	$(CXX) $(CXXFLAGS) $(IMGUI_INC) -c $< -o $@
imgui_widgets.o:       $(IMGUI_DIR)/imgui_widgets.cpp
	$(CXX) $(CXXFLAGS) $(IMGUI_INC) -c $< -o $@
imgui_impl_glfw.o:     $(IMGUI_DIR)/backends/imgui_impl_glfw.cpp
	$(CXX) $(CXXFLAGS) $(IMGUI_INC) -c $< -o $@
imgui_impl_opengl3.o:  $(IMGUI_DIR)/backends/imgui_impl_opengl3.cpp
	$(CXX) $(CXXFLAGS) $(IMGUI_INC) -c $< -o $@
implot.o:              $(IMPLOT_DIR)/implot.cpp
	$(CXX) $(CXXFLAGS) $(IMGUI_INC) -c $< -o $@
implot_items.o:        $(IMPLOT_DIR)/implot_items.cpp
	$(CXX) $(CXXFLAGS) $(IMGUI_INC) -c $< -o $@

# GUI-specific app sources need IMGUI_INC on their compile lines.
# Override the pattern rule for these explicitly.
mission_gui.o: mission_gui.cpp $(APP_HDRS) $(KERNEL_HDRS) gui_plot_panel.h gui_queue.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -I$(OSK_DIR) $(IMGUI_INC) -c $< -o $@
gui_config.o:  gui_config.cpp  $(APP_HDRS) $(KERNEL_HDRS) gui_plot_panel.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -I$(OSK_DIR) $(IMGUI_INC) -c $< -o $@
gui_file_picker.o: gui_file_picker.cpp gui_file_picker.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(IMGUI_INC) -c $< -o $@
gui_subprocess.o: gui_subprocess.cpp gui_subprocess.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@
gui_sweeps.o: gui_sweeps.cpp gui_sweeps.h gui_subprocess.h gui_file_picker.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(IMGUI_INC) -c $< -o $@
gui_screenshot.o: gui_screenshot.cpp gui_screenshot.h third_party/stb/stb_image_write.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -Ithird_party/stb -c $< -o $@
gui_queue.o: gui_queue.cpp gui_queue.h gui_file_picker.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(IMGUI_INC) -c $< -o $@

$(GUI_BIN): mission_gui.o gui_config.o gui_file_picker.o gui_subprocess.o gui_sweeps.o gui_screenshot.o gui_queue.o $(IMGUI_OBJS) $(APP_OBJS) $(KERNEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ $(GUI_LIBS) -o $@

%.o: %.cpp $(APP_HDRS) $(KERNEL_HDRS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -I$(OSK_DIR) -c $< -o $@

$(KERNEL_OBJS): %.o: $(OSK_DIR)/%.cpp $(KERNEL_HDRS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -I$(OSK_DIR) -c $< -o $@

run: run-env run-dyn run-cad run-prop run-aero run-icpt run-ctrl run-guid run-ins run-gps run-star run-rcs

run-env: $(ENV_BIN)
	./$(ENV_BIN)

run-dyn: $(DYN_BIN)
	./$(DYN_BIN)

run-cad: $(CAD_BIN)
	./$(CAD_BIN)

run-prop: $(PROP_BIN)
	./$(PROP_BIN)

run-aero: $(AERO_BIN)
	./$(AERO_BIN)

run-icpt: $(ICPT_BIN)
	./$(ICPT_BIN)

run-ctrl: $(CTRL_BIN)
	./$(CTRL_BIN)

run-guid: $(GUID_BIN)
	./$(GUID_BIN)

run-ins: $(INS_BIN)
	./$(INS_BIN)

run-gps: $(GPS_BIN)
	./$(GPS_BIN)

run-star: $(STAR_BIN)
	./$(STAR_BIN)

run-rcs: $(RCS_BIN)
	./$(RCS_BIN)

run-json: $(JSON_BIN)
	./$(JSON_BIN)

run-adck: $(ADCK_BIN)
	./$(ADCK_BIN)

run-dist: $(DIST_BIN)
	./$(DIST_BIN)

run-mission: $(MISSION_BIN)
	./$(MISSION_BIN)

run-mc: $(MC_BIN)
	./$(MC_BIN)

run-sobol: $(SOBOL_BIN)
	./$(SOBOL_BIN)

# ---- MC viewer ----
# Regenerate mc_viewer.html from the latest MC CSVs and (optionally)
# Sobol CSV.  Run mission_mc and mission_sobol first to produce the
# CSVs that this rule depends on.
mc_viewer.html: monte_carlo_outputs.csv monte_carlo_inputs.csv mc_viewer_template.html $(VIEWER_BIN)
	./$(VIEWER_BIN)

run-mc-viewer: mc_viewer.html
	@echo "Open mc_viewer.html in a browser to view the MC dispersion plots"

mission_log.csv: $(MISSION_BIN)
	./$(MISSION_BIN) >/dev/null

mission_viewer.html: mission_log.csv mission_viewer_template.html $(SVIEWER_BIN)
	./$(SVIEWER_BIN)

run-viewer: mission_viewer.html
	@echo "Open mission_viewer.html in a browser to view the trajectory plots"

clean:
	rm -f $(ENV_BIN) $(DYN_BIN) $(CAD_BIN) $(PROP_BIN) $(AERO_BIN) $(ICPT_BIN) $(CTRL_BIN) $(GUID_BIN) $(INS_BIN) $(GPS_BIN) $(STAR_BIN) $(RCS_BIN) $(JSON_BIN) $(ADCK_BIN) $(DIST_BIN) $(MISSION_BIN) $(MC_BIN) $(SOBOL_BIN) $(WORKER_BIN) $(VIEWER_BIN) $(SVIEWER_BIN) $(MERGE_BIN) $(GUI_BIN) *.o mission_log.csv mission_viewer.html monte_carlo_*.csv monte_carlo_summary.txt mc_viewer.html sobol_indices.csv sobol_summary.txt sobol_indices_reduced.csv sobol_summary_reduced.txt
