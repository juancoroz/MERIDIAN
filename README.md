# MERIDIAN

**Monte Carlo Enabled Rocket Integrated Dynamics & Insight ANalysis**

A 6-degree-of-freedom rocket flight simulator built on an OSK-style
kernel.
Includes a command-line driver, Monte Carlo / Sobol sweep runners, an
interactive HTML viewer for sweep results, and a Dear ImGui-based GUI
with five tabs: Basic, Advanced, Plots, Sweeps, Queue.

MERIDIAN is the project name; the application source under
`rocket_6dof/` is the rocket-simulation codebase, and `osk/`
is the underlying kernel. The directory and binary names use
`rocket_6dof`, `mission`, `mission_gui`, etc. as historical
identifiers; these are the same code as MERIDIAN — just named at the
file-system level for the example they implement.

## About the kernel

The `osk/` directory contains an **OSK-style** simulation kernel —
inspired by the design conventions of the original Object Simulation
Kernel (OSK) created by Ray Sells, but not derived from or affiliated
with that codebase. The original OSK is not publicly released; the
kernel shipped here is an independent implementation that follows the
same broad architectural ideas (state objects, blocks, sim loop,
filer-based output) so users familiar with OSK conventions will
recognize the structure. It is not a drop-in substitute for the
original.

If you have access to Ray Sells' actual OSK and prefer to use it
instead, the application code in `rocket_6dof/` was written
against this kernel's API surface and would need light porting to
match the original.

## About the simulation models

The physics and subsystem models in `rocket_6dof/` —
dynamics, aerodynamics, propulsion, control, guidance, INS, GPS,
star tracker, RCS, TVC, and the coordinate-frame / kinematics
machinery — are based on the modeling approach of:

> Peter H. Zipfel, *Modeling and Simulation of Aerospace Vehicle
> Dynamics*, Third Edition (AIAA Education Series).

and on the accompanying **CADAC++** simulation code that Zipfel
provides alongside the book. The state-variable conventions, the
6-DOF equations of motion, the modular subsystem decomposition, and
much of the frame-transformation bookkeeping follow CADAC++ practice
so that readers familiar with the book and its code will recognize
the structure.

This codebase is an independent reimplementation in the OSK-style
kernel; it is not Zipfel's CADAC++ itself, and is neither derived
from nor affiliated with the author or publisher. If you want the
original CADAC++ reference code, it is distributed with the book.
Anyone using MERIDIAN for serious work should treat Zipfel's book as
the authoritative reference for the underlying models.

## Bundle layout

```
meridian/
├── README.md                     This file (project overview)
├── Makefile                      Top-level driver for the OSK kernel examples
├── osk/                          OSK-style kernel sources
│   └── *.cpp, *.h
├── rocket_6dof/                  6DOF rocket simulator (the main application)
│   ├── Makefile                  Build rules
│   ├── *.cpp, *.h                Application sources
│   ├── mission*.json             Mission configs
│   ├── monte_carlo*.json         Sweep configs
│   ├── aero*.txt, *.asc          Aero deck data
│   ├── *_template.html           HTML viewer templates
│   ├── BUILD_LINUX.md            Linux build instructions
│   └── third_party/              EMPTY -- you must populate this
│       ├── imgui/                Dear ImGui (see below)
│       ├── implot/               ImPlot (see below)
│       └── stb/                  stb_image_write.h (see below)
└── examples/                     OSK kernel demonstration examples
    ├── hello.cpp                 Smallest possible OSK driver
    ├── shm.cpp                   Simple harmonic motion
    ├── ex_1/ … ex_5/             CMD §4.1 - §4.5 (single-stage and multi-stage)
    ├── ex_6_1/ … ex_6_4/         CMD §4.6 (asynchronous event scheduling)
    ├── ex_app2/                  CMD Appendix 2 (RK2 plug-in integrator)
    ├── ex_1_io/                  Filer-driven variant of ex_1
    ├── util_demo/                Filer + Table regression tests
    ├── vmq_demo/                 Vec/Mat/Quat regression tests
    └── sync_test/                update/rpt source-order preservation
```

The `rocket_6dof/` directory is the main deliverable — a full 6DOF
rocket simulator built on the OSK kernel. The `examples/` directory
contains the OSK manual's worked examples: small, self-contained
programs that demonstrate one kernel feature each, useful as a
reference when reading the OSK manual or when building a new
simulation on the kernel.

## Quick start

### Linux

```sh
cd rocket_6dof
# Install build deps (Ubuntu/Debian)
sudo apt install build-essential libglfw3-dev libgl1-mesa-dev

# See BUILD_LINUX.md for third_party setup, then:
make all              # builds CLI binaries: mission, mission_mc, mission_sobol, etc.
make mission_gui      # builds GUI (requires GLFW + OpenGL headers)

./mission mission_launcher.json   # runs the sample mission
./mission_gui                     # launches the GUI
```

To build the OSK kernel examples instead, run `make` from the top
of the bundle (the driver Makefile builds every example in
`examples/`):

```sh
cd meridian/      # the top-level directory
make              # builds examples/hello, examples/shm, examples/ex_1/main, ...
```

Expected baseline output (the "regression test" we use to verify a
working build):

```
Final altitude:   247089.7 m  (247.09 km)
Final dvbi:       9129.46 m/s
```

### Windows

MERIDIAN does not build natively on Windows. The parallel sweep
runner uses POSIX `fork()`/`execvp()`/`waitpid()` directly. Windows
users should build and run inside **WSL** (Windows Subsystem for
Linux) with Ubuntu 22.04 or 24.04 -- everything in this README's
Linux quick-start works there unchanged.

To open the generated HTML viewers (`mission_viewer.html`,
`mc_viewer.html`) in your Windows browser from a WSL shell:

```sh
explorer.exe mission_viewer.html      # uses Windows default browser
```

## Third-party dependencies you must download

The `third_party/` subdirectories ship empty. You must download three
freely-available libraries and place them as documented below.
Versions listed are the ones this codebase was developed and tested
against; newer minor releases will likely work but are unverified.

### 1. Dear ImGui (immediate-mode GUI library)

- Home:    https://github.com/ocornut/imgui
- Version tested: **v1.91.6**
- Direct download:
  https://github.com/ocornut/imgui/archive/refs/tags/v1.91.6.tar.gz

After extracting, place the following files in
`rocket_6dof/third_party/imgui/`:

```
imconfig.h
imgui.cpp
imgui.h
imgui_draw.cpp
imgui_internal.h
imgui_tables.cpp
imgui_widgets.cpp
imstb_rectpack.h
imstb_textedit.h
imstb_truetype.h
backends/imgui_impl_glfw.cpp
backends/imgui_impl_glfw.h
backends/imgui_impl_opengl3.cpp
backends/imgui_impl_opengl3.h
backends/imgui_impl_opengl3_loader.h
```

The simplest path: clone or extract the release, then copy the listed
files. The Makefile compiles ImGui as part of the GUI build target.

### 2. ImPlot (plotting library for Dear ImGui)

- Home:    https://github.com/epezent/implot
- Version tested: **v0.16**
- Direct download:
  https://github.com/epezent/implot/archive/refs/tags/v0.16.tar.gz

After extracting, place the following files in
`rocket_6dof/third_party/implot/`:

```
implot.cpp
implot.h
implot_internal.h
implot_items.cpp
```

### 3. stb_image_write.h (single-header PNG writer)

- Home:    https://github.com/nothings/stb
- Version tested: **v1.16**
- Direct download:
  https://raw.githubusercontent.com/nothings/stb/master/stb_image_write.h

This is a single file. Save it as
`rocket_6dof/third_party/stb/stb_image_write.h`.

Used by the GUI's "Export as PNG" button on the Plots tab.

## What you get

### CLI binaries (built by `make all`)

| Binary               | Purpose                                                |
| -------------------- | ------------------------------------------------------ |
| `mission`            | Single mission run, given a JSON config                |
| `mission_mc`         | Monte Carlo sweep runner                               |
| `mission_sobol`      | Sobol sensitivity analysis runner                      |
| `mission_worker`     | Internal worker for parallel sweeps (called by fork)   |
| `build_mc_viewer`    | Bundles MC/Sobol CSVs into a self-contained HTML file  |
| `build_viewer`       | Older viewer (single-mission)                          |
| `sobol_merge`        | Merges multiple Sobol runs into one CSV                |
| `env_sweep`          | Environment-sweep regression test                      |
| `*_test`             | 16 regression tests covering each module               |

### GUI binary (`make mission_gui`)

A five-tab Dear ImGui application:

1. **Basic** — curated parameter sliders (Launch / Propulsion / Aerodynamics
   with file picker / Control + load relief / Termination). One-click Run.
2. **Advanced** — full editable JSON tree. Edit any numeric or string leaf.
3. **Plots** — five ImPlot panels (altitude, velocity, dynamic pressure,
   angle of attack, Mach vs time). Export to PNG. Load a second CSV and
   overlay it for comparison.
4. **Sweeps** — launch Monte Carlo / Sobol sweeps as subprocesses with a
   live progress bar. Build the interactive viewer. Inline ImPlot-rendered
   histograms and Sobol bar charts.
5. **Queue** — queue multiple mission configs to run sequentially. Per-item
   status, elapsed time, final altitude.

## Known issue: `aero_deck_test`

`aero_deck_test` references the original development layout where the
file `aero_deck_SLV.asc` lived at `../../ssr_unpack/aero_deck_SLV.asc`.
In this bundle the file ships in the working directory, so three of
the 16 regression checks will fail with "AeroDeck: cannot open
../../ssr_unpack/aero_deck_SLV.asc".

Workaround: either patch the test (replace the hard-coded path with
`./aero_deck_SLV.asc`), or create a symlink:

```sh
mkdir -p ../../ssr_unpack
ln -s "$PWD/aero_deck_SLV.asc" ../../ssr_unpack/aero_deck_SLV.asc
```

15 of 16 regression tests pass without the workaround. The aero_deck
functionality itself is fine; only the test scaffolding has hard-coded
paths.

## Verifying the build

Once `make all` and `make mission_gui` succeed:

```sh
# Run the launcher mission. Expected output includes:
#   Final altitude:   247089.7 m  (247.09 km)
#   Final dvbi:       9129.46 m/s
./mission mission_launcher.json | grep "Final"

# Run the 16 regression tests (15 if aero_deck_test path issue applies)
for t in env_sweep dyn_test cad_test prop_test aero_test intercept_test \
         control_test guidance_test ins_test gps_test startrack_test rcs_test \
         json_test distributions_test integration_test; do
    echo -n "$t: "
    ./$t | grep "Total failures"
done

# Launch the GUI (requires DISPLAY)
./mission_gui mission_launcher.json
```

## License

Application code: as provided in source files.
OSK-style kernel: as provided.
Third-party libraries (you download separately):
- Dear ImGui:        MIT License
- ImPlot:            MIT License
- stb_image_write:   Public Domain (Unlicense) / MIT (dual-licensed)

### References

The simulation models are based on Peter H. Zipfel, *Modeling and
Simulation of Aerospace Vehicle Dynamics*, Third Edition (AIAA
Education Series), and the CADAC++ code distributed with that book.
Zipfel's book and CADAC++ remain the property of their author and
publisher; consult the book for the authoritative model definitions.
See "About the simulation models" above.
