# Building MERIDIAN on Linux

Tested on Ubuntu 24.04 with g++ 13. Should work on any reasonably
modern distro with a C++14-capable compiler.

## 1. Install build dependencies

### Ubuntu / Debian

```sh
sudo apt update
sudo apt install build-essential libglfw3-dev libgl1-mesa-dev
```

`build-essential` provides g++, make, and the standard libraries.
`libglfw3-dev` and `libgl1-mesa-dev` are only needed for the GUI build
(`mission_gui`). The CLI binaries (`make all`) don't need them.

#### WSL users (Ubuntu/Debian inside Windows Subsystem for Linux)

The GUI's "Open in browser" button (Sweeps tab) launches the
generated `mc_viewer.html` / `mission_viewer.html` via
`explorer.exe`, which uses the Windows default browser. The launcher
chain tries `explorer.exe` on PATH first, falls back to the absolute
path `/mnt/c/Windows/explorer.exe` if PATH lookup fails (e.g. from
XRDP-launched terminals where WSL's PATH-injection didn't propagate),
then `xdg-open` and `wslview`. No extra packages are required for
the default Ubuntu WSL terminal.

If you launch the GUI from a terminal where Windows binaries aren't
on PATH (typical in XFCE-over-XRDP sessions), the absolute-path
fallback handles it transparently. You can also fix PATH for that
session by adding the following to `~/.profile` (which is sourced
by XRDP login):

```sh
if [ -d /mnt/c/Windows ]; then
    export PATH="$PATH:/mnt/c/Windows:/mnt/c/Windows/System32"
fi
```

### Fedora / RHEL

```sh
sudo dnf install gcc-c++ make glfw-devel mesa-libGL-devel
```

### Arch

```sh
sudo pacman -S base-devel glfw mesa
```

## 2. Set up third-party libraries

The `third_party/` directory ships empty. Populate it with the three
libraries listed in the main README:

### Dear ImGui (v1.91.6)

```sh
cd rocket_6dof
curl -L -o /tmp/imgui-1.91.6.tar.gz \
    https://github.com/ocornut/imgui/archive/refs/tags/v1.91.6.tar.gz
tar -xzf /tmp/imgui-1.91.6.tar.gz -C /tmp

# Copy the files the Makefile expects
cp /tmp/imgui-1.91.6/imconfig.h          third_party/imgui/
cp /tmp/imgui-1.91.6/imgui.cpp           third_party/imgui/
cp /tmp/imgui-1.91.6/imgui.h             third_party/imgui/
cp /tmp/imgui-1.91.6/imgui_draw.cpp      third_party/imgui/
cp /tmp/imgui-1.91.6/imgui_internal.h    third_party/imgui/
cp /tmp/imgui-1.91.6/imgui_tables.cpp    third_party/imgui/
cp /tmp/imgui-1.91.6/imgui_widgets.cpp   third_party/imgui/
cp /tmp/imgui-1.91.6/imstb_rectpack.h    third_party/imgui/
cp /tmp/imgui-1.91.6/imstb_textedit.h    third_party/imgui/
cp /tmp/imgui-1.91.6/imstb_truetype.h    third_party/imgui/
mkdir -p third_party/imgui/backends
cp /tmp/imgui-1.91.6/backends/imgui_impl_glfw.cpp     third_party/imgui/backends/
cp /tmp/imgui-1.91.6/backends/imgui_impl_glfw.h       third_party/imgui/backends/
cp /tmp/imgui-1.91.6/backends/imgui_impl_opengl3.cpp  third_party/imgui/backends/
cp /tmp/imgui-1.91.6/backends/imgui_impl_opengl3.h    third_party/imgui/backends/
cp /tmp/imgui-1.91.6/backends/imgui_impl_opengl3_loader.h third_party/imgui/backends/
```

### ImPlot (v0.16)

```sh
curl -L -o /tmp/implot-0.16.tar.gz \
    https://github.com/epezent/implot/archive/refs/tags/v0.16.tar.gz
tar -xzf /tmp/implot-0.16.tar.gz -C /tmp

cp /tmp/implot-0.16/implot.cpp           third_party/implot/
cp /tmp/implot-0.16/implot.h             third_party/implot/
cp /tmp/implot-0.16/implot_internal.h    third_party/implot/
cp /tmp/implot-0.16/implot_items.cpp     third_party/implot/
```

### stb_image_write (single header)

```sh
curl -L -o third_party/stb/stb_image_write.h \
    https://raw.githubusercontent.com/nothings/stb/master/stb_image_write.h
```

## 3. Build

```sh
# CLI binaries (no GUI deps required)
make all

# GUI binary (requires GLFW + OpenGL)
make mission_gui

# Both
make all mission_gui
```

`make` is single-threaded by default; for a faster build:

```sh
make -j$(nproc) all mission_gui
```

## 4. Verify

```sh
# Run the sample launcher mission. Expected output:
#   Final altitude:   247089.7 m  (247.09 km)
#   Final dvbi:       9129.46 m/s
./mission mission_launcher.json | grep Final

# Run the GUI
./mission_gui mission_launcher.json
```

## Troubleshooting

**`fatal error: GLFW/glfw3.h: No such file or directory`**
Install `libglfw3-dev`. On macOS, `brew install glfw`.

**`undefined reference to glClear` etc.**
GUI build needs `-lGL`. The Makefile includes this already; if you're
on a system that splits GL into multiple libs, you may need to add
`-lOpenGL` or similar.

**`stb_image_write.h: No such file or directory`**
You haven't placed `stb_image_write.h` in `third_party/stb/`. See
"Set up third-party libraries" above.

**`mission_gui` opens a black window**
Check that your OpenGL drivers support GL 3.0+. On a VM or container
without GPU passthrough, try `LIBGL_ALWAYS_SOFTWARE=1 ./mission_gui`
to use software rendering.

**Three `aero_deck_test` failures: "cannot open ../../ssr_unpack/..."**
This is a pre-existing path issue in the test. See the main README
for the symlink workaround.

## Headless / CI

The CLI build (`make all`) works without any display. The GUI build
requires GLFW dev headers but the binary can be launched under `xvfb`
for screenshot-based testing:

```sh
sudo apt install xvfb imagemagick
Xvfb :99 -screen 0 1400x900x24 &
DISPLAY=:99 ./mission_gui mission_launcher.json &
sleep 2
DISPLAY=:99 import -window root /tmp/screenshot.png
```
