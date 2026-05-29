# =====================================================================
# Makefile for examples/ex_1  (Section 4.1 -- Simple Airframe/Autopilot)
#
# Builds an executable `main` from main.cpp + model.cpp linked against
# the OSK kernel.  Run from this directory:
#
#     make            # build
#     make run        # build and run
#     make clean      # remove generated files
#
# Override the compiler or flags from the command line:
#     make CXX=clang++ CXXFLAGS="-std=c++17 -O3"
# =====================================================================

# ---- toolchain --------------------------------------------------------
CXX      ?= g++
CXXFLAGS ?= -std=c++14 -Wall -Wextra -O2

# ---- paths ------------------------------------------------------------
# Kernel lives in ../../osk relative to this example directory, which
# matches the report's convention (Section 2.3.2.1.2):
#     dos> bcc32 hello.cpp ../osk/state.cpp ../osk/block.cpp ../osk/sim.cpp
#     lin> g++   hello.cpp ../osk/state.cpp ../osk/block.cpp ../osk/sim.cpp -o hello
OSK_DIR  := ../../osk

# Kernel sources.  state/block/sim are the Section 3 kernel.  filer
# and table are the Appendix 1 utilities.  ex_1 itself only needs the
# Section 3 ones, but linking all of them costs nothing and keeps this
# makefile consistent with the other examples.
KERNEL_SRCS := \
    $(OSK_DIR)/state.cpp  \
    $(OSK_DIR)/block.cpp  \
    $(OSK_DIR)/sim.cpp    \
    $(OSK_DIR)/filer.cpp  \
    $(OSK_DIR)/table.cpp

# Example sources.
APP_SRCS := main.cpp model.cpp

# Object files: one .o per source, kept in this directory.  We
# substitute the kernel path so kernel .o files land here too rather
# than polluting ../../osk.
APP_OBJS    := $(APP_SRCS:.cpp=.o)
KERNEL_OBJS := $(notdir $(KERNEL_SRCS:.cpp=.o))
OBJS        := $(APP_OBJS) $(KERNEL_OBJS)

# Headers that should trigger a rebuild if changed.
KERNEL_HDRS := $(wildcard $(OSK_DIR)/*.h)
APP_HDRS    := $(wildcard *.h)

TARGET := main

# ---- rules ------------------------------------------------------------
.PHONY: all run clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

# Compile rule for sources living next to this makefile.
%.o: %.cpp $(APP_HDRS) $(KERNEL_HDRS)
	$(CXX) $(CXXFLAGS) -I$(OSK_DIR) -c $< -o $@

# Compile rule for the kernel sources in $(OSK_DIR).  We can't use
# `%.o: $(OSK_DIR)/%.cpp` directly because that would conflict with
# the rule above on a stem match, so each kernel .o gets its own
# explicit rule via a static pattern.
$(KERNEL_OBJS): %.o: $(OSK_DIR)/%.cpp $(KERNEL_HDRS)
	$(CXX) $(CXXFLAGS) -I$(OSK_DIR) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) $(OBJS)
