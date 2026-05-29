# =====================================================================
# Makefile for examples/ex_app2  (Appendix 2 -- Adding an Integrator)
#
# Demonstrates the plug-in pattern by implementing RK2 as a State
# subclass and using it via the templated addIntegrator<State_rk2>()
# form.  Runs the test problem from PDF page A-20.
# =====================================================================

CXX      ?= g++
CXXFLAGS ?= -std=c++14 -Wall -Wextra -O2

OSK_DIR  := ../../osk

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
KERNEL_HDRS := $(wildcard $(OSK_DIR)/*.h)
APP_HDRS    := $(wildcard *.h)

TARGET := main

.PHONY: all run clean

all: $(TARGET)

$(TARGET): main.o $(KERNEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

main.o: main.cpp $(APP_HDRS) $(KERNEL_HDRS)
	$(CXX) $(CXXFLAGS) -I$(OSK_DIR) -c main.cpp -o main.o

$(KERNEL_OBJS): %.o: $(OSK_DIR)/%.cpp $(KERNEL_HDRS)
	$(CXX) $(CXXFLAGS) -I$(OSK_DIR) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) *.o
