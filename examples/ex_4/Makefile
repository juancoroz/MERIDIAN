# =====================================================================
# Makefile for examples/ex_4  (Section 4.4 -- Multiple Connected Models)
#
#     make            # build
#     make run        # build and run
#     make clean
# =====================================================================

CXX      ?= g++
CXXFLAGS ?= -std=c++14 -Wall -Wextra -O2

OSK_DIR  := ../../osk

KERNEL_SRCS := \
    $(OSK_DIR)/state.cpp  \
    $(OSK_DIR)/block.cpp  \
    $(OSK_DIR)/sim.cpp    \
    $(OSK_DIR)/filer.cpp  \
    $(OSK_DIR)/table.cpp

APP_SRCS := main.cpp autopilot.cpp missile.cpp

APP_OBJS    := $(APP_SRCS:.cpp=.o)
KERNEL_OBJS := $(notdir $(KERNEL_SRCS:.cpp=.o))
OBJS        := $(APP_OBJS) $(KERNEL_OBJS)

KERNEL_HDRS := $(wildcard $(OSK_DIR)/*.h)
APP_HDRS    := $(wildcard *.h)

TARGET := main

.PHONY: all run clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

%.o: %.cpp $(APP_HDRS) $(KERNEL_HDRS)
	$(CXX) $(CXXFLAGS) -I$(OSK_DIR) -c $< -o $@

$(KERNEL_OBJS): %.o: $(OSK_DIR)/%.cpp $(KERNEL_HDRS)
	$(CXX) $(CXXFLAGS) -I$(OSK_DIR) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) $(OBJS)
