# =====================================================================
# Makefile for examples/vmq_demo
#
# Standalone unit tests for the Vec, Mat, and Quat utility classes.
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
    $(OSK_DIR)/table.cpp  \
    $(OSK_DIR)/vec.cpp    \
    $(OSK_DIR)/mat.cpp    \
    $(OSK_DIR)/quat.cpp

KERNEL_OBJS := $(notdir $(KERNEL_SRCS:.cpp=.o))
KERNEL_HDRS := $(wildcard $(OSK_DIR)/*.h)

TARGET := vmq_demo

.PHONY: all run clean

all: $(TARGET)

$(TARGET): vmq_demo.o $(KERNEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: %.cpp $(KERNEL_HDRS)
	$(CXX) $(CXXFLAGS) -I$(OSK_DIR) -c $< -o $@

$(KERNEL_OBJS): %.o: $(OSK_DIR)/%.cpp $(KERNEL_HDRS)
	$(CXX) $(CXXFLAGS) -I$(OSK_DIR) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) *.o
