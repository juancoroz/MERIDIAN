# =====================================================================
# Makefile for examples/sync_test
#
# Two tests that verify the kernel respects the source-code ordering
# of operations inside update(), as required by Section 4.2.2 of the
# CMD User's Guide.
#
#     make            # build both tests
#     make run        # build and run both
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

KERNEL_OBJS := $(notdir $(KERNEL_SRCS:.cpp=.o))
KERNEL_HDRS := $(wildcard $(OSK_DIR)/*.h)

TARGETS := order_test sync_test

.PHONY: all run clean

all: $(TARGETS)

order_test: order_test.o $(KERNEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

sync_test:  sync_test.o  $(KERNEL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: %.cpp $(KERNEL_HDRS)
	$(CXX) $(CXXFLAGS) -I$(OSK_DIR) -c $< -o $@

$(KERNEL_OBJS): %.o: $(OSK_DIR)/%.cpp $(KERNEL_HDRS)
	$(CXX) $(CXXFLAGS) -I$(OSK_DIR) -c $< -o $@

run: $(TARGETS)
	@echo "=== order_test ===" && ./order_test
	@echo
	@echo "=== sync_test ===" && ./sync_test

clean:
	rm -f $(TARGETS) *.o
