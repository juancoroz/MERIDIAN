CXX      ?= g++
CXXFLAGS ?= -std=c++14 -Wall -Wextra -O2

KERNEL_SRC = osk/state.cpp osk/block.cpp osk/sim.cpp \
             osk/filer.cpp osk/table.cpp \
             osk/vec.cpp osk/mat.cpp osk/quat.cpp

all: examples/hello examples/shm examples/ex_1/main \
     examples/ex_2/main examples/ex_3/main examples/ex_4/main \
     examples/ex_5/main \
     examples/ex_6_1/main examples/ex_6_2/main \
     examples/ex_6_3/main examples/ex_6_4/main \
     examples/ex_app2/main \
     examples/ex_1_io/main examples/util_demo/util_demo \
     examples/vmq_demo/vmq_demo

examples/hello: examples/hello.cpp $(KERNEL_SRC) osk/*.h
	$(CXX) $(CXXFLAGS) examples/hello.cpp $(KERNEL_SRC) -o $@

examples/shm: examples/shm.cpp $(KERNEL_SRC) osk/*.h
	$(CXX) $(CXXFLAGS) examples/shm.cpp $(KERNEL_SRC) -o $@

examples/ex_1/main: examples/ex_1/main.cpp examples/ex_1/model.cpp \
                    examples/ex_1/model.h $(KERNEL_SRC) osk/*.h
	$(CXX) $(CXXFLAGS) examples/ex_1/main.cpp examples/ex_1/model.cpp \
	    $(KERNEL_SRC) -o $@

examples/ex_2/main: examples/ex_2/main.cpp examples/ex_2/model.cpp \
                    examples/ex_2/model.h $(KERNEL_SRC) osk/*.h
	$(CXX) $(CXXFLAGS) examples/ex_2/main.cpp examples/ex_2/model.cpp \
	    $(KERNEL_SRC) -o $@

examples/ex_3/main: examples/ex_3/main.cpp examples/ex_3/model.cpp \
                    examples/ex_3/model.h $(KERNEL_SRC) osk/*.h
	$(CXX) $(CXXFLAGS) examples/ex_3/main.cpp examples/ex_3/model.cpp \
	    $(KERNEL_SRC) -o $@

examples/ex_4/main: examples/ex_4/main.cpp \
                    examples/ex_4/autopilot.cpp examples/ex_4/autopilot.h \
                    examples/ex_4/missile.cpp examples/ex_4/missile.h \
                    $(KERNEL_SRC) osk/*.h
	$(CXX) $(CXXFLAGS) examples/ex_4/main.cpp \
	    examples/ex_4/autopilot.cpp examples/ex_4/missile.cpp \
	    $(KERNEL_SRC) -o $@

examples/ex_5/main: examples/ex_5/main.cpp \
                    examples/ex_5/autopilot.cpp examples/ex_5/autopilot.h \
                    examples/ex_5/missile.cpp examples/ex_5/missile.h \
                    $(KERNEL_SRC) osk/*.h
	$(CXX) $(CXXFLAGS) examples/ex_5/main.cpp \
	    examples/ex_5/autopilot.cpp examples/ex_5/missile.cpp \
	    $(KERNEL_SRC) -o $@

examples/ex_6_1/main: examples/ex_6_1/main.cpp examples/ex_6_1/model.cpp \
                      examples/ex_6_1/model.h $(KERNEL_SRC) osk/*.h
	$(CXX) $(CXXFLAGS) examples/ex_6_1/main.cpp examples/ex_6_1/model.cpp \
	    $(KERNEL_SRC) -o $@

examples/ex_6_2/main: examples/ex_6_2/main.cpp examples/ex_6_2/model.cpp \
                      examples/ex_6_2/model.h $(KERNEL_SRC) osk/*.h
	$(CXX) $(CXXFLAGS) examples/ex_6_2/main.cpp examples/ex_6_2/model.cpp \
	    $(KERNEL_SRC) -o $@

examples/ex_6_3/main: examples/ex_6_3/main.cpp examples/ex_6_3/model.cpp \
                      examples/ex_6_3/model.h $(KERNEL_SRC) osk/*.h
	$(CXX) $(CXXFLAGS) examples/ex_6_3/main.cpp examples/ex_6_3/model.cpp \
	    $(KERNEL_SRC) -o $@

examples/ex_6_4/main: examples/ex_6_4/main.cpp \
                      examples/ex_6_4/autopilot.cpp examples/ex_6_4/autopilot.h \
                      examples/ex_6_4/missile.cpp examples/ex_6_4/missile.h \
                      $(KERNEL_SRC) osk/*.h
	$(CXX) $(CXXFLAGS) examples/ex_6_4/main.cpp \
	    examples/ex_6_4/autopilot.cpp examples/ex_6_4/missile.cpp \
	    $(KERNEL_SRC) -o $@

examples/ex_1_io/main: examples/ex_1_io/main.cpp examples/ex_1_io/model.cpp \
                       examples/ex_1_io/model.h $(KERNEL_SRC) osk/*.h
	$(CXX) $(CXXFLAGS) examples/ex_1_io/main.cpp examples/ex_1_io/model.cpp \
	    $(KERNEL_SRC) -o $@

examples/ex_app2/main: examples/ex_app2/main.cpp examples/ex_app2/model.h \
                       examples/ex_app2/state_rk2.h \
                       $(KERNEL_SRC) osk/*.h
	$(CXX) $(CXXFLAGS) examples/ex_app2/main.cpp $(KERNEL_SRC) -o $@

examples/util_demo/util_demo: examples/util_demo/util_demo.cpp \
                              $(KERNEL_SRC) osk/*.h
	$(CXX) $(CXXFLAGS) examples/util_demo/util_demo.cpp $(KERNEL_SRC) -o $@

examples/vmq_demo/vmq_demo: examples/vmq_demo/vmq_demo.cpp \
                            $(KERNEL_SRC) osk/*.h
	$(CXX) $(CXXFLAGS) examples/vmq_demo/vmq_demo.cpp $(KERNEL_SRC) -o $@

clean:
	rm -f examples/hello examples/shm \
	      examples/ex_1/main examples/ex_2/main examples/ex_3/main \
	      examples/ex_4/main examples/ex_5/main \
	      examples/ex_6_1/main examples/ex_6_2/main \
	      examples/ex_6_3/main examples/ex_6_4/main \
	      examples/ex_app2/main \
	      examples/ex_1_io/main \
	      examples/util_demo/util_demo \
	      examples/vmq_demo/vmq_demo

.PHONY: all clean
