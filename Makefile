TARGET = main
SRCS = src/iteration_object.cpp src/bfgs.cpp src/main.cpp src/powell.cpp src/dfp.cpp src/math.cpp src/sobolseqgenerator.cpp src/global_min.cpp src/nesterov.cpp src/tinyxml2.cpp src/hessian_free.cpp
OBJS = $(SRCS:.cpp=.o)
CXX ?= gcc
CXXFLAGS ?= -pthread -Ofast -std=c++14 -Wall -Wextra -fmax-errors=2 -I"include"
SUBDIRS := test-methods

.PHONY: all $(SUBDIRS) clean

all: $(TARGET) $(SUBDIRS)
$(TARGET): $(OBJS)
		$(CXX) -o $(TARGET) $(OBJS) $(CXXFLAGS)
$(SUBDIRS):
		$(MAKE) -C $@
 
.c.o:
		$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
		rm -rf $(TARGET) $(OBJS)
