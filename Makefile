TARGET = main
SRCS = src/bfgs.cpp src/main.cpp src/powell.cpp src/dfp.cpp src/math.cpp src/sobolseqgenerator.cpp src/global_min.cpp src/nesterov.cpp src/tinyxml2.cpp src/hessian_free.cpp
OBJS = $(SRCS:.cpp=.o)
CXX ?= gcc
CXXFLAGS ?= -pthread -Ofast -std=c++11 -Wall -Wextra -g -I"include"

.PHONY: all clean

all: $(TARGET)
$(TARGET): $(OBJS)
		$(CXX) -o $(TARGET) $(OBJS) $(CXXFLAGS)
 
.c.o:
		$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
		rm -rf $(TARGET) $(OBJS)