TARGET = main
SRCS = test_methods.cpp ../src/math.cpp ../src/nesterov.cpp ../src/hessian_free.cpp ../src/bfgs.cpp ../src/dfp.cpp ../src/powell.cpp
OBJS = $(SRCS:.cpp=.o)
CXX ?= gcc
CXXFLAGS ?= -pthread -Ofast -std=c++11 -Wall -Wextra -g -I"../include"

.PHONY: all clean

all: $(TARGET)
$(TARGET): $(OBJS)
		$(CXX) -o $(TARGET) $(OBJS) $(CXXFLAGS)
 
.c.o:
		$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
		rm -rf $(TARGET) $(OBJS)