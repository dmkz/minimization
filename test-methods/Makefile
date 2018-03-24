TARGET = main
SRCS = test_methods.cpp ../src/StopCondition.cpp ../src/math.cpp ../src/powell.cpp ../src/hessian_free.cpp ../src/nesterov.cpp ../src/bfgs.cpp ../src/dfp.cpp
OBJS = $(SRCS:.cpp=.o)
CXX ?= gcc
CXXFLAGS ?= -pthread -Ofast -std=c++14 -Wall -Wextra -fmax-errors=2 -I"../include"

.PHONY: all clean

all: $(TARGET)
$(TARGET): $(OBJS)
		$(CXX) -o $(TARGET) $(OBJS) $(CXXFLAGS)
 
.c.o:
		$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
		rm -rf $(TARGET) $(OBJS)