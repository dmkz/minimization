TARGET = main
SRCS = test_NiederreiterBaseTwo.cpp ../src/NiederreiterBaseTwo.cpp
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