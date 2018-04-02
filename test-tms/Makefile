TARGET = main
SRCS = test.cpp ../src/sobolseqgenerator.cpp ../src/NiederreiterBaseTwo.cpp
OBJS = $(SRCS:.cpp=.o)
CXX ?= gcc
CXXFLAGS ?= -Ofast -std=c++14 -Wall -Wextra -fmax-errors=2 -I"../include" -I"."

.PHONY: all clean

all: $(TARGET)
$(TARGET): $(OBJS)
		$(CXX) -o $(TARGET) $(OBJS) $(CXXFLAGS)
 
.c.o:
		$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
		rm -rf $(TARGET) $(OBJS)