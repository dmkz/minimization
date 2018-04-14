TARGET = main
SRCS =          main.cpp    ../src/global_min.cpp      	 ../src/hessian_free.cpp 
SRCS+=    ../src/dfp.cpp      ../src/nesterov.cpp ../src/NiederreiterBaseTwo.cpp  
SRCS+=   ../src/math.cpp ../src/StopCondition.cpp              	 ../src/bfgs.cpp  
SRCS+= ../src/powell.cpp 
OBJS = $(SRCS:.cpp=.o)
CXX ?= gcc
CXXFLAGS ?= -pthread -Ofast -std=c++14 -Wall -Wextra -fmax-errors=2 -I"../include"

.PHONY: all $(SUBDIRS) clean

all: $(TARGET) $(SUBDIRS)
$(TARGET): $(OBJS)
		$(CXX) -o $(TARGET) $(OBJS) $(CXXFLAGS)
$(SUBDIRS):
		$(MAKE) -C $@
 
.c.o:
		$(CXX) $(CXXFLAGS) -c $< -o $@

clean: $(SUBDIRS_CLEAN)
		rm -rf $(TARGET) $(OBJS)
$(SUBDIRS_CLEAN):
		$(MAKE) -C $(@:$(CLEAN_SUFFIX)=$(EMPTY)) clean