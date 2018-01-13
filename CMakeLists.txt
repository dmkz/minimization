cmake_minimum_required(VERSION 3.1)

project(minimization)

include_directories(include)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -Ofast -std=c++11 -Wall -Wextra -g")

add_executable(main src/main.cpp src/sobolseqgenerator.cpp src/math.cpp src/nesterov.cpp src/hessian_free.cpp src/bfgs.cpp src/global_min.cpp)
add_executable(bin/tests/tms_net_tests src/tmsnettests.cpp src/sobolseqgenerator.cpp src/tinyxml2.cpp)

# add_subdirectory (tests)
# enable_testing ()
# add_test (NAME Factorial COMMAND Test)