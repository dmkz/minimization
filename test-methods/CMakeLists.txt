cmake_minimum_required(VERSION 3.1)

project(test-methods)

include_directories(../include)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -Ofast -std=c++11 -Wall -Wextra -g")

add_executable(main test_methods.cpp ../src/math.cpp ../src/nesterov.cpp ../src/hessian_free.cpp ../src/bfgs.cpp)