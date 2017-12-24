cmake_minimum_required(VERSION 3.1)

project(minimization)

include_directories(include)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -g")

add_executable(main src/main.cpp)

add_subdirectory (tests)
enable_testing ()
add_test (NAME Factorial COMMAND Test)
