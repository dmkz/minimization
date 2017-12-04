cmake_minimum_required(VERSION 3.1)

project(minimization)

include_directories(include)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror")

add_executable(main src/main.cpp)
add_executable(main src/tms.cpp)