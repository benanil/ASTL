cmake_minimum_required(VERSION 3.12)
project(ASTL_Project)

# Set the C++ standard to C++14
set(CMAKE_CXX_STANDARD 14)
# Specify the custom build directory
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build)

# Add your source files here
set(SOURCES
    ASTL.cpp
    AdventOfCodeTests.cpp
    Profiler.cpp
)

# Add an executable target
add_executable(astl_test ${SOURCES})