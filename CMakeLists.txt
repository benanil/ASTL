cmake_minimum_required(VERSION 3.12)
project(ASTL_Project)

# Set the C++ standard to C++14
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_WIN32_EXECUTABLE ON)

# Specify the custom build directory
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build)

# Add your source files here
set(SOURCES
    Test/ASTL.cpp
    Profiler.cpp
)

include_directories(${CMAKE_SOURCE_DIR})

# Add an executable target
add_executable(astl_test ${SOURCES})

target_link_libraries(astl_test 
    user32
    gdi32
    opengl32
)

# Copy AOC11.txt and AOC12.txt to build/Test folder during the build process
configure_file(Test/AOC12.txt ${CMAKE_BINARY_DIR}/Test/AOC12.txt COPYONLY)
configure_file(Test/AOC15.txt ${CMAKE_BINARY_DIR}/Test/AOC15.txt COPYONLY)
configure_file(Test/AOC17.txt ${CMAKE_BINARY_DIR}/Test/AOC17.txt COPYONLY)
configure_file(Test/AOC22.txt ${CMAKE_BINARY_DIR}/Test/AOC22.txt COPYONLY)

configure_file(Test/AOC22.txt ${CMAKE_BINARY_DIR}/Test/2021Day1.txt COPYONLY)
configure_file(Test/AOC22.txt ${CMAKE_BINARY_DIR}/Test/2021Day2.txt COPYONLY)
configure_file(Test/AOC22.txt ${CMAKE_BINARY_DIR}/Test/2021Day4.txt COPYONLY)
configure_file(Test/AOC22.txt ${CMAKE_BINARY_DIR}/Test/2021Day5.txt COPYONLY)