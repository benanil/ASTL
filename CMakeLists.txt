cmake_minimum_required(VERSION 3.12)
project(ASTL_Project)

# Set the C++ standard to C++14
set(CMAKE_CXX_STANDARD 14)

# Specify the custom build directory
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build)

# Add your source files here
set(SOURCES
    GLTFParser.cpp
    Test/ASTL.cpp
    Profiler.cpp
    Renderer.cpp
    Window.cpp
)

include_directories(${CMAKE_SOURCE_DIR})

# Add an executable target
add_executable(astl_test ${SOURCES})

# Copy AOC11.txt and AOC12.txt to build/Test folder during the build process
configure_file(Test/AOC12.txt ${CMAKE_BINARY_DIR}/Test/AOC12.txt COPYONLY)
configure_file(Test/AOC15.txt ${CMAKE_BINARY_DIR}/Test/AOC15.txt COPYONLY)
configure_file(Test/AOC17.txt ${CMAKE_BINARY_DIR}/Test/AOC17.txt COPYONLY)
configure_file(Test/AOC22.txt ${CMAKE_BINARY_DIR}/Test/AOC22.txt COPYONLY)

file(COPY ${CMAKE_SOURCE_DIR}/Meshes DESTINATION ${CMAKE_BINARY_DIR}/Meshes)
file(COPY ${CMAKE_SOURCE_DIR}/Textures DESTINATION ${CMAKE_BINARY_DIR}/Textures)
file(COPY ${CMAKE_SOURCE_DIR}/Shaders DESTINATION ${CMAKE_BINARY_DIR}/Shaders)

configure_file(Test/AOC22.txt ${CMAKE_BINARY_DIR}/Test/2021Day1.txt COPYONLY)
configure_file(Test/AOC22.txt ${CMAKE_BINARY_DIR}/Test/2021Day2.txt COPYONLY)
configure_file(Test/AOC22.txt ${CMAKE_BINARY_DIR}/Test/2021Day4.txt COPYONLY)
configure_file(Test/AOC22.txt ${CMAKE_BINARY_DIR}/Test/2021Day5.txt COPYONLY)