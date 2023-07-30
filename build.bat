@echo off

REM Compile the C++ code using g++
g++ -std=c++14 -w -O3 -mavx2 -fno-exceptions ^
ASTL.cpp ^
AdventOfCodeTests.cpp ^
Profiler.cpp ^
-o astl_test