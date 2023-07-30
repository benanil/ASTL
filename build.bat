@echo off

REM Compile the C++ code using g++
g++ -std=c++17 -O3 -mavx2 -lm -fno-exceptions ^
Test/ASTL.cpp ^
Test/AdventOfCodeTests.cpp ^
Profiler.cpp ^
CPURayTrace.cpp ^
-o astl_test