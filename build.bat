@echo off

REM Compile the C++ code using g++
g++ -std=c++17 -O3 -mavx2 -msse4.2 -fno-exceptions ^
Test/ASTL.cpp ^
Test/AdventOfCodeTests.cpp ^
Profiler.cpp ^
CPURayTrace.cpp ^
-o astl_test