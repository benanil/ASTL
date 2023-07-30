@echo off

REM Compile the C++ code using g++
g++ -std=c++14 -w -O3 -mavx2 -fno-exceptions ^
Test/ASTL.cpp ^
Test/AdventOfCodeTests.cpp ^
Profiler.cpp ^
-o astl_test