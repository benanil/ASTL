@echo off

REM Compile the C++ code using g++
REM GLTFParser.cpp
g++ -std=c++14 -O3 -mavx2 -march=native -msse4.2 ^
-s -fno-rtti -fno-stack-protector -fno-exceptions -static-libstdc++ -static-libgcc -fno-unwind-tables ^
Additional/Profiler.cpp ^
Test/ASTL.cpp ^
-o astl_test ^
-mno-needed

