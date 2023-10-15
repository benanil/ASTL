@echo off

REM Compile the C++ code using g++
g++ -std=c++14 -O1 -mavx2 -march=native -msse4.2 ^
-s -flto -fno-rtti -fno-stack-protector -fno-exceptions -fno-unwind-tables ^
Test/ASTL.cpp ^
GLTFParser.cpp ^
-o astl_test