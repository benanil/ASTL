@echo off

REM Compile the C++ code using g++
REM GLTFParser.cpp
g++ -std=c++14 -O1 -mavx2 -march=native -msse4.2 -mwindows -luser32 ^
-s -flto -fno-rtti -fno-stack-protector -fno-exceptions -fno-unwind-tables ^
Test/ASTL.cpp ^
GLTFParser.cpp ^
Window.cpp ^
Renderer.cpp ^
-o astl_test ^
-lopengl32 -lgdi32 info.res icon.res
