@echo off

REM Compile the C++ code using g++
REM GLTFParser.cpp
g++ -std=c++14 -O3 -mavx2 -march=native -msse4.2 -mwindows -luser32 -DAX_USE_WINDOW ^
-s -fno-rtti -fno-stack-protector -fno-exceptions -static-libstdc++ -static-libgcc -fno-unwind-tables ^
Test/ASTL.cpp ^
GLTFParser.cpp ^
Window.cpp ^
Renderer.cpp ^
-o astl_test ^
-lopengl32 -lgdi32 -mno-needed info.res icon.res
