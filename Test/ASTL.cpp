
#include "../Profiler.hpp"
#include "AdventOfCode2021.cpp"
#include "AdventOfCode2022.cpp"
#include <stdio.h>

int main()
{   
    // AdventOfCodeTests(); // 2022
    BeginProfile();
    {
        TimeBlock("all");
        Day1();
        Day2();
        Day4();
    }
    EndAndPrintProfile();
    getchar();
    return 0;
}
  