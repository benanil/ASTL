
#include "../Additional/Profiler.hpp"
#include "AdventOfCode2021.cpp"
#include "AdventOfCode2022.cpp"
#include "AdventOfCode2023.cpp"

int main()
{
    BeginProfile();

    AOC2023Day3();
    AOC2023Day4Part1();

    EndAndPrintProfile();
    // 
    // AdventOfCodeTests21();
    // AdventOfCodeTests22();
    // 
    getchar();
    return 0;
}
