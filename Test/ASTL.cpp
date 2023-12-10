
#include "../Additional/Profiler.hpp"
#include "AdventOfCode2021.cpp"
#include "AdventOfCode2022.cpp"
#include "AdventOfCode2023.cpp"

int main()
{
    BeginProfile();
    
    // AdventOfCodeTests21();
    // AdventOfCodeTests22();
    AOC2023Day1Part1();
    AOC2023Day1Part2();
    AOC2023Day3();
    AOC2023Day4Part1();

    EndAndPrintProfile();
    getchar();
    return 0;
}
