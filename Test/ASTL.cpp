
#include "../Additional/Profiler.hpp"
#include "AdventOfCode2021.cpp"
#include "AdventOfCode2022.cpp"

int main()
{
    BeginProfile();
    
    AdventOfCodeTests21();
    AdventOfCodeTests22();
    
    EndAndPrintProfile();
    getchar();
    return 0;
}
