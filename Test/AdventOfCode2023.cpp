
#include "../Algorithms.hpp"
#include "../Math/Vector.hpp"
#include "../IO.hpp"
#include "../HashMap.hpp"
#include "../HashSet.hpp"
#include "../String.hpp"
#include "../BitSet.hpp"

void AOC2023Day1Part1()
{
    ScopedText text(ReadAllFile("Test/2023Day1.txt"));
    const char *curr = text.text;
    int total = 0;
    while (*curr)
    {
        while (!IsNumber(*curr)) curr++; 

        int firstDigit = *curr++ - '0', lastDigit = -1;
        for (; *curr > '\n'; curr++)
            if (IsNumber(*curr)) lastDigit = *curr - '0';

        if (lastDigit == -1) 
            lastDigit = firstDigit;
        total += firstDigit * 10 + lastDigit;
    }
    printf("\nresult: %i\n", total);
}

inline int StrCmp(const char* a, const char* b)
{
    bool equal = true;
    while (*b)
        equal &= *a++ == *b++;
    return equal;
}

void AOC2023Day1Part2()
{
    ScopedText text(ReadAllFile("Test/2023Day1.txt"));
    const char* curr = text.text;

    int total = 0, firstDigit = -1, lastDigit = -1, findedNumber = -1;
    while (true)
    {
        if (IsNumber(*curr))
            findedNumber = *curr - '0';
        else if (*curr == '\n' || *curr == '\0')
        {
            if (firstDigit == -1) break;
            total += firstDigit * 10 + lastDigit;
            firstDigit = -1, lastDigit = -1;
        }
        else
        // o,t,f,s,e,n
        switch (*curr)
        {
            case 'o': 
                if (StrCmp(curr, "one")) findedNumber = 1;
                break;
            case 't': 
                if (StrCmp(curr, "two"))   findedNumber = 2;
                if (StrCmp(curr, "three")) findedNumber = 3;
                break;
            case 'f':
                if (StrCmp(curr, "four")) findedNumber = 4;
                if (StrCmp(curr, "five")) findedNumber = 5;
                break;
            case 's':
                if (StrCmp(curr, "six"))   findedNumber = 6;
                if (StrCmp(curr, "seven")) findedNumber = 7;
                break;
            case 'e': if (StrCmp(curr, "eight")) findedNumber = 8; break;
            case 'n': if (StrCmp(curr, "nine"))  findedNumber = 9; break;
            default: break;
        }

        if (findedNumber != -1)
        {
            if (firstDigit == -1) firstDigit = findedNumber;
            lastDigit = findedNumber;
            findedNumber = -1;
        }
        curr++;
    }
    printf("second part result: %i \n", total);// 54770
}

void AOC2023Day3()
{
    ScopedText text(ReadAllFile("Test/2023Day3.txt"));
    const char* curr = text.text;

    HashMap<Vector2i, long> map(64u);
    HashMap<Vector2i, char> symbolSet(32u);
    Vector2i currPos = MakeVec2<int>(0,0);
    
    // Generate & Parse map
    while (*curr)
    {
        if (IsNumber(*curr))
        {
            long number = (long)ParsePositiveNumber(curr);
            int numDigits = Log10((unsigned)number) + 1;
            
            for (int i = 0; i < numDigits; i++, currPos.x++)
                map.Insert(currPos, number);
            
            currPos.x--;
            curr--;
        }
        else if (*curr == '\n')
            currPos.y++, currPos.x = -1;
        // detect if its symbol
        else if (!IsWhitespace(*curr) && *curr != '.')
            symbolSet.Insert(currPos, *curr);
        
        currPos.x++;
        curr++;
    }

    long find[2]={ -1,-1 };
    long total = 0, gearTotal = 0, numFind = 0;
    
    const auto checkAndAddFn = [&](Vector2i pos, Vector2i dir) -> void
    {
        if (map.Contains(pos + dir)) 
        {
            pos += dir;
            int number = map[pos];
            int numDigits = Log10((unsigned)number) + 1;
            find[numFind++] = number;
            total += number;
            
            // find starting point
            while (map.Contains(MakeVec2(pos.x-1, pos.y)))
                pos.x--;

            // erase the number since we found it
            for (int i = 0; i < numDigits; i++, pos.x++)
                map.Erase(pos);
        }
    };

    const Vector2i directions[8] = {
        { 1, 0}, {-1, 0}, { 0, 1}, { 0,-1}, { 1, 1}, {-1,-1}, { 1,-1}, {-1, 1}
    };
    // search around symbols, and find unique numbers that are intersects
    for (KeyValuePair<Vector2i, char>& pair : symbolSet)
    {
        find[0] = find[1] = numFind = 0;

        for (int i = 0; i < 8; i++)
            checkAndAddFn(pair.key, directions[i]);

        if (numFind == 2 && pair.value == '*')
            gearTotal += find[0] * find[1];
    }
    
    printf("day 3 result: a=%ld, b=%ld \n", total, gearTotal); // a=553825, b=93994191 
    ASSERT(total == 553825 && gearTotal == 93994191l);
}
            
void AOC2023Day4Part1()
{
    ScopedText text(ReadAllFile("Test/2023Day4.txt"));
    const char* curr = text.text;

    bool winingNumbers[100]{};
    long totalPoint = 0;

    for (int i = 0; i < 218; i++)
    {
        ParsePositiveNumber(curr); // skip game number
        MemsetZero(winingNumbers, sizeof(winingNumbers));

        for (int j = 0; j < 10; j++)
            winingNumbers[ParsePositiveNumber(curr)] = true;
    
        long numWin = 0;
        for (int j = 0; j < 25; j++)
            numWin += winingNumbers[ParsePositiveNumber(curr)];
                
        // const int map[25] = {0, 1, 2, 4, 8, 16, ...};
        totalPoint += (1l << numWin) >> 1; // //map[numWin]
    }
    printf("day4 result: %ld", totalPoint); // 25004
}
