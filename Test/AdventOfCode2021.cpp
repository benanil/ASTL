
#include "../IO.hpp"
#include "../Array.hpp"
#include "../HashSet.hpp"
#include "../HashMap.hpp"
#include "../Math/Vector.hpp"

void Day1()
{
    char *text = ReadAllFile("Test/2021Day1.txt");
    const char *curr = text;

    short numbers[2001];
    short i = 0;

    while (*curr) 
        numbers[i++] = ParseNumber(curr);

    short numIncreased = 0;
    for (short j = 1; j < i; j++)
        numIncreased += numbers[j] > numbers[j-1];

    printf("AOC 2021 day 1 result: %i \n", numIncreased);

    numIncreased = 0;
    for (short j = 0; j < i - 2; j++)
    {
        short a = numbers[j + 0] + numbers[j + 1] + numbers[j + 2];
        short b = numbers[j + 1] + numbers[j + 2] + numbers[j + 3];
        numIncreased += b > a;
    }

    printf("AOC 2021 day 1 result2: %i \n", numIncreased);
    free(text);
}

void Day2()
{
    char* text = ReadAllFile("Test/2021Day2.txt");
    const char *curr = text;
    
    struct Command
    {
        int value      : 31;
        bool isForward : 1;
    };

    Command commands[1001];

    for (int i = 0; i < 1000; i++)
    {
        int wasForward = *curr == 'f';
        int wasUp = *curr == 'u';
        curr += wasForward * 7;
        curr += wasUp * 2;
        curr += *curr == 'd' * 4;
        int num = ParsePositiveNumber(curr);
        num = wasUp ? -num : num;
        
        commands[i].value = num;
        commands[i].isForward = wasForward;

        while (*curr && !IsLower(*curr))
            curr++;
    }
    
    int depthHorizontal[2]{};
    for (int i = 0; i < 1000; i++)
    {
        depthHorizontal[commands[i].isForward] += commands[i].value;
    }
    
    printf("depth: %i, horizontal: %i\n", depthHorizontal[0], depthHorizontal[1]);
    printf("AOC 2021 day 2 result: %i \n", depthHorizontal[0] * depthHorizontal[1]);

    depthHorizontal[0] = 0; 
    depthHorizontal[1] = 0; 
    int aim = 0;
    for (int i = 0; i < 1000; i++)
    {
        bool isForward = commands[i].isForward;
        aim += !isForward * commands[i].value;
        depthHorizontal[1] += isForward * commands[i].value;
        depthHorizontal[0] += isForward * aim * commands[i].value;
    }
    printf("AOC 2021 day 2 result2: %i \n", depthHorizontal[0] * depthHorizontal[1]);
    free(text);
}

// didn't understand the day 3


void Day4()
{
    char* text = ReadAllFile("Test/2021Day4.txt");
    const char* curr = text;
 
    struct alignas(32) Square
    {
        uint number  : 24;
        bool isBingo : 8;
    };

    struct Board
    {
        Square row[5][5];
    };

    StackArray<short, 101> movements;
    StackArray<Board, 101> boards;

    while (*curr != '\n')
    {
        movements.Add((short)ParseNumber(curr));
        curr += *curr != '\n';
    }
    curr++;

    while (*curr)
    {
        int index = boards.Size();
        boards.AddUninitialized(1);
        SmallMemSet(&boards[index], 0, sizeof(Board));

        for (int i = 0; i < 5; i++)
        {
            for (int j = 0; j < 5; j++)
            {
                int num = ParseNumber(curr);
                boards[index].row[i][j].number  = (uint)num;
                boards[index].row[i][j].isBingo = 0;
            }
            while (*curr && (*curr == '\n' || IsWhitespace(*curr)))
                curr++;
        }
    }

    StackHashSet<short, 101> winners;
    int i = 0, j = 0;

    for (; i < movements.Size(); i++)
    {
        short movement = movements[i];

        for (j = 0; j < boards.Size(); j++)
        {
            Board& board = boards[j];

            for (int k = 0; k < 5; k++)
                for (int l = 0; l < 5; l++)
                    board.row[k][l].isBingo |= board.row[k][l].number == movement;
            // check rows
            for (int k = 0; k < 5; k++)
            {
                short cnt = 0;
                for (int l = 0 ; l < 5; l++)
                    cnt += board.row[k][l].isBingo;
                
                // if (cnt == 5) goto result_final; // if part 1
                if (cnt == 5) winners.Insert(j); // if part 2
            }
            
            // check columns
            for (int k = 0; k < 5; k++)
            {
                short cnt = 0;
                for (int l = 0; l < 5; l++)
                    cnt += board.row[l][k].isBingo;
                // if (cnt == 5) goto result_final; // if part 1
                if (cnt == 5) winners.Insert(j); // if part 2
            }

            if (winners.Size() == boards.Size())  // if part 2
                goto result_final;
        }
    }
    printf("nobody wins!\n");
    return;
    result_final:
    {
        int unmarkedSum = 0;
        Board& board = boards[j];
        short movement = movements[i];
    
        for (int k = 0; k < 5; k++)
            for (int l = 0; l < 5; l++)
                unmarkedSum += (board.row[k][l].isBingo == 0) * board.row[k][l].number;
    
        printf("day4 result: %i \n", unmarkedSum * movement);
    }
}

template<typename T>
struct Line2D
{
    Vector2<T> begin, end;
};

// https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
template<typename T>
inline bool LineIntersection(Vector2<T> p1, Vector2<T> p2,
                             Vector2<T> p3, Vector2<T> p4, Vector2<T>* pOut)
{   
    // Store the values for fast access and easy
    // equations-to-code conversion
    T x1 = p1.x, x2 = p2.x, x3 = p3.x, x4 = p4.x;
    T y1 = p1.y, y2 = p2.y, y3 = p3.y, y4 = p4.y;
 
    T d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    // If d is zero, there is no intersection
    if (d == 0) return false;
 
    // Get the x and y
    T pre = (x1*y2 - y1*x2), post = (x3*y4 - y3*x4);
    T x = ( pre * (x3 - x4) - (x1 - x2) * post ) / d;
    T y = ( pre * (y3 - y4) - (y1 - y2) * post ) / d;
 
    // Check if the x and y coordinates are within both lines
    if ( x < MIN(x1, x2) | x > MAX(x1, x2) ||
        x < MIN(x3, x4) | x > MAX(x3, x4) ) return false;
    if ( y < MIN(y1, y2) | y > MAX(y1, y2) ||
        y < MIN(y3, y4) | y > MAX(y3, y4) ) return false;
 
    // Return the point of intersection
    pOut->x = x;
    pOut->y = y;
    return true;
}

// https://stackoverflow.com/questions/17692922/check-is-a-point-x-y-is-between-two-points-drawn-on-a-straight-line
inline bool PointOnLine(Vector2s A, Vector2s B, Vector2s C)
{
    //Check C is within the bounds of the line
    if (C.x > MAX(A.x, B.x) | C.x < MIN(A.x, B.x) || C.y < MIN(A.y, B.y) | C.y > MAX(A.y, B.y))
        return false;
    // Check for when AB is vertical
    if (A.x == B.x) return Abs(A.x - C.x) < 1;
    // Check for when AB is horizontal
    if (A.y == B.y) return Abs(A.y - C.y) < 1;

    // Check istance of the point form the line
    float distFromLine = (float)Abs(((B.x - A.x) * (A.y - C.y))-((A.x - C.x) * (B.y - A.y))) /
                                Sqrt((B.x - A.x) * (B.x - A.x) + (B.y - A.y) * (B.y - A.y));
    return distFromLine < 0.1;
}

void Day5()
{
    char* text = ReadAllFile("Test/2021Day5.txt");
    const char *curr = text; // "0,9 -> 5,9\n" "8,0 -> 0,8\n" "9,4 -> 3,4\n" "2,2 -> 2,1\n" "7,0 -> 7,4\n" "6,4 -> 2,0\n" "0,9 -> 2,9\n" "3,4 -> 1,4\n""0,0 -> 8,8\n""5,5 -> 8,2";
  
    Array<Line2D<short>> lines;
    HashMap<uint32, int32> map;
    int numInt = 0; // intersection count

    while (*curr)
    {
        lines.AddUninitialized(1);
        Line2D<short>& line = lines.Back();
        
        line.begin.x = ParseNumber(curr); curr++; // skip,
        line.begin.y = ParseNumber(curr); curr += 4; // skip " -> "
        line.end.x   = ParseNumber(curr); curr++; // skip,
        line.end.y   = ParseNumber(curr); curr += *curr == '\n';

        for (int i = 0; i < lines.Size()-1; i++)
        {
            Line2D<short>& linei = lines[i];
            // is line horizontal?
            if (line.begin.x == line.end.x && line.begin.x == linei.begin.x && linei.end.x == line.begin.x)
            {
                uint32 x = (uint32)line.begin.x;
                short start = MAX(MIN(line.begin.y, line.end.y), MIN(linei.begin.y, linei.end.y));
                short end   = MIN(MAX(line.begin.y, line.end.y), MAX(linei.begin.y, linei.end.y));
                while (start <= end)
                    numInt += map.Insert(x | (uint32(start++) << 16), 0)->value++ == 0;
            }
            // is line horizontal
            else if (line.begin.y == line.end.y && line.begin.y == linei.begin.y && linei.end.y == line.begin.y)
            {
                uint32 y = uint32(line.begin.y) << 16;
                short start = MAX(MIN(line.begin.x, line.end.x), MIN(linei.begin.x, linei.end.x));
                short end   = MIN(MAX(line.begin.x, line.end.x), MAX(linei.begin.x, linei.end.x));
                while (start <= end)
                    numInt += map.Insert(start++ | y, 0)->value++ == 0;
            }
            else 
            {
                // check is diagonal
                Vector2f lineDir  = ToVector2f(line.end - line.begin).Normalized();
                Vector2f lineiDir = ToVector2f(linei.end - linei.begin).Normalized();
                bool firstEnter = 0;
                equal_branch:
                if (AlmostEqual(lineDir.x, lineiDir.x) && AlmostEqual(lineDir.y, lineiDir.y))
                {
                    Vector2s end = Vector2s::LengthSquared(line.end - line.begin) > 
                                   Vector2s::LengthSquared(linei.end - linei.begin) ? line.end : linei.end;
                    if (PointOnLine(linei.begin, linei.end, line.begin))
                    {
                        Vector2s dir = ToVector2s(lineDir * 2.0f);
                        Vector2s pos = line.begin;

                        for (;pos != end; pos += dir)
                        {
                            uint32 key = pos.x | (uint32(pos.y) << 16);
                            numInt += map.Insert(key, 0)->value++ == 0;
                        }
                        numInt += map.Insert(end.x | (uint32(end.y) << 16), 0)->value++ == 0;
                    }
                    if (PointOnLine(linei.begin, linei.end, line.end))
                    {
                        Vector2s dir = -ToVector2s(lineDir * 2.0f);
                        Vector2s pos = linei.end;

                        for (;pos != end; pos += dir)
                        {
                            uint32 key = pos.x | (uint32(pos.y) << 16);
                            numInt += map.Insert(key, 0)->value++ == 0;
                        }
                        numInt += map.Insert(end.x | (uint32(end.y) << 16), 0)->value++ == 0;
                    }
                    continue;
                }
                
                lineDir = Vector2f::Normalize(ToVector2f(line.begin - line.end));

                if (!firstEnter)
                {
                    Swap(line.begin, line.end);
                    firstEnter = 1;
                    goto equal_branch;
                }
            }
            Vector2s p;
            // if this line intersection works only one point is intersected
            if (LineIntersection(line.begin, line.end, lines[i].begin, lines[i].end, &p))
            {
                uint32 key = uint32(p.x) | (uint32(p.y) << 16);
                KeyValuePair<uint32, int>* point = map.Insert(key, 0); // inserts if already not exist
                numInt += point->value++ == 0;
            }
        }
    }
    
    printf("day 5 result: %i", numInt);
    free(text);
} // 10787 low