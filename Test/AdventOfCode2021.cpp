
#include "../IO.hpp"
#include "../Algorithms.hpp"

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

#include "../Array.hpp"
#include "../HashSet.hpp"

void Day4()
{
    char* text = ReadAllFile("Test/2021Day3.txt");
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