
#include <stdio.h>

#include "ASTL/Algorithms.hpp"
#include "ASTL/Array.hpp"
#include "ASTL/IO.hpp"
#include "ASTL/RedBlackTree.hpp"

static int Day1_1()
{
    ScopedText data = ReadAllText("AdventInput1.txt");
    const char* curr = data;
    
    Array<int> leftList;
    Array<int> rightList;

    while (*curr) {
        leftList.Add(ParsePositiveNumber(curr));
        rightList.Add(ParsePositiveNumber(curr));
    }
    
    ShellSort(leftList.arr, leftList.Size());
    ShellSort(rightList.arr, leftList.Size());

    int result = 0;
    for (int i = 0; i < leftList.Size(); i++) {
        int diff = Abs(leftList[i] - rightList[i]);
        result += diff;
    }

    printf("result is: %d", result);
    getchar();

    return 0;
}

static int day1_2()
{
    ScopedText data = ReadAllText("AdventInput1.txt");
    const char* curr = data;
    
    Array<int> leftList{};
    Map<int, int> rightList{};

    while (*curr)
    {
        leftList.Add(ParsePositiveNumber(curr));
        
        int left = ParsePositiveNumber(curr);
        KeyValuePair keyval = KeyValuePair<int, int>(left, 1);
     
        if (auto* node = rightList.FindNode(keyval))
        {
            node->value.value += 1;
        }
        else
        {
            rightList.Insert(keyval);
        }
    }

    QuickSort(leftList.arr, 0, leftList.Size()-1);

    int result = 0;
    for (int i = 0; i < leftList.Size(); )
    {
        int number = leftList[i];
        auto* node = rightList.FindNode(KeyValuePair<int, int>(number, 1));
        
        if (node != nullptr) {
            result += number * node->value.value;
        }

        while (leftList[i] == number && i < leftList.Size()) {
            i++;
        }
    }
    printf("result is: %d", result);
    getchar();

    return 0;
}

static bool IsSafe(int* numbers, int n)
{
	int prev = numbers[0];
	bool isIncreasing = true;
	bool isDecreasing = true;
	bool anyEqual = false;
    
	for (int i = 1; i < n; i++) 
	{
		int num = numbers[i];
		isIncreasing &= num >= prev;
		isDecreasing &= num <= prev;
		anyEqual |= num == prev;

		if (Abs(num - prev) > 3) {
			isIncreasing = false;
			isDecreasing = false;
		}
		prev = num;
	}
    
	return (isIncreasing || isDecreasing) && !anyEqual;
}

int Day2()
{
	FILE* input = fopen("AdventInput2.txt", "r");
	char line[128];
	int numbers[12];
	int numbers2[12];
	int n = 0;
	int numSafe = 0;

	while (fgets(line, sizeof(line), input))
	{
		const char* curr = line;
		while (*curr > '\n')
		{
			int num;
			if (IsNumber(*curr)) num = *curr - '0', curr++;
			if (IsNumber(*curr)) num = num * 10 + *curr - '0', curr++;
			curr += *curr == ' ';
			curr += *curr == '\r';
			numbers[n++] = num;
		}
		
		if (IsSafe(numbers, n))
		{
			numSafe++;
		}
		#ifdef PART2
		else
		{
			int n2 = 0;
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
					if (j != i) numbers2[n2++] = numbers[j];
				
				if (IsSafe(numbers2, n2)) {
					numSafe++;
					break;
				}
				n2 = 0;
			}
		}
		#endif
		n = 0;
	}
	printf("num safe: %d",  numSafe);
	fclose(input);
	getchar();
	return numSafe;
}

int Day3()
{
	ScopedText text = ReadAllText("AdventInput3.txt");
	const char* curr = text.text;
	long result = 0;
	bool enabled = true;

	while (*curr)
	{
		#if PART2
		if (StartsWith(curr, "don't()")) enabled = false;
		if (StartsWith(curr, "do()")) enabled = true;
		#endif
		
		if (enabled && StartsWith(curr, "mul("))
		{
			int a = -1, b = -1;

			if (IsNumber(*curr)) a = *curr++ - '0'; 
			if (IsNumber(*curr)) a = 10 * a + *curr++ - '0';
			if (IsNumber(*curr)) a = 10 * a + *curr++ - '0';
			
			if (a == -1 || *curr != ',')
				continue;
			curr++; // skip ','

			if (IsNumber(*curr)) b = *curr++ - '0';
			if (IsNumber(*curr)) b = 10 * b + *curr++ - '0';
			if (IsNumber(*curr)) b = 10 * b + *curr++ - '0';
		
			if (b == -1 || *curr != ')')
				continue;
			result += a * b;
		}

		curr += *curr != 0;
	}

	printf("result: %ld\n", result);
	getchar();
	return result;
}

int Day4()
{
    char lines[140][142];
    MemsetZero(lines, sizeof(lines));
    
    FILE* input = fopen("AdventInput4.txt", "r");
    
    int numLines = 0;
    while (fgets(lines[numLines], 142, input))
    {
        numLines++;
    }
    fclose(input);
    
    int appear = 0;
    int numCols = 142;

    // horizontal checks
    for (int i = 0; i < numLines; i++)
    {
        for (int x = 0; x < numCols - 3; x++)
        {
            appear += lines[i][x + 0] == 'X' && lines[i][x + 1] == 'M' && lines[i][x + 2] == 'A' && lines[i][x + 3] == 'S';
            appear += lines[i][x + 0] == 'S' && lines[i][x + 1] == 'A' && lines[i][x + 2] == 'M' && lines[i][x + 3] == 'X';
        }
    }

    // vertical checks
    for (int i = 0; i < numLines - 3; i++)
    {
        for (int x = 0; x < numCols; x++)
        {
            appear += lines[i + 0][x] == 'X' && lines[i + 1][x] == 'M' && lines[i + 2][x] == 'A' && lines[i + 3][x] == 'S';
            appear += lines[i + 0][x] == 'S' && lines[i + 1][x] == 'A' && lines[i + 2][x] == 'M' && lines[i + 3][x] == 'X';
        }
    }

    //  #
    //   #
    for (int i = 0; i < numLines-3; i++)
    {
        for (int x = 0; x < numCols-3; x++)
        {
            appear += lines[i+0][x+0] == 'X' && lines[i+1][x+1] == 'M' && lines[i+2][x+2] == 'A' && lines[i+3][x+3] == 'S';
            appear += lines[i+0][x+0] == 'S' && lines[i+1][x+1] == 'A' && lines[i+2][x+2] == 'M' && lines[i+3][x+3] == 'X';
        }
    }
    //   #
    //  #
    for (int i = 3; i < numLines; i++)
    {
        for (int x = 0; x < numCols - 3; x++)
        {
            appear += lines[i-0][x+0] == 'X' && lines[i-1][x+1] == 'M' && lines[i-2][x+2] == 'A' && lines[i-3][x+3] == 'S';
            appear += lines[i-0][x+0] == 'S' && lines[i-1][x+1] == 'A' && lines[i-2][x+2] == 'M' && lines[i-3][x+3] == 'X';
        }
    }
    printf("num apear: %d\n", appear);

    int numMas = 0;

    for (int i = 1; i < numLines-1; i++)
    {
        for (int x = 1; x < numCols-1; x++)
        {
            if (lines[i][x] == 'A')
            {
                bool a = (lines[i-1][x-1] == 'M' &&
                          lines[i+1][x+1] == 'S') |
                         (lines[i-1][x-1] == 'S' && 
                          lines[i+1][x+1] == 'M');
                
                bool b = a && ((lines[i-1][x+1] == 'S' &&
                                lines[i+1][x-1] == 'M') |
                               (lines[i-1][x+1] == 'M' &&
                                lines[i+1][x-1] == 'S'));
                numMas += a && b;
            }
        }
    }

    printf("num mas: %d", numMas);
    return appear;
}

/////////////     Day 5 missing         /////////////

typedef struct {
    int x, y;
} Vector2;

int OutOfBounds(Vector2 pos) {
    return pos.x < 0 || pos.x >= 131 || 
           pos.y < 0 || pos.y >= 130;
}

void Step(const char (&map)[130][132], Vector2& dir, Vector2& pos)
{
    if (map[pos.y + dir.y][pos.x + dir.x] == '#') 
    {
        // Rotate 90 degrees clockwise (-y, x)
        int temp = dir.x;
        dir.x = -dir.y;
        dir.y = temp;
    }
    else 
    {
        pos.x += dir.x;
        pos.y += dir.y;
    }
}

int Day6() {
    const int width = 132;
    const int height = 130;

    char map[height][width];
    Vector2 dir = {0, -1};
    Vector2 pos = {0, 0};

    MemsetZero(map, sizeof(map));
    FILE* input = fopen("AdventInput5.txt", "r");

    while (fgets(map[pos.y], width, input))
        pos.y++;

    fclose(input);
    
    Vector2 startPos = {0, 0};
    Vector2 startRot = dir;

    for (int y = 0; y < height && pos.x == 0; y++) 
        for (int x = 0; x < width; x++) 
            if (map[y][x] == '^') 
                pos.x = x, pos.y = y,
                startPos.x = x, startPos.y = y;
    // walk around with guardian
    while (!OutOfBounds(pos)) 
    {
        map[pos.y][pos.x] = ','; // mark as visited
        Step(map, dir, pos);
    }

    int numVisited = 0;
    for (int y = 0; y < height; y++) 
        for (int x = 0; x < width; x++) 
            numVisited += map[y][x] == ','; // count num visited

    printf("numVisited: %d\n", numVisited);

    // part 2
    int numInf = 0;
    for (int y = 0; y < height; y++) 
    {
        for (int x = 0; x < width; x++)
        {
            if (map[y][x] == '#' || (x == startPos.x && y == startPos.y))
                continue;

            map[y][x] = '#';
            pos = startPos;
            dir = startRot;

            int numTry = 0;
            while (!OutOfBounds(pos) && numTry < 8000)
            {
                Step(map, dir, pos);
                numTry++;
            }
            numInf += numTry >= 8000;
            map[y][x] = '.'; // revert
        }
    }
    
    printf("numInf: %d\n", numInf);
    return numVisited;
}

// Concatenate two numbers using the concatenation operator 
long long Concat(long long a, long long b) 
{
    long long power = 1;
    while (power <= b) {
        power *= 10;
    }
    return a * power + b;
}

// Recursive function to check if a target value can be achieved
bool Evaluate(long long numbers[16], int count, int index, long long current, long long target) 
{
    if (index == count) {
        return current == target;
    }
    // Add
    if (Evaluate(numbers, count, index + 1, current + numbers[index], target)) {
        return true;
    }

    // Multiply
    if (Evaluate(numbers, count, index + 1, current * numbers[index], target)) {
        return true;
    }
    #ifdef DAY2
    // Concatenate
    if (Evaluate(numbers, count, index + 1, Concat(current, numbers[index]), target)) {
        return true;
    }
    #endif
    return false;
}

int Day7()
{
    FILE *input = fopen("AdventInput7.txt", "r");
    if (!input) {
        perror("Failed to open file");
        return 1;
    }

    char line[128];
    long long totalCalibrationResult = 0;

    while (fgets(line, sizeof(line), input)) 
    {
        const char* curr = line;
        
        long long numbers[16];
        long long target = ParseNumberLong(curr);
        int count = 0;
        while (*curr > '\n')
        {
            numbers[count++] = ParseNumberLong(curr);
            while (*curr > '\n' && !IsNumber(*curr)) curr++;
        }
        curr++; // skip \n

        if (Evaluate(numbers, count, 1, numbers[0], target)) {
            totalCalibrationResult += target;
        }
    }

    fclose(input);
    printf("Total Calibration Result: %lld\n", totalCalibrationResult);
    return 0;
}

int main() 
{
    printf("Final Checksum: %lld\n", checksum);
    return 0;
}