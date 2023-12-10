#include "../Queue.hpp"
#include "../Stack.hpp"
#include "../String.hpp"
#include "../Array.hpp"
#include "../HashMap.hpp"
#include "../HashSet.hpp"
#include "../RedBlackTree.hpp"
#include "../Math/Vector.hpp"

#include <stdio.h>

// #include <chrono>
// #include <thread>

AX_NAMESPACE

template<> struct Hasher<Vector2i>
{
    static __forceinline uint64_t Hash(Vector2i vec)
    {
        uint64_t x = (uint64(vec.y) << 32) | vec.x;
        x *= 0xbf58476d1ce4e5b9ULL;
        x ^= x >> 30ULL;
        return x;
    }
};

template<> struct Hasher<Vector2s>
{
	static __forceinline uint64_t Hash(Vector2s vec)
	{
		uint64_t x = (uint64(vec.y) << 16) | vec.x;
		x ^= x << 30ULL;
		x *= 0xbf58476d1ce4e5b9ULL;
		return x;
	}
};

#pragma warning(disable: 4996) // fopen
#pragma warning(disable: 4554) // | operator

template<typename T>
inline int ManhattanDistance(Vector2<T> a, const Vector2<T>& b)
{
    return Abs(a.x - b.x) + Abs(a.y - b.y);
}

static void Day15() // result 5511201
{
    FILE* file = fopen("Test/AOC15.txt", "r");
    char line[120];
    HashMap<Vector2i, int> sensors{};
    Set<int> beaconXs{};
    Vector2i boundsMIN = MakeVec2(INT32_MAX), boundsMAX = MakeVec2(INT32_MIN);

    while (fgets(line, sizeof(line), file))
    {
        const char* curr = line;

        Vector2i sensorPos;
        sensorPos.x = ParseNumber(curr);
        sensorPos.y = ParseNumber(curr);

        Vector2i beaconPos;
        beaconPos.x = ParseNumber(curr);
        beaconPos.y = ParseNumber(curr);

        Vector2i distance = MakeVec2(Abs(sensorPos.x - beaconPos.x), Abs(sensorPos.y - beaconPos.y)); // ManhattanDistance
        sensors[sensorPos] = distance.x + distance.y;
        beaconXs.Insert(beaconPos.y == 2000000 ? beaconPos.x : INT32_MIN);
        boundsMIN = Min(boundsMIN, beaconPos - distance);
        boundsMAX = Max(boundsMAX, beaconPos + distance);
    }

	int result = 0;
    // check each column if it contains # or not
    for (int j = boundsMIN.x; j <= boundsMAX.x; ++j)
    {
        Vector2i columnPos = MakeVec2(j, 2000000);
        if (sensors.Contains(columnPos)
            || beaconXs.Contains(columnPos.x)) continue; // if this beacon is sensor we will not count this

        // for (auto const& [pos, dist] : sensors)
        for (const KeyValuePair<Vector2i, int>& x : sensors)
        {
            const Vector2i& pos = x.key;
            int dist = x.value;

            int columnToSensor = ManhattanDistance(pos, columnPos);
            if (columnToSensor > 0 && columnToSensor <= dist) { result++; break; }
        }
    }
    printf("Day15 result: %i\n ", result); // result should be 5511201
	ASSERT(result == 5511201);
    fclose(file);
}

struct APoint
{
	union { float costSoFar; float distance; };
	union { Vector2s cameFrom; Vector2s point; };

	APoint() : costSoFar(1e30f) { cameFrom.x = 0; cameFrom.y = 0; }
	APoint(float costSoFar, Vector2s pnt) : costSoFar(costSoFar), cameFrom(pnt) {}

	bool operator < (const APoint& o) const { return distance < o.distance; }
	bool operator > (const APoint& o) const { return distance > o.distance; }
	bool operator == (const APoint& o) const { return point == o.point; }
	bool operator != (const APoint& o) const { return point != o.point; }
};

using DistanceAndPoint = APoint;

static int Day12() // result should be 534
{
    FILE* file = fopen("Test/AOC12.txt", "r");
	TimeBlock("day12");
	Vector2s startPos, targetPos;
	short numRows = 0, numColumns = 0;
	char grid[42][167]; // for now constant size, std::string may used for bigger grids

	// parse input, find num rows, num columns, start and end point 
	while (fgets(grid[numRows], 165, file))
	{
		char* curr = grid[numRows]; numColumns = 0;

		if (*curr == 'S') { startPos.x = numColumns;  startPos.y = numRows; }

		while (*curr > '\n' && !IsWhitespace(*curr))
		{
			if (*curr == 'E') { targetPos.x = numColumns; targetPos.y = numRows; }
			curr++, numColumns++;
		}
		numRows++;
	}
	fclose(file);

	grid[startPos.y][startPos.x] = 'a';
	grid[targetPos.y][targetPos.x] = 'z';

	// implement A* Search Algorithm
	HashMap<Vector2s, APoint> map{64};
	// std::priority_queue<DistanceAndPoint, std::vector<DistanceAndPoint>, std::greater<DistanceAndPoint>> queue{};
	PriorityQueue<DistanceAndPoint, PQCompare_Greater> queue{512};

	queue.Emplace(0.0f, startPos);
	map[startPos] = APoint(0.0f, startPos);

	while (!queue.Empty())
	{
		Vector2s currentPoint = queue.Top().point;
		char height = grid[currentPoint.y][currentPoint.x];
		const auto find = map.ConstFind(currentPoint);
		if (currentPoint == targetPos || find == map.cend()) break;

		APoint current = find->value;

		queue.Pop();

		auto const processNeighbor = [&](Vector2s point)
		{
			if (point.x < 0 || point.y < 0 || point.x > numColumns || point.y > numRows
				|| grid[point.y][point.x] > height + 1) return;

			float newCost = current.costSoFar + 1.0f;
			APoint& neighborPoint = map[point];

			if (neighborPoint.costSoFar == 1e30f || newCost < neighborPoint.costSoFar)
			{
				neighborPoint.costSoFar = newCost;
				neighborPoint.cameFrom = currentPoint;
				float priority = newCost + Vector2s::DistanceSq(point, targetPos);
				queue.Emplace(priority, point);
			}
		};

		processNeighbor(currentPoint + MakeVec2<short>(1, 0));
		processNeighbor(currentPoint + MakeVec2<short>(0, 1));
		processNeighbor(currentPoint + MakeVec2<short>(0, -1));
		processNeighbor(currentPoint + MakeVec2<short>(-1, 0));
	}

	Vector2s currentPoint = targetPos;
	int result = 0;
	while (currentPoint != startPos)
	{
		result++;
		grid[currentPoint.y][currentPoint.x] = '#';
		currentPoint = map[currentPoint].cameFrom;
	}
	printf("Day 12 min steps: %d\n", result);
	ASSERT(result == 534);
	return result;
}

static int Day17() // result should be 3157
{
    char* pattern = ReadAllFile("Test/AOC17.txt");
	TimeBlock("Day17")

	const Vector2s shapes[5][5] = {
		{ {0, 0}, {1, 0}, {2, 0}, {3, 0} }, // ####
		{ {1, 0}, {0, 1}, {1, 1}, {2, 1}, {1, 2} }, // + shape
		{ {0, 0}, {1, 0}, {2, 0}, {2, 1}, {2, 2} }, // L shape
		{ {0, 0}, {0, 1}, {0, 2}, {0, 3} }, // | shape
		{ {0, 0}, {1, 0}, {0, 1}, {1, 1} } // box shape
	};
	const int shapeSizes[5] = { 4,5,5,4,4 }; // num pixels

	Vector2s mapBounds = MakeVec2<short>(7, -1);
	HashSet<Vector2s> blocks;

	auto const CheckColission = [&](const Vector2s* shape, int shapeSize, Vector2s position) 
	{
		bool collided = false;
		for (int i = 0; i < shapeSize; ++i)
		{
			Vector2s pixel = shape[i] + position;
			collided |= blocks.Contains(pixel) | pixel.x >= 7 | pixel.x < 0 | pixel.y < 0;
		}
		return collided;
	};

	const char* curr = pattern;
	int currBlock = 0, cleanHeight = 200;

	for (int g = 0; g < 2022; ++g)
	{
		const int blockIndex = currBlock++ % 5;
		const int shapeSize = shapeSizes[blockIndex];
		const Vector2s* shape = shapes[blockIndex];
		Vector2s position = MakeVec2<short>(2, mapBounds.y + 4);

		auto const visualize = [&]()
		{
			for (int i = 0; i < shapeSize; ++i) 
				blocks.Insert(position + shape[i]);
		
			// system("cls");
			for (int i = MAX(mapBounds.y, (short)7)+4; i >= 0; --i)
			{
				for (int j = 0; j < mapBounds.x; ++j)
				{
					if (blocks.Contains(MakeVec2<short>(j, i))) printf("#");
					else printf(".");
				}
				printf("\n");
			}
		
			for (int i = 0; i < shapeSize; ++i) 
				blocks.Erase(position + shape[i]);
			// using namespace std::chrono_literals;
			// std::this_thread::sleep_for(300ms);
		};

		//visualize();
		while (true)
		{
			short direction = *curr++ - '='; // with ascii table it is clear to understand: '<','=','>'  returns -1 or 1 depending on input axis 
			position.x += direction;
			if (CheckColission(shape, shapeSize, position)) position.x -= direction;
			// 	visualize();
			position.y -= 1; // gravity
			if (CheckColission(shape, shapeSize, position)) { position.y++; break; }
			// visualize();
			if (*curr <= '\n')  curr = pattern; // if we are in the end of the pattern repeat pattern
		}

		for (int i = 0; i < shapeSizes[blockIndex]; ++i)
		{
			Vector2s pixel = shape[i] + position;
			blocks.Insert(pixel); // insert pixels
			mapBounds.y = MAX(mapBounds.y, pixel.y); // check new pixel is higher than top point
		}
		// we don't want to store all of the data that's why we are removing unnecesarry blocks that too below us
		if (mapBounds.y > cleanHeight + 100)
		{
            blocks.EraseIf([cleanHeight](Vector2s v) { return v.y < cleanHeight; });
			cleanHeight += 100;
		}
	}

	// writing to a file is aditional, not must
	// FILE* file = fopen("17Result.txt", "w");	
	// for (short i = MAX<short>(mapBounds.y, 7) + 4; i >= 0; --i)
	// {
	// 	char line[8];
	// 	for (short j = 0; j < 7; ++j)
	// 	{
	// 		if (blocks.Contains(Vector2s(j, i))) line[j] = '#';
	// 		else line[j] = '.';
	// 	}
	// 	line[7] = '\n';
    //     fwrite(line, 1, 8, file);
	// }
	// fclose(file);
	free(pattern);
	printf("Day17 result is: %d\n", mapBounds.y+1); // result should be 3157
	ASSERT(mapBounds.y + 1 == 3157);
	return mapBounds.y;
}


int Day22()
{
	FILE* file = fopen("Test/AOC22.txt", "r");
	char line[5666];
	
	TimeFunction

	HashMap<Vector2i, char> map;
	Vector2i mapBounds{ 0, 0 };
	// parse map
	while (fgets(line, sizeof(line), file))
	{
		if (IsNumber(*line)) break;
		const char* curr = line;
		int oldMAXX = mapBounds.x;
		mapBounds.x = 0;

		while (*curr > '\n')
		{
			if (!IsWhitespace(*curr)) map[mapBounds] = *curr;
			mapBounds.x++, curr++;
		}
		mapBounds.x = MAX(mapBounds.x, oldMAXX);
		mapBounds.y++;
	}
	const char* path = line; // last line of text is our path

	const Vector2i directions[4] = {
		{ 0,-1},
		{ 1, 0}, // we start from here, and depending on the input we go down or up, if we out of bounds we loop this array x % 4
		{ 0, 1},
		{-1, 0}
	};

	int currentDirection = 1;
	Vector2i currentPosition{ 0, 0 };

	auto const FindLoopPosition = [&](Vector2i position, const Vector2i& inverseDir)
	{
		while (map.Contains(position + inverseDir))
		{
			position += inverseDir;
		}
		return map[position] == '#' ? -Vector2i::One() : position;
	};

	auto const Visualize = [&]()
	{
		// system("cls"); // clear screen
		for (int y = MAX(currentPosition.y - 30, 0); y < MIN(currentPosition.y + 30, mapBounds.y); ++y)
		{
			for (int x = 0; x < mapBounds.x; ++x)
			{
				Vector2i pos = MakeVec2(x, y);
				auto const find = map.ConstFind(pos);
				if (pos == currentPosition) printf("X");
				else if (find == map.cend()) printf(" ");
				else printf("%c", find->value);
			}
			printf("\n");
		}
		// using namespace std::chrono_literals;
		// std::this_thread::sleep_for(500ms);
	};

	// find startPosition
	while (!map.Contains(currentPosition)) currentPosition.x++;

	while (*path > '\n')
	{
		int steps = ParseNumber(path);

		while (steps--)
		{
			Vector2i direction = directions[currentDirection];
			Vector2i newPos = currentPosition + direction;
			const auto find = map.Find(newPos);

			if (find == map.end())
			{
				// find loop map position by going inverse of our direction
				Vector2i loopPosition = FindLoopPosition(currentPosition, -direction);
				if (loopPosition == -Vector2i::One()) break;
				else { currentPosition = loopPosition; continue; }
			}

			if (find->value == '#') break;

			currentPosition = newPos;
			// Visualize();
		}

		int direction = *path++;
		if (direction == 'L') currentDirection = currentDirection == 0 ? 3 : currentDirection - 1;
		if (direction == 'R') currentDirection = currentDirection == 3 ? 0 : currentDirection + 1;
	}
	// rotate left for getting index
	currentDirection = currentDirection == 0 ? 3 : currentDirection - 1; // our indices start from 1 to 3

	currentPosition += Vector2i::One(); // result = one indexed array
	int result = 1000 * currentPosition.y + (4 * currentPosition.x) + currentDirection;

	printf("day22 result: %d\n", result); 
	ASSERT(result == 58248);
	return 0;
}

AX_END_NAMESPACE

void AdventOfCodeTests22()
{
	Day15();
	Day12();
	Day17();
	Day22();
}

