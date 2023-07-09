# ASTL

Work in progress. This library is composition of lots of paper, research and development <br>
For compatibility I've used lots of macros in Common.hpp but still we have platform dependent code <br>
Also I've avoided most of the modern C++ features <br>
feel free to contribute, use, modify or oppening PR<br>
I havent used any C++ headers except thread and atomic. (by default threading is not compiled activate with ASTL_MULTI_THREADING)<br>
also I havent used any C headers except <stdint.h>. and in msvc <intrin.h> (for SIMD) <br>
Goal is make this library compile as fast as possible,  easy to use. and easy to read <br>
I have Hash functions and Random number generators in Random.hpp. <br>
Matrix4 and Vector4 uses SIMD extensions(SSE3)<br>
math library is combination of glm and DirectX Math. examples:
```cpp
#include "Math/Vector.hpp"
#include "Math/Matrix.hpp"

static float f = 1.0f; f += 0.01f;
constexpr float distance = 3.14159265f; // this is distance from cube but I did use pi anyways 
Vector3f position(Sin(f) * distance, 0.0f, Cos(f) * distance );
float verticalFOV = 65.0f, nearClip = 0.01f, farClip = 500.0f;

Matrix4 projection = Matrix4::PerspectiveFovRH(verticalFOV * DegToRad, m_NativeWindow->GetWidth(), m_NativeWindow->GetHeight(), nearClip, farClip);
Matrix4 viewProjection = projection * Matrix4::LookAtRH(position, -Vector3f::Normalize(position), Vector3f::Up());

Vector3f baryCentrics = Vector3f(1.0f - hitOut.u - hitOut.v, hitOut.u, hitOut.v);
Matrix3 inverseMat3 = Matrix4::ConvertToMatrix3(hitInstance.inverseTransform);
		
Vector3f n0 = Matrix3::Multiply(inverseMat3, ConvertToFloat3(&triangle.normal0x));
Vector3f n1 = Matrix3::Multiply(inverseMat3, ConvertToFloat3(&triangle.normal1x));
Vector3f n2 = Matrix3::Multiply(inverseMat3, ConvertToFloat3(&triangle.normal2x));
	        
record.normal = Vector3f::Normalize((n0 * baryCentrics.x) + (n1 * baryCentrics.y) + (n2 * baryCentrics.z));

Vector2i x = Vector2i(1, 1) + Vector2(2, 2);
float length =  ToVector2f(x).Length();
```
<br>
HashMap uses Ankerl's algorithm (way faster than std::unordered_map and uses contigues memory)<br>
I have RedBlackTree for Map and Set (ordered lookup tables) using comparison operator instead of hash<br>
ScopedPtr, and SharedPtr. coming soon.<br>
I have Array<T> instead of std::vector<T> <br>
also you can use Bitset1024, Bitset512, Bitset256(SIMD optimized) or Bitset<1234> instead of std::bitset <br>
Here is examples

```cpp
#include "Queue.hpp"
#include "Stack.hpp"
#include "String.hpp"
#include "Array.hpp"
#include "RedBlackTree.hpp"

// Map<uint64_t, String> mapTest{};
Set<uint64_t> set{};
    
Random::PCG pcg;
Random::PCGInitialize(pcg, 12345);

for (int i = 0; i <= 1000; i++)
  set.Insert(Random::PCGNext(pcg) & 1023); // 

Random::PCGInitialize(pcg, 12345);

printf("is contains: %i\n", set.Contains(552));
printf("root data: %llu\n", set.m_root->value);

set.Traverse(Travert);

bool containsAll = 1;
for (int i = 0; i <= 1000; i++)
    containsAll &= map.Contains(Random::PCGNext(pcg) & 1023);

ASSERT(containsAll);

Random::PCGInitialize(pcg, 12345);

for (int i = 0; i <= 1000; i++)
    set.Erase(Random::PCGNext(pcg) & 1023);

set.Traverse(Travert);

String testStr = "floating test: ";
testStr += 1234.567f;
printf("%s\n", testStr.c_str());
String testStr1 = "int test: ";
testStr1 += 1234567;
printf("%s\n", testStr1.c_str());
static const KeyValuePair<int, String> testMapInitializer[13] =
{
    { 23     , "apple"    },
    { 456    , "banana"   },
    { 789    , "cherry"   },
    { 42     , "dog"      },
    { 8675309, "elephant" },
    { 555    , "frog"     },
    { 999    , "guitar"   },
    { 314159 , "house"    },
    { 271828 , "igloo"    },
    { 777    , "jigsaw"   },
    { 888    , "kangaroo" },
    { 666    , "lemon"    },
    { 99999  , "monkey"   }
};
static HashMap<int, String> testMap(testMapInitializer, 13);
        
testMap[4] = "no its five";
String strings[13] = {
    "apple", "banana", "cherry", "dog", "elephant", "frog", "guitar",
    "house", "igloo", "jigsaw", "kangaroo", "lemon", "monkey"
};       

Array<String> strs{};

for (int i = 0; i < 13; ++i)
{
    strs.Add(strings[i]);
}

strs.RemoveAt(3, 3);
for (int i = 0; i < strs.Size(); ++i)
{
    printf("%s \n", strs[i].c_str());
}

for (const auto&[key, value] : testMap)
{
    printf("%i, %s\n", key, value.CStr());
}

HashSet<int> testSet{};
testSet.insert(213444590);
HashMap<int, String> ourMap = testMap;
TestMap("ourMap", ourMap);
String num = ourMap[4];

bool finded1 = ourMap.Find(5ull) != ourMap.end();
String val1 = ourMap.At(4);

ourMap.TryEmplace(33, "33ull");
bool contains1 = ourMap.Contains(33);

ourMap.Erase(33);
contains1 = ourMap.Contains(33ull);

ourMap.Insert(33, "33ull");
contains1 = ourMap.Contains(33ull);
```

here is the test of queue and stack

```cpp
Stack<uint64_t> stack{};
Array<uint64_t> stackArray;
uint64_t xoro[2];
Random::Xoroshiro128PlusInit(xoro);
char chr[9]{};

for (int i = 0; i < 10001; ++i)
{
    uint64_t x = Random::Xoroshiro128Plus(xoro);
    x &= ((1 << 15) - 1);
    stackArray.Add(x);
    stack.Push(x);
}

uint64_t sameCount = 0;

while (stack.Any())
{
    bool same = stack.Top() == stackArray.Back();
    sameCount += same;
    ASSERT(same);
    stack.Pop();
    stackArray.RemoveAt(stackArray.Size() - 1);
}
printf("stack is all same: %i\n", (int)sameCount);

printf("\n\n");

Queue<uint64_t> queue(33);
Array<uint64_t> queueArray(50);
for (int i = 0; i < 33; ++i)
    queue.Enqueue(i), queueArray.EmplaceBack(i);

queue.Enqueue(97); queueArray.EmplaceBack(97);
queue.Enqueue(98); queueArray.EmplaceBack(98);
queue.Enqueue(98); queueArray.EmplaceBack(98);
queue.Enqueue(98); queueArray.EmplaceBack(98);
queue.Enqueue(98); queueArray.EmplaceBack(98);

for (int i = 0; i < 10; i++)
    queue.Dequeue(), queueArray.RemoveAt(0);

Array<uint64_t>::Iterator ait = queueArray.begin();

while (!queue.IsEmpty())
   printf("%i\n", queue.Dequeue() == *ait++);
```

also you can use BVH algorithm(Jacco Bikker's but optimized with SIMD) for Ray casting to the 3d scene<br>
it will return mesh index, hit position hit color etc.<br>
note: I haven't touched multi threading code for a while