
# ASTL standard library

Astl is a data structures, algorithms and math library that is targeting high performance<br>
faster compile times and ease of use, the code is writen in a way that is easier to read<br>
modify and use, there is no other external dependencies C++ headers are not included<br>
not even C headers included, so everything is from scratch. this should work with C++14 and above<br>

only include is immintrin.h that is for simd code,<br>
speaking about simd, math library and other functions has optional simd instructions<br>
to better use cpu hardware<br>

I've been working on this project almost one year, inspired from Stewart Lynch's Sane C++ videos,<br>
although his data structures are not publicly awailable, I've made myself <br>
and improved my algorithm skills a lot<br>

I've tried to use templates as less as possible, but library has templates but with sane way<br>

what pushed me into this is compile times of standard library and its slow performance <br>
with debug mode<br>
## Appendix
All data structures has custom allocator support, but RedBlackTree has Growable fixed Size allocator that Stewart Lynch mentioned one of his videos,<br>
it is fast to allocate and RedBlackTree is fast because of this.
here is the data structures that I have:

* Array
* String
* Stack
* Queue
* PriorityQueue
* HashMap  // using Ankerl's algorithm (way faster than stl)
* HashSet
* Map // using Red Black Tree
* Set

each data structure has an header with same name except Map and Set these structures uses "RedBlackTree.hpp" header.

here is the example usage of data structures
```cpp
#include <ASTL/RedBlackTree.hpp>
#include <ASTL/HashMap.hpp>
#include <ASTL/Math/Vector.hpp>

FILE* file = fopen("TestData/AOC15.txt", "r");
Array<char> line(120, '\n');
HashMap<Vector2i, int> sensors{};
Set<int> beaconXs{};
Vector2i boundsMin = INT32_MAX, boundsMax = INT32_MIN;

while (fgets(line, sizeof(line), file))
{
    const char* curr = line;

    Vector2i sensorPos = ParseVector<Vector2i>(curr);
    Vector2i beaconPos = ParseVector<Vector2i>(curr);

    Vector2i distance = Vector2i(Abs(sensorPos.x - beaconPos.x), Abs(sensorPos.y - beaconPos.y)); // ManhattanDistance
    sensors[sensorPos] = distance.x + distance.y;
    beaconXs.Insert(beaconPos.y == 2'000'000 ? beaconPos.x : INT32_MIN);
    boundsMin = Min(boundsMin, beaconPos - distance);
    boundsMax = Max(boundsMax, beaconPos + distance);
}

int result = 0;
// check each column if it contains # or not
for (int j = boundsMin.x; j <= boundsMax.x; ++j)
{
    Vector2i columnPos = Vector2i(j, 2'000'000);
    // if this beacon is sensor we will not count this
    if (sensors.Contains(columnPos) || beaconXs.Contains(columnPos.x)) continue;

    for (const KeyValuePair<Vector2i, int>& x : sensors)
    {
        const Vector2i& pos = x.key;
        int dist = x.value;

        int columnToSensor = ManhattanDistance(pos, columnPos);
        if (columnToSensor > 0 && columnToSensor <= dist) { result++; break; }
    }
}
printf("Day15 result: %i\n ", result);
fclose(file);
```
data structures and algorithms are tested with many algorithms<br>
you can see more examples in AdventOfCodeTests.cpp<br>


# Math

The Math component of the ASTL library combines elements from glm and XNA Math<br>
while providing a simpler code structure that is easier to read, modify, and compile. <br>
The library offers various math structures to support common mathematical operations. Here are the available structures:<br>

* Vector2f, Vector2d, Vector2i, ...
* Vector3f, Vector3i, ...
* Vector4, ...
* Matrix4 (4x4 matrix)
* Matrix3 (3x3 matrix)
* Quaternion
* Transform
* Camera
The vectors are defined in the "Vector.hpp" header, while Matrix4 and Matrix3 are defined in the "Matrix.hpp" header.<br>
Other math structures have headers with the same name for easy access and organization.<br>

One notable feature of the ASTL math library is its simplicity and compactness.<br>
For example, the total lines of code for all vectors, including SIMD and non-SIMD versions,<br>
amount to only 440 lines. The Stack data structure consists of just 160 lines, and the QueueHeader totals 500 lines, including PriorityQueue.<br>
The Transform class encapsulates position, rotation, and scale properties, and it can be used with Euler angles as well.<br>
The Math.hpp header file contains all the essential math functions, including:<br>

* Pow
* Exp
* Sqrt
* Fract
* Log
* Log2
* Log10
* Floor
* Sign
* CopySign
* Abs
* FAbs
* Min
* Max

Additionally, the library provides trigonometric functions such as:

* Sin
* Cos
* Tan
* ATan
* Atan2
* ASin
* ACos
<br>
Most of these functions are constexpr, and for those that are not,<br>
links are provided so that you can modify them to be constexpr as well.<br>
Although constexpr math functions are officially introduced in C++23, you can still use these functions with C++14.<br>
Example Usage of the Math Library:<br>
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

The provided code demonstrates the usage of the math library.<br>
It showcases various operations involving vectors, matrices, and other math functions.<br>
Feel free to explore and utilize the math library based on your project's requirements.

here how it looks internaly, with bunch of links and comments that'll help you
```cpp
// https://mazzo.li/posts/vectorized-atan2.html
FINLINE constexpr float ATan(float x) {
	const float x_sq = x * x;
	constexpr float a1 =  0.99997726f, a3 = -0.33262347f, a5  = 0.19354346f,
	                a7 = -0.11643287f, a9 =  0.05265332f, a11 = -0.01172120f;
	return x * (a1 + x_sq * (a3 + x_sq * (a5 + x_sq * (a7 + x_sq * (a9 + x_sq * a11)))));
}

// https://en.wikipedia.org/wiki/Sine_and_cosine
// warning: accepts input between -TwoPi and TwoPi  if (Abs(x) > TwoPi) use x = FMod(x + PI, TwoPI) - PI;
FINLINE constexpr float Sin(float x) 
{
	float xx = x * x * x;                // x^3
	float t = x - (xx * 0.16666666666f); // x3/!3  6 = !3 = 1.6666 = rcp(3)
	t += (xx *= x * x) * 0.00833333333f; // x5/!5  120 = !5 = 0.0083 = rcp(5)
	t -= (xx *= x * x) * 0.00019841269f; // x7/!7  5040 = !7
	t += (xx * x * x)  * 2.75573e-06f;   // 362880 = !9
	return t;
}
```

Math library also has half to float, float to half conversion functions
and color packing and unpacking.

# Random

The Random component of the ASTL library provides functionalities for generating random numbers and hash functions. It includes various algorithms and utilities to handle randomness effectively. Here are the key features of the Random.hpp header:

Random Number Generation
Improved hash functions: The library incorporates improved hash functions from the article https://nullprogram.com/blog/2018/07/31/, including MurmurHash, WangHash, and their inverses.
Random seeding: The Random namespace provides functions for seeding the random number generator, such as Random::Seed32 and Random::Seed64.
32-bit random number generation: The Random::PCG structure and related functions, such as Random::PCGInitialize and Random::PCGNext, enable the generation of 32-bit random numbers.
64-bit random number generation: The Random::Xoroshiro128Plus algorithm and the accompanying Random::Xoroshiro128PlusInit function facilitate 64-bit random number generation.

Usage Examples:

```cpp
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

using namespace Random;
uint64_t xoro[2];
Xoroshiro128PlusInit(xoro);
uint64_t* ip = (uint64_t*)p;

for (uint64_t i = 0; i < len / sizeof(uint64_t); i++)
{
    *ip++ = Xoroshiro128Plus(xoro);
}
```

# Memory
As I said everything is from scratch, I've also made: MemCpy, MemSet functions
these are optimized, and uses template arguments alignment and size to determine
alignment at runtime, and to preform single instruction with __movsb you should set size argument otherwise code will determinate alignment at runtime

Also library has Move<T>, and Forward<T> functions instead of std::move and std::forward

this tested code performs same with simd and without simd, and compared to memcpy this is almost same, sometimes faster
```cpp
inline void MemCpyAligned32(void* dst, const void* src, uint64_t sizeInBytes)
{
    uint32_t*       dp  = (uint32_t*)dst;
    const uint32_t* sp  = (const uint32_t*)src;
    const uint32_t* end = (const uint32_t*)((char*)src + (sizeInBytes >> 2));
        
    while (sp < end)
    {
        dp[0] = sp[0];
        dp[1] = sp[1];
        dp[2] = sp[2];
        dp[3] = sp[3];
        dp += 4, sp += 4;
    }
    
    switch (sizeInBytes & 3)
    {
        case 3: *dp++ = *sp++;
        case 2: *dp++ = *sp++;
        case 1: *dp++ = *sp++;
    };
}

// use size for struct/class types such as Vector3 and Matrix4,
// and use MemCpy for big arrays or unknown size arrays
template<int alignment = 0, int size = 0>
inline void MemCpy(void* dst, const void* src, uint64_t sizeInBytes)
{
    if constexpr (size != 0)
#ifdef _MSC_VER
    __movsb((unsigned char*)dst, (unsigned char const*)src, size);
#else
    __builtin_memcpy(dst, src, size);
#endif
    else if (alignment == 8)
        MemCpyAligned64(dst, src, sizeInBytes);
    else if (alignment == 4)
        MemCpyAligned32(dst, src, sizeInBytes);
    else
    {
        // check if both has same alignment
        const uint64_t alignXor = uint64_t(dst) ^ uint64_t(src); 
        if (!(alignXor & 7))
            MemCpyAligned64(dst, src, sizeInBytes);
        else if (!(alignXor & 3))
            MemCpyAligned32(dst, src, sizeInBytes);
        else
        {
            const char* cend = (char*)((char*)src + sizeInBytes);
            const char* scp = (const char*)src;
            char* dcp = (char*)dst;
        
            while (scp < cend) *dcp++ = *scp++;
        }
    }
}

```
unfortunately I've used templates with, allocators but using alignment and size instead of
templates will have less code generated and it will be faster to compile so it is in my todo list.
```cpp
template<typename T> 
struct Allocator 
{
    static constexpr bool IsPOD = false;

    T* Allocate(int count) const; 
    T* AllocateUninitialized(int count) const;
    void Deallocate(T* ptr, int count) const;
    T* Reallocate(T* ptr, int oldCount, int count) const;
};

template<typename T>
struct MallocAllocator
{
    static constexpr bool IsPOD = true;

    T* Allocate(int count) const; 
    T* AllocateUninitialized(int count) const;
    void Deallocate(T* ptr, int count) const;
    T* Reallocate(T* ptr, int oldCount, int count) const;
};

template<typename T>
struct FixedSizeGrowableAllocator // for red black tree
{
    static constexpr bool IsPOD = false;
    static constexpr int InitialSize = NextPowerOf2(512 / Min((int)sizeof(T), 128));

    struct Fragment
    {
        Fragment* next;
        T*        ptr ;
        int64_t   size; // used like a index until we fill the fragment
    };

    int currentCapacity   = 0;
    Fragment* base    = nullptr;
    Fragment* current = nullptr;
    ....
};
```

# Algorithms

"Algorithms.hpp" Header is simple, only 280 line open to sugesstions and improvements.
it includes number parsing, sorting it does the job of <algorithm> header 

Also Bitset.hpp for fast bitsets (SIMD optimized). Bitset1024, Bitset512, Bitset256 or Bitset<1234> instead of std::bitset <br>

```cpp
struct Bitset256
{
	union {
		unsigned long bits[4] = { 0 };
		__m256i sse;
	};
	Bitset256() { Clear(); }
	Bitset256(const char* str) {
		int currBits = 0, i = 0; Clear();  while (*str) bits[currBits] |= (1 << (i & 63)) * (*str++ - '0'), currBits = ++i >> 6;
	}
	Bitset256(ulong r) { bits[0] = r; bits[1] = r; bits[2] = r; bits[3] = r; }
	Bitset256(ulong a, ulong b, ulong c, ulong d) { bits[0] = a; bits[1] = b; bits[2] = c; bits[3] = d; }
	Bitset256(__m256i v) : sse(v) {}
	// to set first n bits use this technique with simd(no need blend):
	// https://gist.github.com/benanil/78ad3600f5e10b9a3f6173afc8565352#file-memcpymemset-cpp-L101
	bool operator[](int idx) const { return !!(bits[idx >> 6] & (1ul << (idx & 63ul))); }
	bool Get(int idx) const { return !!(bits[idx >> 6] & (1ul << (idx & 63ul))); }
	void Set(int idx) { bits[idx >> 6] |= 1ul << (idx & 63); }
	void Reset(int idx) { bits[idx >> 6] &= ~(1ul << (idx & 63)); }

	void Clear() { sse = _mm256_set_epi64x(0ll, 0ll, 0ll, 0ll);  }
	void Flip() { sse = _mm256_xor_si256(sse, _mm256_set1_epi32(0xffffffff)); }
	Bitset256 operator~ () { return { _mm256_xor_si256(sse, _mm256_set1_epi32(0xffffffff)) }; }

	Bitset256 operator &  (const Bitset256 other) const { return _mm256_and_si256(sse, other.sse); }
	Bitset256 operator |  (const Bitset256 other) const { return _mm256_or_si256 (sse, other.sse); }
	Bitset256 operator ^  (const Bitset256 other) const { return _mm256_xor_si256(sse, other.sse); }
	Bitset256 operator &= (const Bitset256 other) { sse = _mm256_and_si256(sse, other.sse); return *this; }
	Bitset256 operator |= (const Bitset256 other) { sse = _mm256_or_si256(sse, other.sse);  return *this; }
	Bitset256 operator ^= (const Bitset256 other) { sse = _mm256_xor_si256(sse, other.sse); return *this; }

	bool All() const {
		return _mm256_movemask_epi8(_mm256_cmpeq_epi64(sse, _mm256_set_epi64x(~0ul, ~0ul, ~0ul, ~0ul))) == ~0ul;
	}
	bool Any() const {
		return _mm256_movemask_epi8(_mm256_cmpgt_epi64(sse, _mm256_setzero_si256())) > 0;
	}
	
	int64 Count() const { return popcount256_epi64(sse); }
};

```

To detect the runtime hardware SIMD (Single Instruction, Multiple Data) support, the code you provided can be formatted as follows:

```cpp
enum CPUIDBits : int
{
    CPUIDBits_SSE    = (1 << 25),
    CPUIDBits_SSE2   = (1 << 26),
    CPUIDBits_SSE3   = (1 << 9),
    CPUIDBits_SSE4_1 = (1 << 19),
    CPUIDBits_SSE4_2 = (1 << 20),
    CPUIDBits_AVX2   = (1 << 5),
    CPUIDBits_AVX512 = (1 << 16),
};

// sometimes you might need runtime extension detection.
// recommended using global variable like this in a cpp file:
// int g_ax_simd_bits = 0; or inline int g_ax_simd_bits = 0; in here.
// then call AX_InitSIMD_CPUID only once in program lifetime
// then use !!(g_ax_simd_bits & CPUIDBits_SSE2) to get the support value
inline int AX_InitSIMD_CPUID()
{
    int info[4];
    AX_CPUID(1, info);
    int mask = 1;
    mask |= info[3] & (CPUIDBits_SSE | CPUIDBits_SSE2);
    mask |= info[2] & (CPUIDBits_SSE3 | CPUIDBits_SSE4_1 | CPUIDBits_SSE4_2);
    AX_CPUID(7, info);
    mask |= info[1] & (CPUIDBits_AVX2 | CPUIDBits_AVX512);
    return mask;
}
```

also you can use BVH algorithm(Jacco Bikker's but optimized with SIMD) for Ray casting to the 3d scene<br>
it will return mesh index, hit position hit color etc. algorithm located in BVH folder<br>
note: I haven't touched multi threading code for a while<br>
<br>
Overall this is an hidden game for C++ and awesome to use with games and performance required programs.<br> 

# Contributing

feel free to contribute

#Contact
If you have any questions, feedback, or suggestions, feel free to reach out:<br>

Email: anilcangulkaya7@gmail.com<br>
Twitter: @anilcanglk12<br>
GitHub: @benanil<br>
<br>
Feel free to reach out regarding any inquiries related to the project.<br>
## Related

Here are some related and helpful links
 
 - [SIMD Sorting by (Wojciech Mula)](https://github.com/WojciechMula/simd-sort)
 - [Optimizing Binary Search - Sergey Slotin - CppCon 2022](https://www.youtube.com/watch?v=1RIPMQQRBWk)
 - [Better Assert by (Chris Wellons)](https://nullprogram.com/blog/2022/06/26/)
 - [Fast dense HashMap-Set by (Martinus Ankerl)](https://github.com/martinus/unordered_dense)
 - [Fast Matrix Inverse](https://lxjk.github.io/2017/09/03/Fast-4x4-Matrix-Inverse-with-SSE-SIMD-Explained.html)