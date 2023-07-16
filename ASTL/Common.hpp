#pragma once

#include "IntFltTypesLimits.hpp"

#if AX_SHARED
#ifdef AX_EXPORT
		#define AX_API __declspec(dllexport)
	#else
		#define AX_API __declspec(dllimport)
	#endif
#else
	#define AX_API
#endif

#ifndef FINLINE
#	ifdef _MSC_VER
#		define FINLINE __forceinline
#define VC_EXTRALEAN 1
#   elif __CLANG__
#       define FINLINE [[clang::always_inline]] 
#	elif __GNUC__
#       define FINLINE  __attribute__((always_inline))
#   endif
#endif

#ifndef VECTORCALL
#   ifdef _MSC_VER
#include <intrin.h>
#		define VECTORCALL __vectorcall
#   elif __CLANG__
#       define VECTORCALL [[clang::vectorcall]] 
#	elif __GNUC__
#       define VECTORCALL  
#   endif
#endif

#if defined(__GNUC__)
#    define AX_PACK(decl) decl __attribute__((__packed__))
#elif defined(_MSC_VER)
#    define AX_PACK(decl) __pragma(pack(push, 1)) decl __pragma(pack(pop))
#else
#error you should define pack function
#endif

#if defined(__GNUC__) || defined(__MINGW32__)
    #define RESTRICT __restrict__
#elif defined(_MSC_VER)
    #define RESTRICT __restrict
#else
    #define RESTRICT
#endif

#if defined( __GNUC__ ) || defined(__INTEGRITY)
#   define AX_ALIGNED(_x)          __attribute__ ((aligned(_x)))
#elif defined( _WIN32) && (_MSC_VER)                                                                                   
#	define AX_ALIGNED(_x)          __declspec(align(_x))      
#else
#   warning  Need to implement some method to align data here
#	define  CL_ALIGNED(_x)
#endif

#ifndef AXGLOBALCONST
#	if _MSC_VER
#		define AXGLOBALCONST extern const __declspec(selectany)
#	elif defined(__GNUC__) && !defined(__MINGW32__)
#		define AXGLOBALCONST extern const __attribute__((weak))
#   else 
#       define AXGLOBALCONST 
#	endif
#endif

#if defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__clang__)
#    define AX_LIKELY(x) __builtin_expect(x, 1)  
#    define AX_UNLIKELY(x) __builtin_expect(x, 0)
#else
#    define AX_LIKELY(x) (x)   
#    define AX_UNLIKELY(x) (x) 
#endif

// https://nullprogram.com/blog/2022/06/26/
#ifdef _DEBUG
#  if __GNUC__
#    define ASSERT(c) if (!(c)) __builtin_trap()
#  elif _MSC_VER
#    define ASSERT(c) if (!(c)) __debugbreak()
#  else
#    define ASSERT(c) if (!(c)) *(volatile int *)0 = 0
#  endif
#else
#  define ASSERT(c)
#endif

#if defined(__has_builtin)
#   define AX_COMPILER_HAS_BUILTIN(x) __has_builtin(x)
#else
#   define AX_COMPILER_HAS_BUILTIN(x) 0
#endif

#if AX_COMPILER_HAS_BUILTIN(__builtin_assume)
#   define AX_ASSUME(x) __builtin_assume(x)
#elif defined(_MSC_VER)
#   define AX_ASSUME(x) __assume(x)
#else
#   define AX_ASSUME(x) (void)(x)
#endif

#if AX_COMPILER_HAS_BUILTIN(__builtin_unreachable)
#   define AX_UNREACHABLE() __builtin_unreachable()
#elif _MSC_VER
#   define AX_UNREACHABLE() __assume(0)
#else
#   define AX_UNREACHABLE() 
#endif

#if AX_COMPILER_HAS_BUILTIN(__builtin_prefetch)
#   define AX_PREFETCH(x) __builtin_prefetch(x)
#elif defined(_MSC_VER)
#   define AX_PREFETCH(x) _mm_prefetch(x, _MM_HINT_T0)
#else
#   define AX_PREFETCH(x)
#endif

// https://gist.github.com/boxmein/7d8e5fae7febafc5851e
// https://en.wikipedia.org/wiki/CPUID
// example usage:
// void get_cpu_model(char *cpu_model) { // return example: "AMD Ryzen R5 1600"
//     int* cpumdl = (int*)cpu_model;
//     AX_CPUID(0x80000002, cpumdl); cpumdl += 4;
//     AX_CPUID(0x80000003, cpumdl); cpumdl += 4;
//     AX_CPUID(0x80000004, cpumdl); 
// }
// int arr[4];
// AX_CPUID(1, arr);
// int numCores = (arr[1] >> 16) & 0xff; // virtual cores included
#if defined(__clang__) || defined(__GNUC__)
#   include <cpuid.h>
#   define AX_CPUID(num, regs) __cpuid(num, regs[0], regs[1], regs[2], regs[3])
#else
#   define AX_CPUID(num, regs) __cpuid(regs, num)
#endif

/* Architecture Detection */
// detection code from mini audio
// you can define AX_NO_SSE2 or AX_NO_AVX2 in order to disable this extensions
#if defined(__x86_64__) || defined(_M_X64)
#   define AX_X64
#elif defined(__i386) || defined(_M_IX86)
#   define AX_X86
#elif defined(__arm__) || defined(_M_ARM) || defined(__arm64) || defined(__arm64__) || defined(__aarch64__) || defined(_M_ARM64)
#   define AX_ARM
#endif

/* Intrinsics Support */
#if (defined(AX_X64) || defined(AX_X86)) && !defined(__COSMOPOLITAN__)
    #if defined(_MSC_VER) && !defined(__clang__)
        #if _MSC_VER >= 1400 && !defined(AX_NO_SSE2)   /* 2005 */
            #define AX_SUPPORT_SSE
        #endif
        #if _MSC_VER >= 1700 && !defined(AX_NO_AVX2)   /* 2012 */
            #define AX_SUPPORT_AVX2
        #endif
    #else
        #if defined(__SSE2__) && !defined(AX_NO_SSE2)
            #define AX_SUPPORT_SSE
        #endif
        #if defined(__AVX2__) && !defined(AX_NO_AVX2)
            #define AX_SUPPORT_AVX2
        #endif
    #endif

    /* If at this point we still haven't determined compiler support for the intrinsics just fall back to __has_include. */
    #if !defined(__GNUC__) && !defined(__clang__) && defined(__has_include)
        #if !defined(AX_SUPPORT_SSE)   && !defined(AX_NO_SSE2)   && __has_include(<emmintrin.h>)
            #define AX_SUPPORT_SSE
        #endif
        #if !defined(AX_SUPPORT_AVX2)   && !defined(AX_NO_AVX2)   && __has_include(<immintrin.h>)
            #define AX_SUPPORT_AVX2
        #endif
    #endif

    #if defined(AX_SUPPORT_AVX2) || defined(AX_SUPPORT_AVX)
        #include <immintrin.h>
    #elif defined(AX_SUPPORT_SSE)
        #include <emmintrin.h>
    #endif
#endif

template<typename T>
FINLINE constexpr T PopCount(T x)
{
    // according to intel intrinsic, popcnt instruction is 3 cycle (equal to mulps, addps) 
    // throughput is even double of mulps and addps which is 1.0 (%100)
    // https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html
#ifdef AX_SUPPORT_SSE
	if      constexpr (sizeof(T) == 4) return _mm_popcnt_u32(x);
    else if constexpr (sizeof(T) == 8) return _mm_popcnt_u64(x);
#elif defined(__GNUC__) && !defined(__MINGW32__)
	if      constexpr (sizeof(T) == 4) return __builtin_popcount(x);
	else if constexpr (sizeof(T) == 8) return __builtin_popcountl(x);
#else
	if constexpr (sizeof(T) == 4)
	{
		x =  x - ((x >> 1) & 0x55555555);        // add pairs of bits
		x = (x & 0x33333333) + ((x >> 2) & 0x33333333);  // quads
		x = (x + (x >> 4)) & 0x0F0F0F0F;        // groups of 8
		return (x * 0x01010101) >> 24;          // horizontal sum of bytes	
	}
	else if (sizeof(T) == 8) // standard popcount; from wikipedia
	{
		x -= ((x >> 1) & 0x5555555555555555ull);
		x = (x & 0x3333333333333333ull) + (x >> 2 & 0x3333333333333333ull);
		return ((x + (x >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
	}
#endif
}

template<typename T>
FINLINE constexpr T TrailingZeroCount(T x) 
{
#ifdef _MSC_VER
	if      constexpr (sizeof(T) == 4) return _tzcnt_u32(x);
	else if constexpr (sizeof(T) == 8) return _tzcnt_u64(x);
#elif defined(__GNUC__) && !defined(__MINGW32__)
	if      constexpr (sizeof(T) == 4) return __builtin_ctz(x);
	else if constexpr (sizeof(T) == 8) return __builtin_ctzll(x);
#else
	return PopCount((x & -x) - 1);
#endif
}

template<typename T>
FINLINE constexpr T LeadingZeroCount(T x)
{
#ifdef _MSC_VER
	if      constexpr (sizeof(T) == 4) return _lzcnt_u32(x);
	else if constexpr (sizeof(T) == 8) return _lzcnt_u64(x);
#elif defined(__GNUC__) && !defined(__MINGW32__)
	if      constexpr (sizeof(T) == 4) return __builtin_clz(x);
	else if constexpr (sizeof(T) == 8) return __builtin_clzll(x);
#else
	x |= (x >> 1);
	x |= (x >> 2);
	x |= (x >> 4);
	x |= (x >> 8);
	x |= (x >> 16);
	return 32 - PopCount(x);
#endif
}

template<typename To, typename From>
FINLINE constexpr To BitCast(const From& _Val) {
	return __builtin_bit_cast(To, _Val);
}

template<typename T> FINLINE constexpr T Min(T a, T b) { return a < b ? a : b; }
template<typename T> FINLINE constexpr T Max(T a, T b) { return a > b ? a : b; }
template<typename T> FINLINE constexpr T Clamp(T x, T a, T b) { return Max(a, Min(b, x)); }
template<typename T> FINLINE constexpr T Abs(T x) { return x < T(0) ? -x : x; }

constexpr int NextPowerOf2(int x)
{
    x--;
    x |= x >> 1; x |= x >> 2; x |= x >> 4;
    x |= x >> 8; x |= x >> 16;
    return ++x;
}

constexpr int64_t NextPowerOf2(int64_t x)
{
    x--;
    x |= x >> 1; x |= x >> 2;  x |= x >> 4;
    x |= x >> 8; x |= x >> 16; x |= x >> 32;
    return ++x;
}

template<> FINLINE constexpr float Abs(float x)
{ 
    int ix = BitCast<int>(x) & 0x7FFFFFFF; // every bit except sign mask
    return BitCast<float>(ix);
}

template<> FINLINE constexpr double Abs(double x)
{ 
    uint64 ix = BitCast<uint64>(x) & (~(1ull << 63ull));// every bit except sign mask
    return BitCast<double>(ix);
}

// maybe we should move this to Algorithms.hpp
template<typename T>
inline uint PointerDistance(const T* begin, const T* end)
{
    uint result = 0;
    while (begin++ < end) result++;
	  return result;
}

inline constexpr int CalculateArrayGrowth(int _size)
{
    const int addition = _size / 2;
    if (AX_UNLIKELY(_size > (INT32_MAX - addition))) {
        return INT32_MAX; // growth would overflow
    }
    return _size + addition; // growth is sufficient
}

// enum for skipping initialization
enum EForceInit
{
	ForceInit
};

template<typename A, typename B>
struct Pair
{
    A first; 
    B second;

    Pair(){}
    Pair(A x, B y) : first((A&&)x), second((B&&)y) {}

    bool operator == (const Pair& other) {
        return first == other.first && second == other.second;
    }

    bool operator != (const Pair& other) {
        return first != other.first || second != other.second;
    }
};

template<typename KeyT, typename ValueT>
struct KeyValuePair
{
    KeyT   key{};
    ValueT value{};

    KeyValuePair() {}
    KeyValuePair(KeyT ky, ValueT val) : key((KeyT&&)ky), value((ValueT&&)val) {}
    
    bool operator == (const KeyValuePair& other) {
        return key == other.key && value == other.value;
    }

    bool operator != (const KeyValuePair& other) {
        return key != other.key || value != other.value;
    }
};
