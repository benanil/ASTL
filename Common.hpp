
// This file is mostly included from other files, it does a lot of things here is what it does
// Compiler Spesiffic features   : VectorCall, force inline, ASSERT, UNREACHABLE,
// Determinate CPU Architecture  : AX_ARM, AX_X86 AX_CPUID, AX_SUPPORT_SSE, AX_SUPPORT_NEON, AVX.. 
// Determinate Operating System  : defines the PLATFORM_XXX macro
// CPP Version Macros            : current cpp version
// Memory Operations             : SmallMemCpy, SmallMemSet, unaligned load
// Bit Operations                : PopCount, ByteSwap, TrailingZeroCount, LeadingZeroCount 
// Basic Math Logical Operations : Min, Max, Clamp, Abs...
// Utilities                     : ArraySize, PointerDistance, Pair, KeyValuePair...
// Almost all macros are here, most of them are agressively inlined for debug and runtime performance

#pragma once

#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <float.h>

typedef unsigned char      uint8;
typedef unsigned short     uint16;
typedef unsigned int       uint32;
typedef unsigned long long uint64;
typedef char      int8;
typedef short     int16;
typedef int       int32;
typedef long long int64;

typedef uint8_t  uchar;
typedef uint16_t ushort;
typedef uint32_t uint;

//------------------------------------------------------------------------
// Compiler Spesiffic features

#ifdef AX_EXPORT
    #define AX_API __declspec(dllexport)
#else
    #define AX_API __declspec(dllimport)
#endif

#ifdef _MSC_VER
    // do nothing   it already has __forceinline
#elif __CLANG__
    #define __forceinline [[clang::always_inline]] 
#elif __GNUC__
    #ifndef __forceinline
        #define __forceinline inline __attribute__((always_inline))
    #endif
#endif

#ifdef _MSC_VER
    #include <intrin.h>
    #define VECTORCALL __vectorcall
#elif __CLANG__
    #define VECTORCALL [[clang::vectorcall]] 
#elif __GNUC__
    #define VECTORCALL  
#endif

#if defined(__cplusplus) &&  __cplusplus >= 201103L
   #define AX_THREAD_LOCAL       thread_local
#elif defined(__GNUC__) && __GNUC__ < 5
   #define AX_THREAD_LOCAL       __thread
#elif defined(_MSC_VER)
   #define AX_THREAD_LOCAL       __declspec(thread)
#elif defined (__STDC_VERSION__) && __STDC_VERSION__ >= 201112L && !defined(__STDC_NO_THREADS__)
   #define AX_THREAD_LOCAL       _Thread_local
#endif

#if defined(__GNUC__)
    #define AX_PACK(decl) decl __attribute__((__packed__))
#elif defined(_MSC_VER)
    #define AX_PACK(decl) __pragma(pack(push, 1)); decl __pragma(pack(pop));
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

#if _MSC_VER
    #define AXGLOBALCONST extern const __declspec(selectany)
#elif defined(__GNUC__) && !defined(__MINGW32__)
    #define AXGLOBALCONST extern const __attribute__((weak))
#else 
    #define AXGLOBALCONST extern const 
#endif

#if defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__clang__)
    #define AX_LIKELY(x) __builtin_expect(x, 1)  
    #define AX_UNLIKELY(x) __builtin_expect(x, 0)
#else
    #define AX_LIKELY(x) (x)   
    #define AX_UNLIKELY(x) (x) 
#endif

// https://nullprogram.com/blog/2022/06/26/
#if defined(_DEBUG) || defined(Debug)
    #if __GNUC__
        #define ASSERT(c) if (!(c)) __builtin_trap()
    #elif _MSC_VER
        #define ASSERT(c) if (!(c)) __debugbreak()
    #else
        #define ASSERT(c) if (!(c)) *(volatile int *)0 = 0
    #endif
#else
    #define ASSERT(c)
#endif

// https://nullprogram.com/blog/2022/06/26/
#if defined(_DEBUG) || defined(Debug)
    #if __GNUC__
        #define ASSERTR(c, r) if (!(c)) { __builtin_trap(); r; }
    #elif _MSC_VER
        #define ASSERTR(c, r) if (!(c)) { __debugbreak(); r; }
    #else
        #define ASSERTR(c, r) if (!(c)) { *(volatile int *)0 = 0; r }
    #endif
#else
    #define ASSERTR(c, r) if (!(c)) { r; }
#endif

#if defined(__has_builtin)
    #define AX_COMPILER_HAS_BUILTIN(x) __has_builtin(x)
#else
    #define AX_COMPILER_HAS_BUILTIN(x) 0
#endif

#if AX_COMPILER_HAS_BUILTIN(__builtin_assume)
    #define AX_ASSUME(x) __builtin_assume(x)
#elif defined(_MSC_VER)
    #define AX_ASSUME(x) __assume(x)
#else
    #define AX_ASSUME(x) (void)(x)
#endif

#if AX_COMPILER_HAS_BUILTIN(__builtin_unreachable)
    #define AX_UNREACHABLE() __builtin_unreachable()
#elif _MSC_VER
    #define AX_UNREACHABLE() __assume(0)
#else
    #define AX_UNREACHABLE() 
#endif

#if AX_COMPILER_HAS_BUILTIN(__builtin_prefetch)
    #define AX_PREFETCH(x) __builtin_prefetch(x)
#elif defined(_MSC_VER)
    #define AX_PREFETCH(x) _mm_prefetch(x, _MM_HINT_T0)
#else
    #define AX_PREFETCH(x)
#endif


//------------------------------------------------------------------------
// Determinate CPU Architecture

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
//  #include <cpuid.h>
//  #define AX_CPUID(num, regs)       __cpuid(num, regs[0], regs[1], regs[2], regs[3])
//  #define AX_CPUID2(num, sub, regs) __cpuid_count(num, sub, regs[0], regs[1], regs[2], regs[3])
    #define AX_CPUID(num, regs)       
#else
    #define AX_CPUID(num, regs) __cpuid(regs, num)
#endif

/* Architecture Detection */
// detection code from mini audio
// you can define AX_NO_SSE2 or AX_NO_AVX2 in order to disable this extensions

#if defined(__x86_64__) || defined(_M_X64)
    #define AX_X64
#elif defined(__i386) || defined(_M_IX86)
    #define AX_X86
#elif defined(_M_ARM) || defined(_M_ARM64) || defined(_M_HYBRID_X86_ARM64) || defined(_M_ARM64EC) || __arm__ || __aarch64__
    #define AX_ARM
#endif

#if defined(AX_ARM)
    #if defined(_MSC_VER) && !defined(__clang__) && (defined(_M_ARM64) || defined(_M_HYBRID_X86_ARM64) || defined(_M_ARM64EC) || defined(__aarch64__))
        #include <arm64_neon.h>
    #else
        #include <arm_neon.h>
    #endif
#endif

// write AX_NO_SSE2 or AX_NO_AVX2 to disable vector instructions

/* Intrinsics Support */
#if (defined(AX_X64) || defined(AX_X86)) && !defined(AX_ARM)
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

//------------------------------------------------------------------------
// Determinate Operating System

// Check windows
#if defined(_WIN32) || defined(_WIN64)
    #define PLATFORM_WINDOWS 1
#endif

// Check unix
#if defined(unix) || defined(__unix__) || defined(__unix) || defined(__APPLE__)
    #define PLATFORM_UNIX 1
#endif

// Check Linux
#if defined(linux) || defined(__linux)
    #define PLATFORM_LINUX 1
#endif

// Check macos
#if defined(__APPLE__)
    #define PLATFORM_UNIX 1
    #define PLATFORM_MACOSX 1
#endif

#if defined(__ANDROID__)
    #define PLATFORM_ANDROID
#endif

//------------------------------------------------------------------------
// CPP version macros

#if defined(__clang__) || defined(__GNUC__)
    #define AX_CPP_VERSION __cplusplus
    #define AX_CPP14 201402L
    #define AX_CPP17 201703L
    #define AX_CPP20 202002L
#elif defined(_MSC_VER)
    #define AX_CPP_VERSION __cplusplus
    #define AX_CPP14 1900
    #define AX_CPP17 1910
    #define AX_CPP20 1920
#endif

#if AX_CPP_VERSION < AX_CPP14
    // below c++ 14 does not support constexpr functions
    #define __constexpr
#else
    #define __constexpr constexpr
#endif

#define inline_constexpr __forceinline __constexpr

#if AX_CPP_VERSION >= AX_CPP17
    #define if_constexpr if constexpr
#else
    #define if_constexpr if
#endif

#if AX_CPP_VERSION < AX_CPP14
    // below c++ 14 does not support constexpr
    #define __const const
#else
    #define __const constexpr
#endif

#ifndef AX_NO_UNROLL
    #if defined(__clang__)
        #define AX_NO_UNROLL _Pragma("clang loop unroll(disable)") _Pragma("clang loop vectorize(disable)")
    #elif defined(__GNUC__) >= 8
        #define AX_NO_UNROLL _Pragma("GCC unroll 0")
    #elif defined(_MSC_VER)
        #define AX_NO_UNROLL __pragma(loop(no_vector))
    #else
        #define AX_NO_UNROLL
    #endif
#endif

//------------------------------------------------------------------------
// Memory Operations, memcpy, memset, unaligned load

#ifdef _MSC_VER
    #define SmallMemCpy(dst, src, size) __movsb((unsigned char*)(dst), (unsigned char*)(src), size);
#else
    #define SmallMemCpy(dst, src, size) __builtin_memcpy(dst, src, size);
#endif

#ifdef _MSC_VER
    #define SmallMemSet(dst, val, size) __stosb((unsigned char*)(dst), val, size);
#else
    #define SmallMemSet(dst, val, size) __builtin_memset(dst, val, size);
#endif

#define MemsetZero(dst, size) SmallMemSet(dst, 0, size)

#if defined(_MSC_VER) && !defined(__clang__)
    #if defined(_M_IX86) //< The __unaligned modifier isn't valid for the x86 platform.
        #define UnalignedLoad64(ptr) *((uint64_t*)(ptr))
    #else
        #define UnalignedLoad64(ptr) *((__unaligned uint64_t*)(ptr))
    #endif
#else
    __forceinline uint64_t UnalignedLoad64(void const* ptr) {
        __attribute__((aligned(1))) uint64_t const *result = (uint64_t const *)ptr;
        return *result;
    }
#endif

#if defined(_MSC_VER) && !defined(__clang__)
    #if defined(_M_IX86) //< The __unaligned modifier isn't valid for the x86 platform.
        #define UnalignedLoad32(ptr) *((uint32_t*)(ptr))
    #else
        #define UnalignedLoad32(ptr) *((__unaligned uint32_t*)(ptr))
    #endif
#else
    __forceinline uint64_t UnalignedLoad32(void const* ptr) {
        __attribute__((aligned(1))) uint32_t const *result = (uint32_t const *)ptr;
        return *result;
    }
#endif

#define UnalignedLoadWord(x) (sizeof(unsigned long long) == 8 ? UnalignedLoad64(x) : UnalignedLoad32(x))

//------------------------------------------------------------------------
// Namespace Begin

// #define AX_USE_NAMESPACE
#ifdef AX_USE_NAMESPACE
    #define AX_NAMESPACE namespace ax {
    #define AX_END_NAMESPACE }
#else
    #define AX_NAMESPACE
    #define AX_END_NAMESPACE
#endif

AX_NAMESPACE

//------------------------------------------------------------------------
// Bit Operations

#if defined(_MSC_VER)     /* Visual Studio */
    #define AX_BSWAP32(x) _byteswap_ulong(x)
#elif (defined (__GNUC__) && (__GNUC__ * 100 + __GNUC_MINOR__ >= 403)) \
|| (defined(__clang__) && __has_builtin(__builtin_bswap32))
    #define AX_BSWAP32(x) __builtin_bswap32(x)
#else
inline uint32_t AX_BSWAP32(uint32_t x) {
    return ((in << 24) & 0xff000000 ) |
           ((in <<  8) & 0x00ff0000 ) |
           ((in >>  8) & 0x0000ff00 ) |
           ((in >> 24) & 0x000000ff );
}
#endif

#if defined(_MSC_VER) 
    #define ByteSwap32(x) _byteswap_uint64(x)
#elif (defined (__GNUC__) && (__GNUC__ * 100 + __GNUC_MINOR__ >= 403)) \
|| (defined(__clang__) && __has_builtin(__builtin_bswap32))
    #define ByteSwap64(x) __builtin_bswap64(x)
#else
inline uint64_t ByteSwap(uint64_t x) {
    return ((x << 56) & 0xff00000000000000ULL) |
           ((x << 40) & 0x00ff000000000000ULL) |
           ((x << 24) & 0x0000ff0000000000ULL) |
           ((x << 8)  & 0x000000ff00000000ULL) |
           ((x >> 8)  & 0x00000000ff000000ULL) |
           ((x >> 24) & 0x0000000000ff0000ULL) |
           ((x >> 40) & 0x000000000000ff00ULL) |
           ((x >> 56) & 0x00000000000000ffULL);
}
#endif

#define ByteSwapWord(x) (sizeof(unsigned long long) == 8 ? ByteSwap64(x) : ByteSwap32(x))

// according to intel intrinsic, popcnt instruction is 3 cycle (equal to mulps, addps) 
// throughput is even double of mulps and addps which is 1.0 (%100)
// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html

#if defined(__ARM_NEON__)
    #define PopCount32(x) vcnt_u8((int8x8_t)x)
#elif defined(AX_SUPPORT_SSE)
    #define PopCount32(x) _mm_popcnt_u32(x)
    #define PopCount64(x) _mm_popcnt_u64(x)
#elif defined(__GNUC__) || !defined(__MINGW32__)
    #define PopCount32(x) __builtin_popcount(x)
    #define PopCount64(x) __builtin_popcountl(x)
#else

inline uint32_t PopCount32(uint32_t x) {
    x =  x - ((x >> 1) & 0x55555555);        // add pairs of bits
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333);  // quads
    x = (x + (x >> 4)) & 0x0F0F0F0F;        // groups of 8
    return (x * 0x01010101) >> 24;          // horizontal sum of bytes	
}

// standard popcount; from wikipedia
inline uint64_t PopCount64(uint64_t x) {
    x -= ((x >> 1) & 0x5555555555555555ull);
    x = (x & 0x3333333333333333ull) + (x >> 2 & 0x3333333333333333ull);
    return ((x + (x >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}
#endif

#ifdef _MSC_VER
    #define TrailingZeroCount32(x) _tzcnt_u32(x)
    #define TrailingZeroCount64(x) _tzcnt_u64(x)
#elif defined(__GNUC__) || !defined(__MINGW32__)
    #define TrailingZeroCount32(x) __builtin_ctz(x)
    #define TrailingZeroCount64(x) __builtin_ctzll(x)
#else
    #define TrailingZeroCount32(x) PopCount64((x & -x) - 1u)
    #define TrailingZeroCount64(x) PopCount64((x & -x) - 1ull)
#endif

#define TrailingZeroCountWord(x) (sizeof(unsigned long long) == 8 ? TrailingZeroCount64(x) : TrailingZeroCount32(x))

#ifdef _MSC_VER
    #define LeadingZeroCount32(x) _lzcnt_u32(x)
    #define LeadingZeroCount64(x) _lzcnt_u64(x)
#elif defined(__GNUC__) || !defined(__MINGW32__)
    #define LeadingZeroCount32(x) __builtin_clz(x)
    #define LeadingZeroCount64(x) __builtin_clzll(x)
#else
template<typename T> inline T LeadingZeroCount64(T x)
{
    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);
    if (sizeof(T) == 8) x |= (x >> 32); 
    return (sizeof(T) * 8) - PopCount(x);
}
#endif


#define LeadingZeroCountWord(x) (sizeof(unsigned long long) == 8 ? LeadingZeroCount64(x) : LeadingZeroCount32(x))

template<typename T>
inline_constexpr T NextSetBit(T* bits)
{
    *bits &= ~T(1);
    T tz = sizeof(T) == 8 ? (T)TrailingZeroCount64((uint64_t)*bits) : (T)TrailingZeroCount32((uint32_t)*bits);
    *bits >>= tz;
    return tz;
}

#define EnumHasBit(_enum, bit) !!(_enum & bit)

template<typename To, typename From>
inline_constexpr To BitCast(const From& _Val) 
{
#if defined(_MSC_VER) && AX_CPP_VERSION < AX_CPP17
  return *reinterpret_cast<const To*>(&_Val);
#else
  return __builtin_bit_cast(To, _Val);
#endif
}

//------------------------------------------------------------------------
// Basic Math Logical Operations, min, max, clamp, abs, 

#ifndef MIN
    #if AX_CPP_VERSION >= AX_CPP17
    template<typename T> inline_constexpr T MIN(T a, T b) { return a < b ? a : b; }
    template<typename T> inline_constexpr T MAX(T a, T b) { return a > b ? a : b; }
    #else
        // using macro if less than 17 because we want this to be constexpr
        #ifndef MIN
            #define MIN(a, b) ((a) < (b) ? (a) : (b))
            #define MAX(a, b) ((a) > (b) ? (a) : (b))
        #endif
    #endif
#endif

template<typename T> 
inline_constexpr T Min3(T a, T b, T c) {
    T res = a < b ? a : b;
    return res < c ? c : res;
}

template<typename T> 
inline_constexpr T Max3(T a, T b, T c) {
    T res = a > b ? a : b;
    return res > c ? res : c;
}

template<typename T> 
inline_constexpr T Clamp(T x, T a, T b) {
    return MAX(a, MIN(b, x));
}

inline_constexpr float Clamp01(float x) {
    return MAX(0.0f, MIN(1.0f, x));
}

inline_constexpr int64_t Abs(int64_t x) 
{
    int64_t temp = x >> 63;
    return (x ^ temp) - temp;
}

inline_constexpr int Abs(int x)
{
    int temp = x >> 31;
    return (x ^ temp) - temp;
}

inline_constexpr float Abs(float x)
{
    int ix = BitCast<int>(x) & 0x7FFFFFFF; // every bit except sign mask
    return BitCast<float>(ix);
}

inline_constexpr double Abs(double x)
{
    uint64_t  ix = BitCast<uint64_t >(x) & (~(1ull << 63ull));// every bit except sign mask
    return BitCast<double>(ix);
}

template<typename T>
inline_constexpr bool IsPowerOfTwo(T x) { 
    return (x != 0) && ((x & (x - 1)) == 0); 
}

inline_constexpr int NextPowerOf2(int x) {
    x--;
    x |= x >> 1; x |= x >> 2; x |= x >> 4;
    x |= x >> 8; x |= x >> 16;
    return ++x;
}

inline_constexpr int64_t NextPowerOf2(int64_t x) {
    x--;
    x |= x >> 1; x |= x >> 2;  x |= x >> 4;
    x |= x >> 8; x |= x >> 16; x |= x >> 32;
    return ++x;
}

//------------------------------------------------------------------------
// Utilities

// change it as is mobile maybe?
inline_constexpr bool IsAndroid()
{
    #ifdef __ANDROID__
    return true;
    #else
    return false;
    #endif
}

template<typename T, int N>
__constexpr int ArraySize(const T (&)[N]) { return N; }

// maybe we should move this to Algorithms.hpp
template<typename T>
inline uint PointerDistance(const T* begin, const T* end)
{
    return uint((char*)end - (char*)begin) / sizeof(T);
}

inline int CalculateArrayGrowth(int _size)
{
    const int addition = _size >> 1;
    if (_size + addition < 0) {
        return INT32_MAX; // growth would overflow
    }
    return _size + addition; // growth is sufficient
}

#if (defined(__GNUC__) || defined(__clang__)) || \
    (defined(_MSC_VER) && (AX_CPP_VERSION >= AX_CPP17))
    #define StringLength(s) (int)__builtin_strlen(s)
#else
typedef unsigned long long unsignedLongLong;
// http://www.lrdev.com/lr/c/strlen.c
inline int StringLength(char const* s)
{
    char const* p = s;
    const unsignedLongLong m = 0x7efefefefefefeffull; 
    const unsignedLongLong n = ~m;

    for (; (unsignedLongLong)p & (sizeof(unsignedLongLong) - 1); p++) 
        if (!*p)
            return (int)(unsignedLongLong)(p - s);

    for (;;) 
    {
        // memory is aligned from now on
        unsignedLongLong i = *(const unsignedLongLong*)p;

        if (!(((i + m) ^ ~i) & n)) {
            p += sizeof(unsignedLongLong);
        }
        else
        {
            for (i = sizeof(unsignedLongLong); i; p++, i--) 
                if (!*p) return (int)(unsignedLongLong)(p - s);
        }
    }
}
#endif

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

    bool operator > (const KeyValuePair& other) {
        return key > other.key;
    }

    bool operator < (const KeyValuePair& other) {
        return key < other.key;
    }

    bool operator == (const KeyValuePair& other) {
        return key == other.key && value == other.value;
    }

    bool operator != (const KeyValuePair& other) {
        return key != other.key || value != other.value;
    }
};

// I hate this but I will use this for selecting allocator in data structures
template <bool B, typename T, typename F>
struct Conditional {
    using Type = T;
};

template <typename T, typename F>
struct Conditional<false, T, F> {
    using Type = F;
};

template <bool B, typename T, typename F>
using ConditionalT = typename Conditional<B, T, F>::Type;

AX_END_NAMESPACE