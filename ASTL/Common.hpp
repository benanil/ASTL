#pragma once
#define _CRT_NON_CONFORMING_SWPRINTFS 1
#define _CRT_SECURE_NO_WARNINGS 1

#include <stdint.h>

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
#       include <intrin.h>
#   elif __CLANG__
#       define FINLINE [[clang::always_inline]] 
#	elif __GNUC__
#       define FINLINE  __attribute__((always_inline))
#   endif
#endif

#ifndef VECTORCALL
#   ifdef _MSC_VER
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
#    define ax_assert(c) if (!(c)) __builtin_trap()
#  elif _MSC_VER
#    define ax_assert(c) if (!(c)) __debugbreak()
#  else
#    define ax_assert(c) if (!(c)) *(volatile int *)0 = 0
#  endif
#else
#  define ax_assert(c)
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

#define EXPAND(x) x
#define FOR_EACH_1(what, x, ...) what(x)
#define FOR_EACH_2(what, x, ...) what(x); EXPAND(FOR_EACH_1(what, __VA_ARGS__))
#define FOR_EACH_3(what, x, ...) what(x); EXPAND(FOR_EACH_2(what, __VA_ARGS__))
#define FOR_EACH_4(what, x, ...) what(x); EXPAND(FOR_EACH_3(what, __VA_ARGS__))
#define FOR_EACH_5(what, x, ...) what(x); EXPAND(FOR_EACH_4(what, __VA_ARGS__))
#define FOR_EACH_6(what, x, ...) what(x); EXPAND(FOR_EACH_5(what, __VA_ARGS__))

#define FOR_EACH_NARG(...) FOR_EACH_NARG_(__VA_ARGS__, FOR_EACH_RSEQ_N())
#define FOR_EACH_NARG_(...) EXPAND(FOR_EACH_ARG_N(__VA_ARGS__))
#define FOR_EACH_ARG_N(_1, _2, _3, _4, _5, _6, N, ...) N
#define FOR_EACH_RSEQ_N() 6, 5, 4, 3, 2, 1, 0

#define GLUE(x, y) x##y
#define FOR_EACH_(N, what, ...) EXPAND(GLUE(FOR_EACH_, N)(what, __VA_ARGS__))
#define FOR_EACH(what, ...) FOR_EACH_(FOR_EACH_NARG(__VA_ARGS__), what, __VA_ARGS__)
#define DELETE_ALL(...) FOR_EACH(delete, __VA_ARGS__)
#define FREE_ALL(...) FOR_EACH(free, __VA_ARGS__)

#define ENUM_FLAGS(ENUMNAME, ENUMTYPE) \
inline ENUMNAME& operator |= (ENUMNAME& a, ENUMNAME b) noexcept { return (ENUMNAME&)(((ENUMTYPE&)a) |= ((ENUMTYPE)b)); } \
inline ENUMNAME& operator &= (ENUMNAME& a, ENUMNAME b) noexcept { return (ENUMNAME&)(((ENUMTYPE&)a) &= ((ENUMTYPE)b)); } \
inline ENUMNAME& operator ^= (ENUMNAME& a, ENUMNAME b) noexcept { return (ENUMNAME&)(((ENUMTYPE&)a) ^= ((ENUMTYPE)b)); } \
inline constexpr ENUMNAME operator | (ENUMNAME a, ENUMNAME b) noexcept { return ENUMNAME(((ENUMTYPE)a) | ((ENUMTYPE)b));		} \
inline constexpr ENUMNAME operator & (ENUMNAME a, ENUMNAME b) noexcept { return ENUMNAME(((ENUMTYPE)a) & ((ENUMTYPE)b));		} \
inline constexpr ENUMNAME operator ~ (ENUMNAME a)			  noexcept { return ENUMNAME(~((ENUMTYPE)a));						} \
inline constexpr ENUMNAME operator ^ (ENUMNAME a, ENUMNAME b) noexcept { return ENUMNAME(((ENUMTYPE)a) ^ (ENUMTYPE)b);		} 

template<typename T>
_NODISCARD FINLINE constexpr T PopCount(T x) noexcept
{
#ifdef _MSC_VER
	if      constexpr (sizeof(T) == 4) return __popcnt(x);
	else if constexpr (sizeof(T) == 8) return __popcnt64(x);
#elif defined(__GNUC__) && !defined(__MINGW32__)
	if      constexpr (sizeof(T) == 4) return __builtin_popcount(x);
	else if constexpr (sizeof(T) == 8) return __builtin_popcountl(x);
#else
	if constexpr (sizeof(T) == 4)
	{
		i = i - ((i >> 1) & 0x55555555);        // add pairs of bits
		i = (i & 0x33333333) + ((i >> 2) & 0x33333333);  // quads
		i = (i + (i >> 4)) & 0x0F0F0F0F;        // groups of 8
		return (i * 0x01010101) >> 24;          // horizontal sum of bytes	
	}
	else if (sizeof(T) == 8) // standard popcount; from wikipedia
	{
		i -= ((i >> 1) & 0x5555555555555555ull);
		i = (i & 0x3333333333333333ull) + (i >> 2 & 0x3333333333333333ull);
		return ((i + (i >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
	}
#endif
}

template<typename T>
_NODISCARD FINLINE constexpr T TrailingZeroCount(T x) noexcept
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
_NODISCARD FINLINE constexpr T LeadingZeroCount(T x) noexcept
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
_NODISCARD FINLINE constexpr To BitCast(const From& _Val) noexcept {
	return __builtin_bit_cast(To, _Val);
}

typedef int32_t int32;
typedef int64_t int64;

typedef uint16_t ushort;
typedef uint32_t uint  ;
typedef uint64_t ulong ;

typedef uint8_t uint8;
typedef uint16_t uint16;
typedef uint32_t uint32;
typedef uint64_t uint64;

template<typename T> FINLINE constexpr T Max(T a, T b) { return a > b ? a : b; }
template<typename T> FINLINE constexpr T Min(T a, T b) { return a < b ? a : b; }
template<typename T> FINLINE constexpr T Clamp(T x, T a, T b) { return Max(a, Min(b, x)); }
FINLINE constexpr _NODISCARD int Abs(int x) { return x < 0 ? -x : x; }

// maybe we should move this to Algorithms.hpp
template<typename T>
inline uint PointerDistance(const T* begin, const T* end)
{
    uint result = 0;
    while (begin++ < end) result++;
	  return result;
}

/*for qsort*/ template<typename T>
inline int QLess(const void* a, const void* b) { return *(T*)a < *(T*)b; }
/*for qsort*/ template<typename T>
inline int QGreater(const void* a, const void* b) { return *(T*)a > *(T*)b; }

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

    bool operator == (const Pair& other) {
        return first == other.first && second == other.second;
    }

    bool operator != (const Pair& other) {
        return first != other.first || second != other.second;
    }

    Pair(){}
    Pair(A x, B y) : first((A&&)x), second((B&&)y) {}
};

template<typename KeyT, typename ValueT>
struct KeyValuePair
{
    KeyT key{};
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