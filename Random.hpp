
// this header file contains random and hash functions
// generate random numbers: PCG for 32 bit. and Xoroshiro128Plus for 64bit
// Hash Chunk of Data     : WYHash::Hash or MurmurHash64

// uint     WangHash(uint x);
// uint     WangHashInverse(uint x);
// uint64_t MurmurHash(uint64_t x);
// uint64_t MurmurHashInverse(uint64_t x);
// uint     Seed32() 
// uint64_t Seed64()
// float    NextFloat01(uint32 next) 
// float    RepatMinMax(uint32 next, float min, float max) 
// double   NextDouble01(uint64_t next) 
// double   RepatMinMax(uint64_t next, double min, double max)      // usage:
// uint32   RepatMinMax(uint32 next, uint32 min, uint32 max)        // RandomNextFloat01(PCGNext(pcg))
// uint64_t RepatMinMax(uint64_t next, uint64_t min, uint64_t max)  // RepatMINMAX(PCGNext(pcg), 120, 200);
// int      RepatMinMax(int next, int _min, int _max)               // RepatMINMAX(Xoroshiro128Plus(xoro), 120ull, 200ull);
// uint32   PCGNext(PCG& pcg)
// uint     PCG2Next(uint& rng_state)
// void     PCGInitialize(PCG& pcg, uint64_t initstate, uint64_t seed)
// void     PCGInitialize(PCG& pcg, uint64_t seed) 
// void     Xoroshiro128PlusInit(uint64_t  s[2])
// void     Xoroshiro128PlusSeed(uint64_t  s[2], uint64_t  seed)
// uint64_t Xoroshiro128Plus(uint64_t  s[2])
// void     Suffle(T* begin, uint64_t len)
// uint32_t MurmurHash32(const uint8_t* key, size_t len, uint32_t seed)
// uint64_t MurmurHash64(const void * key, int len, uint64_t seed)
// uint64_t WYHash::Hash(void const* key, size_t len)
// uint64_t WYHash::Hash(uint64_t x) 
//#string hash functions
// uint64_t StringToHash64(const char* str, uint64_t len)
// uint     StringToHash(const char* str, uint hash = 0)
// uint     PathToHash(const char* str)

#pragma once

#include "Math/Math.hpp"
#include "Algorithms.hpp"

AX_NAMESPACE 

// Not WangHash actually we can say skeeto hash.
// developed and highly optimized by Chris Wellons
// https://github.com/skeeto/hash-prospector https://nullprogram.com/blog/2018/07/31/
__constexpr __forceinline uint WangHash(uint x) { 
	x ^= x >> 16u; x *= 0x7feb352du;
	x ^= x >> 15u; x *= 0x846ca68bu;
	return x ^ (x >> 16u);
}

// given Wang hash returns input value: 
// WangHash(x) = 234525;
// x = InverseWangHash(234525);
__constexpr __forceinline uint WangHashInverse(uint x)  {
	x ^= x >> 16u; x *= 0x7feb352du;
	x ^= x >> 15u; x *= 0x846ca68bu;
	return x ^ (x >> 16u);
}

__constexpr __forceinline uint64_t MurmurHash(uint64_t x) {
	x ^= x >> 30ULL; x *= 0xbf58476d1ce4e5b9ULL;
	x ^= x >> 27ULL; x *= 0x94d049bb133111ebULL;
	return x ^ (x >> 31ULL);
}

__constexpr __forceinline uint64_t MurmurHashInverse(uint64_t x) {
	x ^= x >> 31ULL ^ x >> 62ULL; x *= 0x319642b2d24d8ec3ULL;
	x ^= x >> 27ULL ^ x >> 54ULL; x *= 0x96de1b173f119089ULL;
	return x ^ (x >> 30ULL ^ x >> 60ULL);
}

// todo find way of random seed generation
#ifndef _MSCVER
#include <immintrin.h> // intrin.h is included defaultly with msvc
#endif

namespace Random
{
	// these random seeds slower than PCG and MTwister but good choice for random seed
	// also seeds are cryptographic 
	__forceinline uint Seed32() {
		uint32 result;
		_rdseed32_step(&result); // or faster __rdtsc
		return result;
	}

	__forceinline uint64_t Seed64() {
		uint64_t result;
		_rdseed64_step(&result);// or faster __rdtsc
		return result;
	}

	__forceinline float NextFloat01(uint32 next) {
		return float(next >> 8) / 16777216.0f;
	}
	
	__forceinline float RepatMinMax(uint32 next, float min, float max) {
		return min + (NextFloat01(next) * Abs(min - max));
	}
	
	__forceinline double NextDouble01(uint64_t next) {
		return (next & 0x001FFFFFFFFFFFFF) / 9007199254740992.0;
	}
	
	__forceinline double RepatMinMax(uint64_t next, double min, double max) {
		return min + (NextDouble01(next) * Abs(min - max));
	}
	
	__forceinline uint32 RepatMinMax(uint32 next, uint32 min, uint32 max) { return min + (next % (max - min)); }
	__forceinline uint64_t RepatMinMax(uint64_t next, uint64_t min, uint64_t max) { return min + (next % (max - min)); }
	__forceinline int    RepatMinMax(int next, int _min, int _max) { return _min + (next % (_max - _min)); }

	// https://www.pcg-random.org/index.html
	// we can also add global state in a cpp file
	// compared to m_MT chace friendly
	struct PCG 
	{
		uint64_t state = 0x853c49e6748fea9bULL;
		uint64_t inc = 0xda3e39cb94b95bdbULL;
	};            
	
	// usage:
	// RandomNextFloat01(PCGNext(pcg))
	// RepatMINMAX(PCGNext(pcg), 120, 200);
	// RepatMINMAX(Xoroshiro128Plus(xoro), 120ull, 200ull);

	__forceinline uint32 PCGNext(PCG& pcg)
	{
		uint64_t oldstate = pcg.state;
		pcg.state = oldstate * 6364136223846793005ULL + (pcg.inc | 1);
		uint64_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
		uint64_t rot = oldstate >> 59u;
#pragma warning(disable : 4146) // unary minus warning fix
		// if you get unary minus error disable sdl checks from msvc settings
		return uint32((xorshifted >> rot) | (xorshifted << ((-rot) & 31)));
	}

	__forceinline uint PCG2Next(uint& rng_state)
	{
		uint state = rng_state;
		rng_state = state * 747796405u + 2891336453u;
		uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
		return (word >> 22u) ^ word;
	}

	__forceinline void PCGInitialize(PCG& pcg, uint64_t initstate, uint64_t seed)
	{
		pcg.inc = (seed << 1u) | 1u;
		pcg.state = initstate;
		PCGNext(pcg);
	}

	__forceinline void PCGInitialize(PCG& pcg, uint64_t seed) 
	{
		pcg.state = 0x853c49e6748fea9bULL;
		pcg.inc = seed << 1 | 1u;
	}

	__forceinline void Xoroshiro128PlusInit(uint64_t  s[2])
	{
		s[0] += Seed64(); s[1] += Seed64();
		s[0] |= 1; // non zero
	}
	
	__forceinline void Xoroshiro128PlusSeed(uint64_t  s[2], uint64_t  seed)
	{
		s[0] |= 1; // non zero
		s[0] = MurmurHash(seed); 
		s[1] = MurmurHash(s[0] ^ (seed * 1099511628211ULL));
	}
	
	// concise hashing function. https://nullprogram.com/blog/2017/09/21/
	__forceinline uint64_t  Xoroshiro128Plus(uint64_t  s[2])
	{
		uint64_t  s0 = s[0];
		uint64_t  s1 = s[1];
		uint64_t  result = s0 + s1;
		s1 ^= s0;
		s[0] = ((s0 << 55) | (s0 >> 9)) ^ s1 ^ (s1 << 14);
		s[1] = (s1 << 36) | (s1 >> 28);
		return result;
	}

	// too see alternative random number generator look at Aditional.hpp for mersene twister pseudo random number generators

	template<typename T>
	inline void Suffle(T* begin, uint64_t len)
	{
		uint64_t  xoro[2];
		Xoroshiro128PlusInit(xoro);
		const uint64_t  halfLen = len / 2;

		// swap %60 of the array
		for (uint64_t i = 0; i < (halfLen + (halfLen / 3)); ++i)
		{
			Swap(begin[Xoroshiro128Plus(xoro) % len],
			     begin[Xoroshiro128Plus(xoro) % len]);
		}
	}
} // namespace Random end

inline uint32_t murmur_32_scramble(uint32_t k) {
	k *= 0xcc9e2d51u;
	k = (k << 15u) | (k >> 17u);
	k *= 0x1b873593u;
	return k;
}

// https://en.wikipedia.org/wiki/MurmurHash
inline uint32_t MurmurHash32(const uint8_t* key, size_t len, uint32_t seed)
{
	uint32_t h = seed;
	uint32_t k;
	/* Read in groups of 4. */
	for (size_t i = len >> 2; i; i--) {
		// Here is a source of differing results across endiannesses.
		// A swap here has no effects on hash properties though.
		SmallMemCpy(&k, key, sizeof(uint32_t));
		key += sizeof(uint32_t);
		h ^= murmur_32_scramble(k);
		h = (h << 13u) | (h >> 19u);
		h = h * 5u + 0xe6546b64u;
	}
	/* Read the rest. */
	k = 0u;
	for (size_t i = len & 3u; i; i--) {
		k <<= 8u;
		k |= key[i - 1u];
	}
	// A swap is *not* necessary here because the preceding loop already
	// places the low bytes in the low places according to whatever endianness
	// we use. Swaps only apply when the memory is copied in a chunk.
	h ^= murmur_32_scramble(k) ^ len; 
	/* Finalize. */
	h ^= h >> 16u;
	h *= 0x7feb352du;
	h ^= h >> 15u;
	h *= 0x846ca68bu;
	h ^= h >> 16u;
	return h;
}

inline uint64_t MurmurHash64(const void * key, int len, uint64_t seed)
{
	const uint64_t m = 0xc6a4a7935bd1e995ULL;
	const int      r = 47;
	uint64_t       h = seed ^ (len * m);

	const uint64_t * data = (const uint64_t *)key;
	const uint64_t * end = data + (len >> 3);

	while (data != end)
	{
		uint64_t k;
		SmallMemCpy(&k, data++, sizeof(uint64));
		k *= m;
		k ^= k >> r;
		k *= m;

		h ^= k;
		h *= m;
	}

	const unsigned char * data2 = (const unsigned char*)data;
	uint64_t  d;
	SmallMemCpy(&d, data, len & 7);
	h ^= d;
	h *= 0xbf58476d1ce4e5b9ULL;
	h ^= h >> 27ULL;
	h *= 0x94d049bb133111ebULL;
	return h ^ (h >> 31ULL);
}

// https://github.com/martinus/unordered_dense/blob/main/include/ankerl/unordered_dense.h
// This is a stripped-down implementation of wyhash: https://github.com/wangyi-fudan/wyhash
// No big-endian support (because different values on different machines don't matter),
// hardcodes seed and the secret, reformattes the code, and clang-tidy fixes.
namespace WYHash
{
	__forceinline void mum(uint64_t * RESTRICT a, uint64_t * RESTRICT b)
	{
#if defined(__SIZEOF_INT128__)
		__uint128_t r = *a;
		r *= *b;
		*a = static_cast<uint64_t >(r);
		*b = static_cast<uint64_t >(r >> 64U);
#elif defined(_MSC_VER) && defined(_M_X64)
		*a = _umul128(*a, *b, b);
#else
		uint64_t  ha = *a >> 32U;
		uint64_t  hb = *b >> 32U;
		uint64_t  la = static_cast<uint32_t>(*a);
		uint64_t  lb = static_cast<uint32_t>(*b);
		uint64_t  hi{}, lo{};
		uint64_t  rh = ha * hb;
		uint64_t  rm0 = ha * lb;
		uint64_t  rm1 = hb * la;
		uint64_t  rl = la * lb;
		uint64_t  t = rl + (rm0 << 32U);
		auto c = static_cast<uint64_t >(t < rl);
		lo = t + (rm1 << 32U);
		c += static_cast<uint64_t >(lo < t);
		hi = rh + (rm0 >> 32U) + (rm1 >> 32U) + c;
		*a = lo;
		*b = hi;
#endif
	}

	__forceinline uint64_t mix(uint64_t a, uint64_t b) { mum(&a, &b); return a ^ b; }
	__forceinline uint64_t r8(const uint8* p) { return *(uint64*)p; }
	__forceinline uint64_t r4(const uint8* p) { return (uint64)*(uint32*)p; }
	// reads 1, 2, or 3 bytes
	__forceinline uint64_t r3(const uint8* p, size_t k) 
	{
		return (uint64_t (p[0]) << 16U) | (uint64_t (p[k >> 1U]) << 8U) | p[k - 1]; 
	}

	// alternative algorithm for this is xxHash which is really efficient
	// https://github.com/Cyan4973/xxHash

	[[nodiscard]] __forceinline uint64_t Hash(void const* key, size_t len)
	{
        const size_t secret[4] = { 0xa0761d6478bd642full,
                                   0xe7037ed1a0b428dbull,
                                   0x8ebc6af09c88c6e3ull,
                                   0x589965cc75374cc3ull };

		uint8 const* p = (uint8 const*)key;
		uint64_t  seed = secret[0];
		uint64_t  a{}, b{};

		if (AX_LIKELY(len <= 16)) {
			if (AX_LIKELY(len >= 4)) {
				a = (r4(p) << 32U) | r4(p + ((len >> 3U) << 2U));
				b = (r4(p + len - 4) << 32U) | r4(p + len - 4 - ((len >> 3U) << 2U));
			} else if (AX_LIKELY(len > 0)) {
				a = r3(p, len);
				b = 0;
			} else {
				a = 0;
				b = 0;
			}
		} else {
			size_t i = len;
			if (AX_UNLIKELY(i > 48)) {
				uint64_t  see1 = seed;
				uint64_t  see2 = seed;
				do {
					seed = mix(r8(p) ^ secret[1], r8(p + 8) ^ seed);
					see1 = mix(r8(p + 16) ^ secret[2], r8(p + 24) ^ see1);
					see2 = mix(r8(p + 32) ^ secret[3], r8(p + 40) ^ see2);
					p += 48;
					i -= 48;
				} while (AX_LIKELY(i > 48));
				seed ^= see1 ^ see2;
			}
			while (AX_UNLIKELY(i > 16)) {
				seed = mix(r8(p) ^ secret[1], r8(p + 8) ^ seed);
				i -= 16;
				p += 16;
			}
			a = r8(p + i - 16);
			b = r8(p + i - 8);
		}

		return mix(secret[1] ^ len, mix(a ^ secret[1], b ^ seed));
	}

	[[nodiscard]] __forceinline uint64_t Hash(uint64_t x) 
	{
		return mix(x, UINT64_C(0x9E3779B97F4A7C15));
	}
}

__constexpr inline uint64_t StringToHash64(const char* str, uint64_t len)
{
	const uint64_t m = 0xc6a4a7935bd1e995ULL;
	const uint64_t r = 47;
	uint64_t h       = 0x9E3779B97F4A7C15ull ^ (len * m);
	uint64_t shift   = 0ull;

	while (len >= 10)
	{
		uint64_t k = 0ull;
    
		while (shift < 60)
		{
			k     |= ulong(ToUpper(*str++) - '0') << shift;
			shift += 6ull;
		}
		// fill missing 4 bits, 10 is random shift to choose for last 4 bit
		k |= (~0xFFFFFFFFFFFFFFF8) & (k << 10ull);

		k *= m;
		k ^= k >> r;
		k *= m;

		h ^= k;
		h *= m;
		shift = 0ull;
		len  -= 10;
	}

	uint64_t d = 0ull;
	while (*str)
	{
		d     |= uint64_t(ToUpper(*str++) - '0') << shift;
		shift += 6ull;
	}

	h ^= d;
	h *= 0xbf58476d1ce4e5b9ULL;
	h ^= h >> 27ULL;
	h *= 0x94d049bb133111ebULL;
	return h ^ (h >> 31ULL);
}

__constexpr inline uint StringToHash(const char* str, uint hash = 0)
{
	while (*str)
		hash = *str++ + (hash << 6u) + (hash << 16u) - hash;
	return hash;
}

__constexpr inline uint PathToHash(const char* str)
{
	uint hash = 0u, idx = 0u, shift = 0u;
	while (str[idx] && idx < 4u)
		hash |= uint(str[idx]) << shift, shift += 8u, idx++;
	return StringToHash(str + idx, WangHash(hash));
}

template<typename T> struct  Hasher
{
	static __forceinline uint64_t Hash(const T& x)
	{
		if_constexpr (sizeof(T) == 4) return uint64(WangHash(BitCast<uint32>(x))) * 0x9ddfea08eb382d69ull;
		else if      (sizeof(T) == 8) return MurmurHash(BitCast<uint64>(x));
		else                          return WYHash::Hash(&x, sizeof(T));
	}
};

AX_END_NAMESPACE 