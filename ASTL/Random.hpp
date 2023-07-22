// this header file contains random and hash functions
// recommended use PCG for 32 bit random numbers. and Xoroshiro128Plus for 64bit
// for big objects WYHash recommended, 

#pragma once

#include "Math/Math.hpp"

// Not WangHash actually we can say skeeto hash.
// developed and highly optimized by Chris Wellons
// https://github.com/skeeto/hash-prospector https://nullprogram.com/blog/2018/07/31/
constexpr FINLINE uint WangHash(uint x) { 
	x ^= x >> 16u; x *= 0x7feb352du;
	x ^= x >> 15u; x *= 0x846ca68bu;
	return x ^ (x >> 16u);
}

// given Wang hash returns input value: 
// WangHash(x) = 234525;
// x = InverseWangHash(234525);
constexpr FINLINE uint WangHashInverse(uint x)  {
	x ^= x >> 16u; x *= 0x7feb352du;
	x ^= x >> 15u; x *= 0x846ca68bu;
	return x ^ (x >> 16u);
}

constexpr FINLINE uint64 MurmurHash(uint64 x) {
	x ^= x >> 30ULL; x *= 0xbf58476d1ce4e5b9ULL;
	x ^= x >> 27ULL; x *= 0x94d049bb133111ebULL;
	return x ^ (x >> 31ULL);
}

constexpr FINLINE uint64 MurmurHashInverse(uint64 x) {
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
	FINLINE uint Seed32() {
		uint32 result;
		_rdseed32_step(&result); // or faster __rdtsc
		return result;
	}

	FINLINE uint64 Seed64() {
		uint64 result;
		_rdseed64_step(&result);// or faster __rdtsc
		return result;
	}

	FINLINE float NextFloat01(uint32 next) {
		return float(next >> 8) / 16777216.0f;
	}
	
	FINLINE	float RepatMinMax(uint32 next, float min, float max) {
		return min + (NextFloat01(next) * Abs(min - max));
	}
	
	FINLINE double NextDouble01(uint64 next) {
		return (next & 0x001FFFFFFFFFFFFF) / 9007199254740992.0;
	}
	
	FINLINE double RepatMinMax(uint64 next, double min, double max) {
		return min + (NextDouble01(next) * Abs(min - max));
	}
	
	FINLINE uint32 RepatMinMax(uint32 next, uint32 min, uint32 max) { return min + (next % (max - min)); }
	FINLINE uint64 RepatMinMax(uint64 next, uint64 min, uint64 max) { return min + (next % (max - min)); }
	FINLINE int    RepatMinMax(int next, int _min, int _max) { return _min + (next % (_max - _min)); }

	// https://www.pcg-random.org/index.html
	// we can also add global state in a cpp file
	// compared to m_MT chace friendly
	struct PCG 
	{
		uint64 state = 0x853c49e6748fea9bULL;
		uint64 inc = 0xda3e39cb94b95bdbULL;
	};            
	
	// usage:
	// RandomNextFloat01(PCGNext(pcg))
	// RepatMinMax(PCGNext(pcg), 120, 200);
	// RepatMinMax(Xoroshiro128Plus(xoro), 120ull, 200ull);

	FINLINE uint32 PCGNext(PCG& pcg)
	{
		uint64 oldstate = pcg.state;
		pcg.state = oldstate * 6364136223846793005ULL + (pcg.inc | 1);
		uint64 xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
		uint64 rot = oldstate >> 59u;
#pragma warning(disable : 4146) // unary minus warning fix
		// if you get unary minus error disable sdl checks from msvc settings
		return uint32((xorshifted >> rot) | (xorshifted << ((-rot) & 31)));
	}

	FINLINE uint PCG2Next(uint& rng_state)
	{
		uint state = rng_state;
		rng_state = state * 747796405u + 2891336453u;
		uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
		return (word >> 22u) ^ word;
	}

	FINLINE void PCGInitialize(PCG& pcg, uint64 initstate, uint64 seed)
	{
		pcg.inc = (seed << 1u) | 1u;
		pcg.state = initstate;
		PCGNext(pcg);
	}

	FINLINE void PCGInitialize(PCG& pcg, uint64 seed) 
	{
		pcg.state = 0x853c49e6748fea9bULL;
		pcg.inc = seed << 1 | 1u;
	}

	FINLINE void Xoroshiro128PlusInit(uint64_t s[2])
	{
		s[0] += Seed64(); s[1] += Seed64();
		s[0] |= 1; // non zero
	}
	
	FINLINE void Xoroshiro128PlusSeed(uint64_t s[2], uint64_t seed)
	{
		s[0] |= 1; // non zero
		s[0] = MurmurHash(seed); 
		s[1] = MurmurHash(s[0] ^ (seed * 1099511628211ULL));
	}
	
	// concise hashing function. https://nullprogram.com/blog/2017/09/21/
	FINLINE uint64_t Xoroshiro128Plus(uint64_t s[2])
	{
		uint64_t s0 = s[0];
		uint64_t s1 = s[1];
		uint64_t result = s0 + s1;
		s1 ^= s0;
		s[0] = ((s0 << 55) | (s0 >> 9)) ^ s1 ^ (s1 << 14);
		s[1] = (s1 << 36) | (s1 >> 28);
		return result;
	}

	// too see alternative random number generator look at Aditional.hpp for mersene twister pseudo random number generators

	template<typename T>
	inline void Suffle(T* begin, uint64 len)
	{
		uint64_t xoro[2];
		Xoroshiro128PlusInit(xoro);
		const uint64_t halfLen = len / 2;

		// swap %60 of the array
		for (uint64 i = 0; i < (halfLen + (halfLen / 3)); ++i)
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

// I recommend to use murmur hash with small strings
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

inline uint64 MurmurHash64(const void * key, int len, uint64 seed)
{
	const uint64 m = 0xc6a4a7935bd1e995ULL;
	const int    r = 47;
	uint64       h = seed ^ (len * m);

	const uint64 * data = (const uint64 *)key;
	const uint64 * end = data + (len >> 3);

	while (data != end)
	{
		uint64 k;
		SmallMemCpy(&k, data++, sizeof(uint64));
		k *= m;
		k ^= k >> r;
		k *= m;

		h ^= k;
		h *= m;
	}

	const unsigned char * data2 = (const unsigned char*)data;
	uint64_t d;
	SmallMemCpy(&d, data, len & 7);
	h ^= d;
	h *= m;
	h ^= h >> r;
	h *= 0x94d049bb133111ebULL;
	h ^= h >> 31ULL;
	return h;
}

// https://github.com/martinus/unordered_dense/blob/main/include/ankerl/unordered_dense.h
// This is a stripped-down implementation of wyhash: https://github.com/wangyi-fudan/wyhash
// No big-endian support (because different values on different machines don't matter),
// hardcodes seed and the secret, reformattes the code, and clang-tidy fixes.
namespace WYHash
{
	FINLINE void mum(uint64_t* RESTRICT a, uint64_t* RESTRICT b)
	{
#if defined(__SIZEOF_INT128__)
		__uint128_t r = *a;
		r *= *b;
		*a = static_cast<uint64_t>(r);
		*b = static_cast<uint64_t>(r >> 64U);
#elif defined(_MSC_VER) && defined(_M_X64)
		*a = _umul128(*a, *b, b);
#else
		uint64_t ha = *a >> 32U;
		uint64_t hb = *b >> 32U;
		uint64_t la = static_cast<uint32_t>(*a);
		uint64_t lb = static_cast<uint32_t>(*b);
		uint64_t hi{}, lo{};
		uint64_t rh = ha * hb;
		uint64_t rm0 = ha * lb;
		uint64_t rm1 = hb * la;
		uint64_t rl = la * lb;
		uint64_t t = rl + (rm0 << 32U);
		auto c = static_cast<uint64_t>(t < rl);
		lo = t + (rm1 << 32U);
		c += static_cast<uint64_t>(lo < t);
		hi = rh + (rm0 >> 32U) + (rm1 >> 32U) + c;
		*a = lo;
		*b = hi;
#endif
	}

	FINLINE uint64 mix(uint64 a, uint64 b) { mum(&a, &b); return a ^ b; }
	FINLINE uint64 r8(const uint8* p) { return *(uint64*)p; }
	FINLINE uint64 r4(const uint8* p) { return (uint64)*(uint32*)p; }
	// reads 1, 2, or 3 bytes
	FINLINE uint64 r3(const uint8* p, size_t k) 
	{
		return (uint64_t(p[0]) << 16U) | (uint64_t(p[k >> 1U]) << 8U) | p[k - 1]; 
	}

	// alternative algorithm for this is xxHash which is really efficient
	// https://github.com/Cyan4973/xxHash

	[[nodiscard]] FINLINE uint64 Hash(void const* key, size_t len)
	{
        const size_t secret[4] = { 0xa0761d6478bd642full,
                                   0xe7037ed1a0b428dbull,
                                   0x8ebc6af09c88c6e3ull,
                                   0x589965cc75374cc3ull };

		uint8 const* p = (uint8 const*)key;
		uint64_t seed = secret[0];
		uint64_t a{}, b{};

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
				uint64_t see1 = seed;
				uint64_t see2 = seed;
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

	[[nodiscard]] FINLINE uint64 Hash(uint64 x) 
	{
		return mix(x, UINT64_C(0x9E3779B97F4A7C15));
	}
}

// use WYHash if you don't need constexpr
constexpr inline ulong StringToHash64(const char* str, ulong hash = 0)
{
	while (*str)
		hash = *str++ + (hash << 6ull) + (hash << 16ull) - hash;
	return hash;
}

constexpr inline ulong PathToHash64(const char* str)
{
	ulong hash = 0u, idx = 0, shift = 0;
	while (str[idx] && idx < 8)
		hash |= ulong(str[idx]) << shift, shift += 8ull, idx++;
	return StringToHash64(str + idx, MurmurHash(hash));
}

constexpr inline uint StringToHash(const char* str, uint hash = 0)
{
	while (*str)
		hash = *str++ + (hash << 6u) + (hash << 16u) - hash;
	return hash;
}

constexpr inline uint PathToHash(const char* str)
{
	uint hash = 0u, idx = 0u, shift = 0u;
	while (str[idx] && idx < 4u)
		hash |= uint(str[idx]) << shift, shift += 8u, idx++;
	return StringToHash(str + idx, WangHash(hash));
}

template<typename T> struct  Hasher
{
	static FINLINE uint64 Hash(const T& x)
	{
		if constexpr (sizeof(T) == 4) return uint64(WangHash(BitCast<uint32>(x))) * 0x9ddfea08eb382d69ull;
		else if      (sizeof(T) == 8) return MurmurHash(BitCast<uint64>(x));
		else                          return WYHash::Hash(&x, sizeof(T));
	}
};