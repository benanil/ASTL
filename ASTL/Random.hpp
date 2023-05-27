// this header file contains random and hash functions
// recommended use PCG for 32 bit random numbers. and Xoroshiro128Plus for 64bit
// for big objects WYHash recommended, 

#pragma once

#include "Math/Math.hpp"

constexpr FINLINE uint FNVHash(const void* ptr, uint64 length) {
	const char* str = (const char*)ptr;
	uint hash = 2166136261u;
	
	for (uint i = 0; i < length; str++, i++) {
		hash *= 16777619u; hash ^= *str;
	}
	return hash;
}

constexpr FINLINE uint64 FNVHash64(const void* ptr, uint64 length) {
	const char* str = (const char*)ptr;
	uint64 hash = 14695981039346656037ull;
	
	for (uint64 i = 0; i < length; str++, i++) {
		hash *= 1099511628211ULL; hash ^= *str;
	}
	return hash;
}

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

// kind of slower than Wang but more accurate
constexpr FINLINE uint32_t TripleHash(uint32_t x) {
	x ^= x >> 17u; x *= 0xed5ad4bbu;
	x ^= x >> 11u; x *= 0xac4c1b51u;
	x ^= x >> 15u; x *= 0x31848babu;
	return x ^ (x >> 14u);
}

constexpr FINLINE uint32_t TripleHashInverse(uint32_t x) {
	x ^= x >> 14u ^ x >> 28u; x *= 0x32b21703u;
	x ^= x >> 15u ^ x >> 30u; x *= 0x469e0db1u;
	x ^= x >> 11u ^ x >> 22u; x *= 0x79a85073u;
	return x ^ (x >> 17u);
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

// #include <immintrin.h>
extern int __cdecl _rdseed32_step(unsigned int *);
#if defined(_M_X64)
extern int __cdecl _rdseed64_step(unsigned __int64 *);
#endif /* defined (_M_X64) */

namespace Random
{
	// these random seeds slower than PCG and MTwister but good choice for random seed
	// also seeds are cryptographic 
	// todo make seeds platform independent
	FINLINE uint Seed32() {
		uint32 result;
		_rdseed32_step(&result);
		return result;
	}

	FINLINE uint64 Seed64() {
		uint64 result;
		_rdseed64_step(&result);
		return result;
	}

	FINLINE float NextFloat01(uint32 next) {
		return float(next >> 8) / 16777216.0f;
	}
	
	FINLINE	float RepatMinMax(uint32 next, float min, float max) {
		return min + (NextFloat01(next) * FAbs(min - max));
	}
	
	FINLINE double NextDouble01(uint64 next) {
		return (next & 0x001FFFFFFFFFFFFF) / 9007199254740992.0;
	}
	
	FINLINE double RepatMinMax(uint64 next, double min, double max) {
		return min + (NextDouble01(next) * FAbs(min - max));
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

	// Copyright (c) 2011, 2013 Mutsuo Saito, Makoto Matsumoto,
	// Hiroshima University and The University of Tokyo. All rights reserved.
	// generated from paper: http://www.math.sci.hiroshima-u.ac.jp/~m-mat/m_MT/ARTICLES/mt.pdf
	// also I don't recommend using more than one instance in a theread
	class MTwister64
	{
		static constexpr int N = 624, M = 367;
		uint64 m_MT[N];
		int m_Index = N + 1;

	public:
		// any non zero integer can be used as a seed
		MTwister64(uint64 seed = 4357ul)
		{
			m_MT[0] = seed & ~0ul;
			for (m_Index = 1; m_Index < N; ++m_Index)
				m_MT[m_Index] = (69069 * m_MT[m_Index - 1]) & ~0ul;
		}

		uint32 Next()
		{
			if (m_Index >= N) GenerateNumbers();
			uint64 x = m_MT[m_Index++];
			x ^= x >> 11;
			x ^= x << 7 & 0x9d2c5680ul;
			x ^= x << 15 & 0xefc60000ul;
			x ^= x >> 18;
			return int(x >> 16);
		}

		uint64 Next64() { return uint32(Next() >> 16); }

	private:

		void GenerateNumbers()
		{
			static const uint64 mag01[2] = { 0x0, 0x9908b0dful };
			int kk = 0; uint64 y;

			while (kk < N - M) // unrolled for chace line optimizations
			{
				y = (m_MT[kk] & 0x80000000ul) | (m_MT[kk + 1] & 0x7ffffffful);
				m_MT[kk] = m_MT[kk + M] ^ (y >> 1) ^ mag01[y & 0x1];
				kk++;
				y = (m_MT[kk] & 0x80000000ul) | (m_MT[kk + 1] & 0x7ffffffful);
				m_MT[kk] = m_MT[kk + M] ^ (y >> 1) ^ mag01[y & 0x1];
				kk++;
				y = (m_MT[kk] & 0x80000000ul) | (m_MT[kk + 1] & 0x7ffffffful);
				m_MT[kk] = m_MT[kk + M] ^ (y >> 1) ^ mag01[y & 0x1];
				kk++;
			}
			kk--;

			while (kk < N - 1)
			{
				y = (m_MT[kk] & 0x80000000ul) | (m_MT[kk + 1] & 0x7ffffffful);
				m_MT[kk] = m_MT[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1];
				++kk;
				y = (m_MT[kk] & 0x80000000ul) | (m_MT[kk + 1] & 0x7ffffffful);
				m_MT[kk] = m_MT[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1];
				++kk;
				y = (m_MT[kk] & 0x80000000ul) | (m_MT[kk + 1] & 0x7ffffffful);
				m_MT[kk] = m_MT[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1];
				++kk;
			}

			y = (m_MT[N - 1] & 0x80000000ul) | (m_MT[0] & 0x7ffffffful);
			m_MT[N - 1] = m_MT[M - 1] ^ (y >> 1) ^ mag01[y & 0x1];

			m_Index = 0;
		}
	};

	class MTwister
	{
		static constexpr int SIZE = 624, PERIOD = 397;
		static constexpr int DIFF = SIZE - PERIOD;
		static constexpr uint32 MAGIC = 0x9908b0df;
		uint32 m_MT[SIZE];
		int m_Index = SIZE;

	public:
		// value doesn't matter
		MTwister(uint32 value = 4586u)
		{
			m_MT[0] = value; m_Index = SIZE;
			for (uint32 i = 1; i < SIZE; ++i)
				m_MT[i] = 0x6c078965 * (m_MT[i - 1] ^ m_MT[i - 1] >> 30) + i;
		}

		uint32 Next()
		{
			if (AX_UNLIKELY(m_Index >= SIZE)) 
				GenerateNumbers();

			uint32 y = m_MT[m_Index++];
			y ^= y >> 11u;
			y ^= y << 7u & 0x9d2c5680u;
			y ^= y << 15u & 0xefc60000u;
			return y ^ (y >> 18u);
		}

		uint64 Next64()
		{
			if (AX_UNLIKELY(m_Index + 1 >= SIZE)) 
				GenerateNumbers();

			uint64 y = m_MT[m_Index++] & (uint64(m_MT[m_Index++]) << 32ul);
			y ^= y >> 11ul;
			y ^= y << 7ul  & (0x9d2c5680ull & (0x9d2c5680ull << 32ull)); 
			y ^= y << 15ul & (0xefc60000ull & (0xefc60000ull << 32ull)); 
			return y ^ (y >> 18ul);
		}

	private:
		void GenerateNumbers()
		{
			size_t i = 0; uint32_t y;

			while (i < DIFF) {
				y = (0x80000000 & m_MT[i]) | (0x7FFFFFFF & m_MT[i + 1]);
				m_MT[i] = m_MT[i + PERIOD] ^ (y >> 1) ^ (((int(y) << 31) >> 31) & MAGIC);
				++i;
				y = (0x80000000 & m_MT[i]) | (0x7FFFFFFF & m_MT[i + 1]);
				m_MT[i] = m_MT[i + PERIOD] ^ (y >> 1) ^ (((int(y) << 31) >> 31) & MAGIC);
				++i;
			}

			while (i < SIZE - 1)
			{
				y = (0x80000000 & m_MT[i]) | (0x7FFFFFFF & m_MT[i + 1]);
				m_MT[i] = m_MT[i - DIFF] ^ (y >> 1) ^ (((int(y) << 31) >> 31) & MAGIC);
				++i;
				y = (0x80000000 & m_MT[i]) | (0x7FFFFFFF & m_MT[i + 1]);
				m_MT[i] = m_MT[i - DIFF] ^ (y >> 1) ^ (((int(y) << 31) >> 31) & MAGIC);
				++i;
				y = (0x80000000 & m_MT[i]) | (0x7FFFFFFF & m_MT[i + 1]);
				m_MT[i] = m_MT[i - DIFF] ^ (y >> 1) ^ (((int(y) << 31) >> 31) & MAGIC);
				++i;
			}

			// i = 623, last step rolls over
			y = (0x80000000 & m_MT[SIZE - 1]) | (0x7FFFFFFF & m_MT[0]);
			m_MT[SIZE - 1] = m_MT[PERIOD - 1] ^ (y >> 1) ^ (((int32_t(y) << 31) >> 31) & MAGIC);

			m_Index = 0;
		}
	};
} // namespace Random end


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
		else if constexpr (sizeof(T) == 8) return MurmurHash(BitCast<uint64>(x));
		else return WYHash::Hash(&x, sizeof(T));
	}
};