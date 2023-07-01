#include "Common.hpp"

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
// usage example, FOR_EACH(fclose, aFile, bFile, cFile);
#define FOR_EACH_(N, what, ...) EXPAND(GLUE(FOR_EACH_, N)(what, __VA_ARGS__))
#define FOR_EACH(what, ...) FOR_EACH_(FOR_EACH_NARG(__VA_ARGS__), what, __VA_ARGS__)
#define DELETE_ALL(...) FOR_EACH(delete, __VA_ARGS__)
#define FREE_ALL(...) FOR_EACH(free, __VA_ARGS__)


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
		y ^= y << 7ul & (0x9d2c5680ull & (0x9d2c5680ull << 32ull));
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