#pragma once

#include "Common.hpp"

AX_NAMESPACE

template<int numBits> struct Bitset
{
public:
	uint64_t bits[MAX(numBits/64, 1)] = { 0 };
	const int size = (numBits / 64) + 1;

	void Set(int idx) {
		int arrIndex = idx / 64;
		int bitIndex = idx & 63;
		bits[arrIndex] |= 1ul << bitIndex;
	}

	bool Get(int idx) {
		int arrIndex = idx / 64;
		int bitIndex = idx & 63;
		return (bits[arrIndex] >> bitIndex) & 1;
	}

	void Reset(int idx) {
		int arrIndex = idx / 64;
		int bitIndex = idx & 63;
		bits[arrIndex] &= ~(1ul << bitIndex);
	}

	void Flip() {
		for (int i = 0; i < size - 1; ++i)
			bits[i] = ~bits[i];
		int lastBits = (1ul << (numBits & 63)) - 1;
		bits[size - 1] ^= lastBits;
	}

	bool All() {
		for (int i = 0; i < size - 1; ++i)
			if (bits[i] != ~0ul) return false;
		int lastBits = (1ul << (numBits & 63)) - 1;
		return (bits[size - 1] & lastBits) == lastBits;
	}

	bool One() {
		for (int i = 0; i < size; ++i)
			if (bits[i] > 0ul) return true;
		return false;
	}

	int Count() {
		int sum = 0;
		for (int i = 0; i < size - 1; ++i)
			sum += __builtin_popcountl(bits[i]);

		int lastBits = numBits & 63;
		for (int i = 0; i < lastBits; ++i)
			sum += ((1ul << i) & bits[size - 1]) > 0;
		return sum;
	}
};

template<typename T>
inline __constexpr void FillN(T* ptr, int len, T val) {
	for (int i = 0; i < len; ++i)
		ptr[i] = val;
}

struct Bitset128
{
	ulong bits[2] = { 0 };

	Bitset128() { bits[0] = 0ul; bits[1] = 0ul; }
	Bitset128(const char* str) {
		int i = 0;
		while (*str) {
			bits[i > 63] |= (1 << (i & 63)) * (*str++ - '0');
		}
	}
	Bitset128(ulong repeat) { bits[0] = repeat; bits[1] = repeat;  }
	Bitset128(ulong a, ulong b) { bits[0] = a; bits[1] = b; }

	bool Get(int idx) const { return (bits[idx > 63] >> idx) & 1; }
	void Set(int idx)       { bits[idx > 63] |= 1ul << (idx & 63); }
	void Reset(int idx)     { bits[idx > 63] &= ~(1ul << (idx & 63)); }

	Bitset128 operator &  (const Bitset128 other) const { return { bits[0] & other.bits[0], bits[1] & other.bits[1] }; }
	Bitset128 operator |  (const Bitset128 other) const { return { bits[0] | other.bits[0], bits[1] | other.bits[1] }; }
	Bitset128 operator ^  (const Bitset128 other) const { return { bits[0] ^ other.bits[0], bits[1] ^ other.bits[1] }; }
	Bitset128 operator &= (const Bitset128 other) { bits[0] &= other.bits[0], bits[1] &= other.bits[1]; return *this; }
	Bitset128 operator |= (const Bitset128 other) { bits[0] |= other.bits[0], bits[1] |= other.bits[1]; return *this; }
	Bitset128 operator ^= (const Bitset128 other) { bits[0] ^= other.bits[0], bits[1] ^= other.bits[1]; return *this; }

	void Clear() { bits[0] = 0ul, bits[1] = 0ul;  }
	void Flip()  { bits[0] = ~bits[0], bits[1] = ~bits[1]; }
	
	bool All()  const { return bits[0] == ~0ul && bits[1] == ~0ul; }
	bool Any()  const { return bits[0] + bits[1] > 0; }
	int Count() const { return PopCount(bits[0]) + PopCount(bits[1]); }
};


AX_END_NAMESPACE