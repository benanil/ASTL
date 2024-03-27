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
		bits[arrIndex] |= 1ull << bitIndex;
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

	// >>
	void ShiftRight(int n)
	{
		// shift that integer and the integers that is left of it
		int numInt = size - 1;
		while (n > 0)
		{
			int numShift = MIN(n, 64);
			uint64_t mask = ~0ull >> (64 - numShift);
			int carry = 0;
			int i = numInt;

			while (i >= 0)
			{
				int tmp = bits[i] & mask;
				bits[i] >>= numShift;
				bits[i] |= carry << (64 - numShift);
				carry = tmp;
				i--;
			}
			n -= 64;
			numInt -= 1;
		}
	}
	
	// <<
	void ShiftLeft(int n)
	{
		int numInt = 0;
		while (n > 0)
		{
			int numShift = MIN(n, 64);
			uint64_t mask = ~0ull << (64 - numShift);
			int carry = 0;
			int i = numInt;

			while (i < size)
			{
				int tmp = bits[i] & mask;
				bits[i] <<= numShift;
				bits[i] |= carry >> (64 - numShift);
				carry = tmp;
				i++;
			}
			n -= 64;
			numInt += 1;
		}
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

	Bitset operator & (const Bitset other) const {
		Bitset res;
		for (int i = 0; i < size; i++)
			res.bits[i] = bits[i] & other.bits[i];
		return res; 
	}

	Bitset operator | (const Bitset other) const {
		Bitset res;
		for (int i = 0; i < size; i++)
			res.bits[i] = bits[i] | other.bits[i];
		return res; 
	}

	Bitset operator ^ (const Bitset other) const {
		Bitset res;
		for (int i = 0; i < size; i++)
			res.bits[i] = bits[i] ^ other.bits[i];
		return res; 
	}

	Bitset operator &= (const Bitset other) {
		for (int i = 0; i < size; i++)
			bits[i] &= other.bits[i];
		return *this;
	}

	Bitset operator |= (const Bitset other) {
		for (int i = 0; i < size; i++)
			bits[i] |= other.bits[i];
		return *this;
	}

	Bitset operator ^= (const Bitset other) {
		for (int i = 0; i < size; i++)
			bits[i] ^= other.bits[i];
		return *this;
	}
};

template<typename T>
inline __constexpr void FillN(T* ptr, int len, T val) {
	for (int i = 0; i < len; ++i)
		ptr[i] = val;
}

struct Bitset128
{
	uint64_t bits[2] = { 0 };

	Bitset128() { bits[0] = 0ul; bits[1] = 0ul; }
	Bitset128(const char* str) {
		int i = 0;
		while (*str) {
			bits[i > 63] |= (1 << (i & 63)) * (*str++ - '0');
		}
	}
	Bitset128(uint64_t repeat) { bits[0] = repeat; bits[1] = repeat;  }
	Bitset128(uint64_t a, uint64_t b) { bits[0] = a; bits[1] = b; }

	bool Get(int idx) const { return (bits[idx > 63] >> idx) & 1; }
	void Set(int idx)       { bits[idx > 63] |= 1ull << (idx & 63); }
	void Reset(int idx)     { bits[idx > 63] &= ~(1ull << (idx & 63)); }

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
	int Count() const { return (int)PopCount64(bits[0]) + (int)PopCount64(bits[1]); }
};


AX_END_NAMESPACE