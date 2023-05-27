#pragma once	
#include "Common.hpp"
#include "Math/SIMDCommon.hpp"

template<int numBits> struct Bitset
{
	static constexpr int size = (numBits / 64) + 1;

	ulong bits[size] = { 0 };
	Bitset(const char* str) {
		int currBit = 0, i = 0;
		while (*str) {
			bits[currBit++] |= (1 << (i & 63)) * (*str++ - '0'), currBit = ++i >> 6;
		}
	}

	Bitset() {}

	bool Get(int idx) const { ax_assert(idx < numBits && idx >= 0); return !!(bits[idx >> 6] & (1ul << idx & 63)); }
	void Set(int idx)       { ax_assert(idx < numBits && idx >= 0); bits[idx >> 6] |= 1ul << (idx & 63); }
	void Reset(int idx)     { ax_assert(idx < numBits && idx >= 0); bits[idx >> 6] &= ~(1ul << (idx & 63)); }

	void Clear() { for (int i = 0; i < size; ++i) bits[i] =      0ul; }
	void Flip()  { for (int i = 0; i < size; ++i) bits[i] = ~bits[i]; }

	bool All() const {
		for (int i = 0; i < size - 1; ++i)
			if (bits[i] != ~0ul) return false;
		int lastBits = (1ul << (numBits & 63)) - 1;
		return (bits[size - 1] & lastBits) == lastBits;
	}

	bool Any() const {
		for (int i = 0; i < size-1; ++i)
			if (bits[i] > 0ul) return true;
		int lastBits = (1ul << (numBits & 63)) - 1;
		return (bits[size - 1] & lastBits) > 0ul;
	}

	int Count() const {
		int sum = 0;
		for (int i = 0; i < size-1; ++i)
			sum += AXPopCount64(bits[i]);
        
		int lastBits = numBits & 63;
		for (int i = 0; i < lastBits; ++i)
			sum += ((1ul << i) & bits[size-1]) > 0;
		return sum;
	}
	//  String toString() {
	//  String str; str.Resize(numBits); for (int i = 0, currBit = 0; i < size; ++i) for (int j = 0; j < 64 && currBit < numBits; ++j, ++currBit) str[currBit] = (bits[i] & (1ul << j) ? '1' : '0');
	// 	return str; }
};

template<typename T>
inline constexpr void FillN(T* ptr, int len, T val) {
	for (int i = 0; i < len; ++i)
		ptr[i] = val;
}

struct Bitset128
{
	unsigned long bits[2] = { 0 };

	Bitset128() { bits[0] = 0ul; bits[1] = 0ul; }
	Bitset128(const char* str) {
		int i = 0;
		while (*str) {
			bits[i > 63] |= (1 << (i & 63)) * (*str++ - '0');
		}
	}
	Bitset128(unsigned long repeat) { bits[0] = repeat; bits[1] = repeat;  }
	Bitset128(ulong a, ulong b) { bits[0] = a; bits[1] = b; }
	bool Get(int idx) const { return !!(bits[idx > 63] & (1ul << (idx & 63ul))); }
	void Set(int idx)   { bits[idx > 63] |= 1ul << (idx & 63); }
	void Reset(int idx) { bits[idx > 63] &= ~(1ul << (idx & 63)); }

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

#define _AND(a, b) _mm256_and_si256(a, b)
#define _ROR(a, b) _mm256_or_si256 (a, b)
#define _XOR(a, b) _mm256_xor_si256(a, b)

struct Bitset512
{
	union {
		unsigned long bits[8] = { 0 };
		__m256i sse[2];
	};

	Bitset512() {  }
	Bitset512(const char* str) {
		int currBits = 0, i = 0; Clear();  while (*str) bits[currBits] |= (1 << (i & 63)) * (*str++ - '0'), currBits = ++i >> 6;
	}
	Bitset512(unsigned long repeat) { FillN(bits, 8, repeat); }
	Bitset512(__m256i a, __m256i b) { sse[0] = a; sse[1] = b; }
	Bitset512 operator~ () { return Bitset512(_mm256_xor_si256(sse[0], _mm256_set1_epi32(0xffffffff)), _mm256_xor_si256(sse[1], _mm256_set1_epi32(0xffffffff))); }

	Bitset512 operator &  (const Bitset512 o) const { return { _mm256_and_si256(sse[0], o.sse[0]), _mm256_and_si256(sse[1], o.sse[1]) }; }
	Bitset512 operator |  (const Bitset512 o) const { return { _mm256_or_si256 (sse[0], o.sse[0]), _mm256_or_si256 (sse[1], o.sse[1]) }; }
	Bitset512 operator ^  (const Bitset512 o) const { return { _mm256_xor_si256(sse[0], o.sse[0]), _mm256_xor_si256(sse[1], o.sse[1]) }; }
	Bitset512 operator &= (const Bitset512 o) { sse[0] = _mm256_and_si256(sse[0], o.sse[0]); sse[1] = _mm256_and_si256(sse[1], o.sse[1]); return *this; }
	Bitset512 operator |= (const Bitset512 o) { sse[0] = _mm256_or_si256 (sse[0], o.sse[0]); sse[1] = _mm256_or_si256 (sse[1], o.sse[1]); return *this; }
	Bitset512 operator ^= (const Bitset512 o) { sse[0] = _mm256_xor_si256(sse[0], o.sse[0]); sse[1] = _mm256_xor_si256(sse[1], o.sse[1]); return *this; }

	void And(Bitset512& o) { sse[0] = _mm256_and_si256(sse[0], o.sse[0]); sse[1] = _mm256_and_si256(sse[1], o.sse[1]); }
	void Or (Bitset512& o) { sse[0] = _mm256_or_si256 (sse[0], o.sse[0]); sse[1] = _mm256_or_si256 (sse[1], o.sse[1]); }
	void Xor(Bitset512& o) { sse[0] = _mm256_xor_si256(sse[0], o.sse[0]); sse[1] = _mm256_xor_si256(sse[1], o.sse[1]); }

	bool Get(int idx) const { return !!(bits[idx >> 6] & (1ul << (idx & 63ul))); }
	void Set(int idx) { bits[idx >> 6] |= 1ul << (idx & 63); }
	void Reset(int idx) { bits[idx >> 6] &= ~(1ul << (idx & 63)); }

	void Clear() { sse[0] = sse[1] = _mm256_setzero_si256(); }
	void Flip()  {
		sse[0] = _mm256_xor_si256(sse[0], _mm256_set1_epi32(0xffffffff)); 
		sse[1] = _mm256_xor_si256(sse[1], _mm256_set1_epi32(0xffffffff));
	}

	bool All() const {
		const __m256i full = _mm256_set_epi64x(~0ul, ~0ul, ~0ul, ~0ul);
		return _mm256_movemask_epi8(_mm256_cmpeq_epi64(sse[0], full)) == ~0ul && _mm256_movemask_epi8(_mm256_cmpeq_epi64(sse[1], full)) == ~0ul;
	}
	bool Any() const {
		const __m256i zero = _mm256_setzero_si256();
		return _mm256_movemask_epi8(_mm256_cmpgt_epi64(sse[0], zero)) > 0 || _mm256_movemask_epi8(_mm256_cmpgt_epi64(sse[1], zero)) > 0;
	}

	int Count() const {
		return hsum_256_epi64(_mm256_add_epi64(popcnt256si(sse[0]), popcnt256si(sse[1])));
	}
};

struct Bitset1024
{
	union {
		unsigned long bits[16] = { 0 };
		struct { Bitset512 b1, b2; };
		struct { __m256i v[4]; };
	};
	
	Bitset1024() { Clear(); }
	Bitset1024(const char* str) {
		int currBits = 0, i = 0; Clear();  while (*str) bits[currBits] |= (1 << (i & 63)) * (*str++ - '0'), currBits = ++i >> 6;
	}
	Bitset1024(unsigned long repeat) { FillN(bits, 16, repeat); }
	Bitset1024(Bitset512 a, Bitset512 b) { b1 = (Bitset512&&)a; b2 = (Bitset512&&)b; }

	bool Get(int idx)  const { return !!(bits[idx >> 6] & (1ul << (idx & 63ul))); }
	void Set(int idx)   { bits[idx >> 6] |= 1ul << (idx & 63); }
	void Reset(int idx) { bits[idx >> 6] &= ~(1ul << (idx & 63)); }
	Bitset1024 operator~  () { return {~b1, ~b2}; }

	Bitset1024 operator &  (const Bitset1024& o) const { Bitset1024 r; r.v[0] = _AND(v[0], o.v[0]); r.v[1] = _AND(v[1], o.v[1]); r.v[2] = _AND(v[2], o.v[2]); r.v[3] = _AND(v[3], o.v[3]); return r; }
	Bitset1024 operator |  (const Bitset1024& o) const { Bitset1024 r; r.v[0] = _ROR(v[0], o.v[0]); r.v[1] = _ROR(v[1], o.v[1]); r.v[2] = _ROR(v[2], o.v[2]); r.v[3] = _ROR(v[3], o.v[3]); return r; }
	Bitset1024 operator ^  (const Bitset1024& o) const { Bitset1024 r; r.v[0] = _XOR(v[0], o.v[0]); r.v[1] = _XOR(v[1], o.v[1]); r.v[2] = _XOR(v[2], o.v[2]); r.v[3] = _XOR(v[3], o.v[3]); return r; }

	Bitset1024& operator &= (const Bitset1024& o) { v[0] = _AND(v[0], o.v[0]); v[1] = _AND(v[1], o.v[1]); v[2] = _AND(v[2], o.v[2]); v[3] = _AND(v[3], o.v[3]); return *this; }
	Bitset1024& operator |= (const Bitset1024& o) { v[0] = _ROR(v[0], o.v[0]); v[1] = _ROR(v[1], o.v[1]); v[2] = _ROR(v[2], o.v[2]); v[3] = _ROR(v[3], o.v[3]); return *this; }
	Bitset1024& operator ^= (const Bitset1024& o) { v[0] = _XOR(v[0], o.v[0]); v[1] = _XOR(v[1], o.v[1]); v[2] = _XOR(v[2], o.v[2]); v[3] = _XOR(v[3], o.v[3]); return *this; }

	void And(Bitset1024& o) { v[0] = _AND(v[0], o.v[0]); v[1] = _AND(v[1], o.v[1]); v[2] = _AND(v[2], o.v[2]); v[3] = _AND(v[3], o.v[3]); }
	void Or(Bitset1024&  o) { v[0] = _ROR(v[0], o.v[0]); v[1] = _ROR(v[1], o.v[1]); v[2] = _ROR(v[2], o.v[2]); v[3] = _ROR(v[3], o.v[3]); }
	void Xor(Bitset1024& o) { v[0] = _XOR(v[0], o.v[0]); v[1] = _XOR(v[1], o.v[1]); v[2] = _XOR(v[2], o.v[2]); v[3] = _XOR(v[3], o.v[3]); }

	void Clear() { v[0] = v[1] = v[2] = v[3] = _mm256_setzero_si256(); }
	void Flip()  {
		const __m256i full = _mm256_set1_epi32(0xffffffff);
		v[0] = _mm256_xor_si256(v[0], full); v[1] = _mm256_xor_si256(v[1], full);
		v[2] = _mm256_xor_si256(v[2], full); v[3] = _mm256_xor_si256(v[3], full);
	}
	bool All() const {
		const __m256i full = _mm256_set_epi64x(~0ul, ~0ul, ~0ul, ~0ul);
		return _mm256_movemask_epi8(_mm256_cmpeq_epi64(v[0], full)) == ~0ul && _mm256_movemask_epi8(_mm256_cmpeq_epi64(v[1], full)) == ~0ul && _mm256_movemask_epi8(_mm256_cmpeq_epi64(v[2], full)) == ~0ul && _mm256_movemask_epi8(_mm256_cmpeq_epi64(v[3], full)) == ~0ul;
	}
	bool Any() const {
		const __m256i zero = _mm256_setzero_si256();
		return _mm256_movemask_epi8(_mm256_cmpgt_epi64(v[0], zero)) > 0 || _mm256_movemask_epi8(_mm256_cmpgt_epi64(v[1], zero)) > 0 || _mm256_movemask_epi8(_mm256_cmpgt_epi64(v[2], zero)) > 0 || _mm256_movemask_epi8(_mm256_cmpgt_epi64(v[3], zero)) > 0;
	}
	int64 Count() const
	{
		return hsum_256_epi64( // horizontal_add(popcnt(v0,v1,v2,v3))
			_mm256_add_epi64(_mm256_add_epi64(popcnt256si(v[0]), popcnt256si(v[1])), _mm256_add_epi64(popcnt256si(v[2]), popcnt256si(v[3])))
		);
	}
};

struct Bitset2048
{
	union {
		unsigned long bits[32] = { 0 };
		struct { Bitset1024 b1, b2; };
		struct { __m256i v[8]; };
	};
	
	Bitset2048() { Clear(); }
	Bitset2048(const char* str) {
		int currBits = 0, i = 0; Clear();  while (*str) bits[currBits] |= (1 << (i & 63)) * (*str++ - '0'), currBits = ++i >> 6;
	}
	Bitset2048(unsigned long repeat) { FillN(bits, 32, repeat); }
	Bitset2048(Bitset1024 a, Bitset1024 b) { b1 = a; b2 = b; }

	Bitset2048 operator~  () { return {~b1, ~b2}; }

	Bitset2048 operator &  (const Bitset2048& o) const { Bitset2048 r; r.v[0] = _AND(v[0], o.v[0]); r.v[1] = _AND(v[1], o.v[1]); r.v[2] = _AND(v[2], o.v[2]); r.v[3] = _AND(v[3], o.v[3]); r.v[4] = _AND(v[4], o.v[4]); r.v[5] = _AND(v[5], o.v[5]); r.v[6] = _AND(v[6], o.v[6]); r.v[7] = _AND(v[7], o.v[7]); return r; }
	Bitset2048 operator |  (const Bitset2048& o) const { Bitset2048 r; r.v[0] = _ROR(v[0], o.v[0]); r.v[1] = _ROR(v[1], o.v[1]); r.v[2] = _ROR(v[2], o.v[2]); r.v[3] = _ROR(v[3], o.v[3]); r.v[4] = _ROR(v[4], o.v[4]); r.v[5] = _ROR(v[5], o.v[5]); r.v[6] = _ROR(v[6], o.v[6]); r.v[7] = _ROR(v[7], o.v[7]); return r; }
	Bitset2048 operator ^  (const Bitset2048& o) const { Bitset2048 r; r.v[0] = _XOR(v[0], o.v[0]); r.v[1] = _XOR(v[1], o.v[1]); r.v[2] = _XOR(v[2], o.v[2]); r.v[3] = _XOR(v[3], o.v[3]); r.v[4] = _XOR(v[4], o.v[4]); r.v[5] = _XOR(v[5], o.v[5]); r.v[6] = _XOR(v[6], o.v[6]); r.v[7] = _XOR(v[7], o.v[7]); return r; }

	Bitset2048& operator &=  (const Bitset2048& o) { v[0] = _AND(v[0], o.v[0]); v[1] = _AND(v[1], o.v[1]); v[2] = _AND(v[2], o.v[2]); v[3] = _AND(v[3], o.v[3]); v[4] = _AND(v[4], o.v[4]); v[5] = _AND(v[5], o.v[5]); v[6] = _AND(v[6], o.v[6]); v[7] = _AND(v[7], o.v[7]); return *this; }
	Bitset2048& operator |=  (const Bitset2048& o) { v[0] = _ROR(v[0], o.v[0]); v[1] = _ROR(v[1], o.v[1]); v[2] = _ROR(v[2], o.v[2]); v[3] = _ROR(v[3], o.v[3]); v[4] = _ROR(v[4], o.v[4]); v[5] = _ROR(v[5], o.v[5]); v[6] = _ROR(v[6], o.v[6]); v[7] = _ROR(v[7], o.v[7]); return *this; }
	Bitset2048& operator ^=  (const Bitset2048& o) { v[0] = _XOR(v[0], o.v[0]); v[1] = _XOR(v[1], o.v[1]); v[2] = _XOR(v[2], o.v[2]); v[3] = _XOR(v[3], o.v[3]); v[4] = _XOR(v[4], o.v[4]); v[5] = _XOR(v[5], o.v[5]); v[6] = _XOR(v[6], o.v[6]); v[7] = _XOR(v[7], o.v[7]); return *this; }

	void And(Bitset2048& o) { v[0] = _AND(v[0], o.v[0]); v[1] = _AND(v[1], o.v[1]); v[2] = _AND(v[2], o.v[2]); v[3] = _AND(v[3], o.v[3]); v[4] = _AND(v[4], o.v[4]); v[5] = _AND(v[5], o.v[5]); v[6] = _AND(v[6], o.v[6]); v[7] = _AND(v[7], o.v[7]); }
	void Or (Bitset2048& o) { v[0] = _ROR(v[0], o.v[0]); v[1] = _ROR(v[1], o.v[1]); v[2] = _ROR(v[2], o.v[2]); v[3] = _ROR(v[3], o.v[3]); v[4] = _ROR(v[4], o.v[4]); v[5] = _ROR(v[5], o.v[5]); v[6] = _ROR(v[6], o.v[6]); v[7] = _ROR(v[7], o.v[7]); }
	void Xor(Bitset2048& o) { v[0] = _XOR(v[0], o.v[0]); v[1] = _XOR(v[1], o.v[1]); v[2] = _XOR(v[2], o.v[2]); v[3] = _XOR(v[3], o.v[3]); v[4] = _XOR(v[4], o.v[4]); v[5] = _XOR(v[5], o.v[5]); v[6] = _XOR(v[6], o.v[6]); v[7] = _XOR(v[7], o.v[7]); }

	bool Get(int idx) const { return !!(bits[idx >> 6] & (1ul << (idx & 63ul))); }
	void Set(int idx) { bits[idx >> 6] |= 1ul << (idx & 63); }
	void Reset(int idx) { bits[idx >> 6] &= ~(1ul << (idx & 63)); }

	void Clear() { v[0] = v[1] = v[2] = v[3] = v[4] = v[5] = v[6] = v[7] = _mm256_setzero_si256(); }
	void Flip() {
		const __m256i full = _mm256_set1_epi32(0xffffffff);
		v[0] = _mm256_xor_si256(v[0], full); v[1] = _mm256_xor_si256(v[1], full); v[2] = _mm256_xor_si256(v[2], full); v[3] = _mm256_xor_si256(v[3], full);
		v[4] = _mm256_xor_si256(v[4], full); v[5] = _mm256_xor_si256(v[5], full); v[6] = _mm256_xor_si256(v[6], full); v[7] = _mm256_xor_si256(v[7], full);
	}
	bool All()  const { return b1.All() && b2.All(); }
	bool Any()  const { return b1.Any() && b2.Any(); }
	int Count() const { return b1.Count() + b2.Count(); }
};

struct Bitset4096
{
	union {
		unsigned long bits[64] = { 0 };
		struct { Bitset2048 b1, b2; };
		struct { __m256i v[16]; };
	};
	
	Bitset4096() { Clear(); }
	Bitset4096(const char* str) {
		int currBits = 0, i = 0; Clear();  while (*str) bits[currBits] |= (1 << (i & 63)) * (*str++ - '0'), currBits = ++i >> 6;
	}
	Bitset4096(unsigned long repeat) { FillN(bits, 32, repeat); }
	Bitset4096(Bitset2048 a, Bitset2048 b) { b1 = a; b2 = b; }

	Bitset4096 operator~  () { return {~b1, ~b2}; }

	Bitset4096 operator &  (const Bitset4096& o) const { Bitset4096 r; for (int i = 0; i < 16; ++i) r.v[i] = _mm256_and_si256(v[i], o.v[i]); return r; }
	Bitset4096 operator |  (const Bitset4096& o) const { Bitset4096 r; for (int i = 0; i < 16; ++i) r.v[i] = _mm256_or_si256(v[i], o.v[i]); return r; }
	Bitset4096 operator ^  (const Bitset4096& o) const { Bitset4096 r; for (int i = 0; i < 16; ++i) r.v[i] = _mm256_xor_si256(v[i], o.v[i]); return r; }

	Bitset4096& operator &=  (const Bitset4096& o) { for (int i = 0; i < 16; ++i) v[i] = _mm256_and_si256(v[i], o.v[i]); return *this; }
	Bitset4096& operator |=  (const Bitset4096& o) { for (int i = 0; i < 16; ++i) v[i] = _mm256_or_si256(v[i], o.v[i]); return *this; }
	Bitset4096& operator ^=  (const Bitset4096& o) { for (int i = 0; i < 16; ++i) v[i] = _mm256_xor_si256(v[i], o.v[i]); return *this; }

	void And(Bitset4096& o) { for(int i = 0; i < 16; ++i) v[i] = _mm256_and_si256(v[i], o.v[i]); }
	void Or (Bitset4096& o) { for(int i = 0; i < 16; ++i) v[i] = _mm256_or_si256(v[i], o.v[i]); }
	void Xor(Bitset4096& o) { for(int i = 0; i < 16; ++i) v[i] = _mm256_xor_si256(v[i], o.v[i]); }

	bool Get(int idx) const { return !!(bits[idx >> 6] & (1ul << (idx & 63ul))); }
	void Set(int idx) { bits[idx >> 6] |= 1ul << (idx & 63); }
	void Reset(int idx) { bits[idx >> 6] &= ~(1ul << (idx & 63)); }

	void Clear() { b1.Clear(), b2.Clear(); }
	void Flip()  { b1.Flip(), b2.Flip(); }
	bool All()  const { return b1.All() && b2.All(); }
	bool Any()  const { return b1.Any() && b2.Any(); }
	int Count() const { return b1.Count() + b2.Count(); }
};

#undef _AND
#undef _ROR
#undef _XOR
