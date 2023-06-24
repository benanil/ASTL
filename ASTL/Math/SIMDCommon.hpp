#pragma once
#include "Math.hpp"

enum CPUIDBits : int
{
    CPUIDBits_SSE    = (1 << 25),
    CPUIDBits_SSE2   = (1 << 26),
    CPUIDBits_SSE3   = (1 << 9),
    CPUIDBits_SSE4_1 = (1 << 19),
    CPUIDBits_SSE4_2 = (1 << 20),
    CPUIDBits_AVX2   = (1 << 5),
    CPUIDBits_AVX512 = (1 << 16),
};

// for runtime SIMD extension detection.
// sometimes you might need runtime extension detection.
// recommended using global variable like this in a cpp file:
// int g_ax_simd_bits = 0; or inline int g_ax_simd_bits = 0; in here.
// then call AX_InitSIMD_CPUID only once in program lifetime
// then use !!(g_ax_simd_bits & CPUIDBits_SSE2) to get the support value
inline int AX_InitSIMD_CPUID()
{
    int info[4];
    AX_CPUID(1, info);
    int mask = 1;
    mask |= info[3] & (CPUIDBits_SSE | CPUIDBits_SSE2);
    mask |= info[2] & (CPUIDBits_SSE3 | CPUIDBits_SSE4_1 | CPUIDBits_SSE4_2);
    AX_CPUID(7, info);
    mask |= info[1] & (CPUIDBits_AVX2 | CPUIDBits_AVX512);
    return mask;
}

#ifdef AX_SUPPORT_SSE2

AX_ALIGNED(16) struct Vector4UI
{
	union
	{
		struct { uint x, y, z, w; };
		__m128 vec;
	};

	inline operator __m128 () const noexcept { return vec; }
	inline operator __m128i() const noexcept { return _mm_castps_si128(vec); }
	inline operator __m128d() const noexcept { return _mm_castps_pd(vec); }

	Vector4UI() : x(0), y(0), z(0) {}
	Vector4UI(uint _x, uint _y, uint _z, uint _w) : x(_x), y(_y), z(_z), w(_w) {}
};

AX_ALIGNED(16) struct Vector432F
{
	union
	{
		struct { float x, y, z, w; };
		__m128 vec;
	};

	inline operator __m128 () const noexcept { return vec; }
	inline operator __m128i() const noexcept { return _mm_castps_si128(vec); }
	inline operator __m128d() const noexcept { return _mm_castps_pd(vec); }

	Vector432F() : x(0), y(0), z(0), w(0) {}
	constexpr Vector432F(float _x, float _y, float _z, float _w) : x(_x), y(_y), z(_z), w(_w) {}
};

constexpr uint32_t AX_SELECT_0 = 0x00000000;
constexpr uint32_t AX_SELECT_1 = 0xFFFFFFFF;

AXGLOBALCONST Vector4UI g_XMSelect1000 = { AX_SELECT_1, AX_SELECT_0, AX_SELECT_0, AX_SELECT_0 };
AXGLOBALCONST Vector4UI g_XMSelect1100 = { AX_SELECT_1, AX_SELECT_1, AX_SELECT_0, AX_SELECT_0 };
AXGLOBALCONST Vector4UI g_XMSelect1110 = { AX_SELECT_1, AX_SELECT_1, AX_SELECT_1, AX_SELECT_0 };
AXGLOBALCONST Vector4UI g_XMSelect1011 = { AX_SELECT_1, AX_SELECT_0, AX_SELECT_1, AX_SELECT_1 };

AXGLOBALCONST Vector432F g_XMIdentityR0 = { 1.0f, 0.0f, 0.0f, 0.0f };
AXGLOBALCONST Vector432F g_XMIdentityR1 = { 0.0f, 1.0f, 0.0f, 0.0f };
AXGLOBALCONST Vector432F g_XMIdentityR2 = { 0.0f, 0.0f, 1.0f, 0.0f };
AXGLOBALCONST Vector432F g_XMIdentityR3 = { 0.0f, 0.0f, 0.0f, 1.0f };

AXGLOBALCONST Vector4UI g_XMMaskXY = { 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000, 0x00000000 };
AXGLOBALCONST Vector4UI g_XMMask3  = { 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000 };
AXGLOBALCONST Vector4UI g_XMMaskX  = { 0xFFFFFFFF, 0x00000000, 0x00000000, 0x00000000 };
AXGLOBALCONST Vector4UI g_XMMaskY  = { 0x00000000, 0xFFFFFFFF, 0x00000000, 0x00000000 };
AXGLOBALCONST Vector4UI g_XMMaskZ  = { 0x00000000, 0x00000000, 0xFFFFFFFF, 0x00000000 };
AXGLOBALCONST Vector4UI g_XMMaskW  = { 0x00000000, 0x00000000, 0x00000000, 0xFFFFFFFF };

AXGLOBALCONST Vector432F g_XMOne         = { 1.0f, 1.0f, 1.0f, 1.0f };
AXGLOBALCONST Vector432F g_XMOne3        = { 1.0f, 1.0f, 1.0f, 0.0f };
AXGLOBALCONST Vector432F g_XMZero        = { 0.0f, 0.0f, 0.0f, 0.0f };
AXGLOBALCONST Vector432F g_XMTwo         = { 2.0f, 2.0f, 2.0f, 2.0f };
AXGLOBALCONST Vector432F g_XMFour        = { 4.0f, 4.0f, 4.0f, 4.0f };
AXGLOBALCONST Vector432F g_XMSix         = { 6.0f, 6.0f, 6.0f, 6.0f };
AXGLOBALCONST Vector432F g_XMNegativeOne = { -1.0f, -1.0f, -1.0f, -1.0f };
AXGLOBALCONST Vector432F g_XMOneHalf     = { 0.5f, 0.5f, 0.5f, 0.5f };

FINLINE __m128 VECTORCALL SSESelect(const __m128 V1, const __m128 V2, const __m128& Control)
{
	return _mm_or_ps(_mm_andnot_ps(Control, V1), _mm_and_ps(V2, Control));
}

inline constexpr int _mm_shuffle(int fp3, int fp2, int fp1, int fp0)
{
	return (((fp3) << 6) | ((fp2) << 4) | ((fp1) << 2) | ((fp0)));
}

FINLINE __m128 VECTORCALL SSESplatX(const __m128 V1) { return _mm_permute_ps(V1, _mm_shuffle(0, 0, 0, 0)); }
FINLINE __m128 VECTORCALL SSESplatY(const __m128 V1) { return _mm_permute_ps(V1, _mm_shuffle(1, 1, 1, 1)); }
FINLINE __m128 VECTORCALL SSESplatZ(const __m128 V1) { return _mm_permute_ps(V1, _mm_shuffle(2, 2, 2, 2)); }
FINLINE __m128 VECTORCALL SSESplatW(const __m128 V1) { return _mm_permute_ps(V1, _mm_shuffle(3, 3, 3, 3)); }

FINLINE float VECTORCALL Min3(__m128 ab)
{
	__m128 xy = _mm_min_ps(SSESplatX(ab), SSESplatY(ab));
	return _mm_cvtss_f32(_mm_min_ps(xy, SSESplatZ(ab)));
}

FINLINE float VECTORCALL Max3(__m128 ab)
{
	__m128 xy = _mm_max_ps(SSESplatX(ab), SSESplatY(ab));
	return _mm_cvtss_f32(_mm_max_ps(xy, SSESplatZ(ab)));
}

FINLINE void VECTORCALL SSEStoreVector3(float* f, __m128 vec)
{
	_mm_store_ss(f + 0, vec);
	_mm_store_ss(f + 1, _mm_permute_ps(vec, _mm_shuffle(1,1,1,1)));
	_mm_store_ss(f + 2, _mm_permute_ps(vec, _mm_shuffle(2,2,2,2)));
}

FINLINE float VECTORCALL SSEVectorGetX(__m128 V) {
	return _mm_cvtss_f32(V);
}

FINLINE float VECTORCALL SSEVectorGetY(__m128 V) {
	return _mm_cvtss_f32(_mm_shuffle_ps(V, V, _mm_shuffle(1, 1, 1, 1)));
}

FINLINE float VECTORCALL SSEVectorGetZ(__m128 V) {
	return _mm_cvtss_f32(_mm_shuffle_ps(V, V, _mm_shuffle(2, 2, 2, 2)));
}

FINLINE float VECTORCALL SSEVectorGetW(__m128 V) {
	return _mm_cvtss_f32(_mm_shuffle_ps(V, V, _mm_shuffle(3, 3, 3, 3)));
}

FINLINE __m128 VECTORCALL SSEVectorLength(const __m128 V)
{
	return _mm_sqrt_ps(_mm_dp_ps(V, V, 0x7f));
}

FINLINE __m128 VECTORCALL SSEVectorNormalize(const __m128 V)
{
	return _mm_mul_ps(_mm_rsqrt_ps(_mm_dp_ps(V, V, 0x7f)), V);
}

FINLINE float VECTORCALL SSEVectorLengthf(const __m128 v)
{
	return _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(v, v, 0x71)));
}

FINLINE __m128 VECTORCALL SSEVector3Cross(const __m128 vec0, const __m128 vec1)
{
	__m128 tmp0 = _mm_shuffle_ps(vec0, vec0, _MM_SHUFFLE(3,0,2,1));
	__m128 tmp1 = _mm_shuffle_ps(vec1, vec1, _MM_SHUFFLE(3,1,0,2));
	__m128 tmp2 = _mm_mul_ps(tmp0, vec1);
	__m128 tmp3 = _mm_mul_ps(tmp0, tmp1);
	__m128 tmp4 = _mm_shuffle_ps(tmp2, tmp2, _MM_SHUFFLE(3,0,2,1));
	return _mm_sub_ps(tmp3, tmp4);
}

FINLINE __m128 VECTORCALL SSEVector3Dot(const __m128 V1, const __m128 V2)
{
	__m128 vDot = _mm_mul_ps(V1, V2);
	__m128 vTemp = _mm_permute_ps(vDot, _mm_shuffle(2, 1, 2, 1));
	vDot = _mm_add_ss(vDot, vTemp);
	vTemp = _mm_permute_ps(vTemp, _mm_shuffle(1, 1, 1, 1));
	vDot = _mm_add_ss(vDot, vTemp);
	// Splat x
	return _mm_permute_ps(vDot, _mm_shuffle(0, 0, 0, 0));
}

FINLINE __m128i VECTORCALL Multiply64Bit(__m128i ab, __m128i cd)
{
	__m128i b = _mm_srli_epi64(ab, 32);
	__m128i bc = _mm_mul_epu32(b, cd);
	__m128i d = _mm_srli_epi64(cd, 32);
	__m128i ad = _mm_mul_epu32(ab, d);
	__m128i high = _mm_add_epi64(bc, ad);
	high = _mm_slli_epi64(high, 32);
	return _mm_add_epi64(high, _mm_mul_epu32(ab, cd));
}

FINLINE int VECTORCALL hsum_128_epi32avx(__m128i x)
{
	__m128i hi64 = _mm_unpackhi_epi64(x, x); // 3-operand non-destructive AVX lets us save a byte without needing a movdqa
	__m128i sum64 = _mm_add_epi32(hi64, x);
	__m128i hi32 = _mm_shuffle_epi32(sum64, _mm_shuffle(2, 3, 0, 1));    // Swap the low two elements
	__m128i sum32 = _mm_add_epi32(sum64, hi32);
	return _mm_cvtsi128_si32(sum32);       // movd
}

FINLINE double VECTORCALL hsum_128_pdavx(__m128d x)
{
	__m128d hi64 = _mm_unpackhi_pd(x, x); // 3-operand non-destructive AVX lets us save a byte without needing a movdqa
	__m128d sum64 = _mm_add_pd(hi64, x);
	__m128d hi32 = _mm_shuffle_pd(sum64, sum64, _mm_shuffle(2, 3, 0, 1));    // Swap the low two elements
	__m128d sum32 = _mm_add_pd(sum64, hi32);
	return _mm_cvtsd_f64(sum32);       // movd
}

FINLINE float VECTORCALL hsum_ps_sse3(__m128 v) {
	__m128 shuf = _mm_movehdup_ps(v);        // broadcast elements 3,1 to 2,0
	__m128 sums = _mm_add_ps(v, shuf);
	shuf = _mm_movehl_ps(shuf, sums); // high half -> low half
	sums = _mm_add_ss(sums, shuf);
	return _mm_cvtss_f32(sums);
}

#endif // AX_SUPPORT_SSE2

#ifdef AX_SUPPORT_AVX2

FINLINE __m256i VECTORCALL SSESelect(const __m256i V1, const __m256i V2, const __m256i& Control)
{
	__m256i vTemp1 = _mm256_andnot_epi32(Control, V1);
	__m256i vTemp2 = _mm256_and_epi32(V2, Control);
	return _mm256_or_epi32(vTemp1, vTemp2);
}


// https://stackoverflow.com/questions/17863411/sse-multiplication-of-2-64-bit-integers
FINLINE __m256i VECTORCALL Multiply64Bit(__m256i ab, __m256i cd)
{
	__m256i b = _mm256_srli_epi64(ab, 64);
	__m256i bc = _mm256_mul_epu32(b, cd);
	__m256i d = _mm256_srli_epi64(cd, 64);
	__m256i ad = _mm256_mul_epu32(ab, d);
	__m256i high = _mm256_add_epi64(bc, ad);
	high = _mm256_slli_epi64(high, 64);
	return _mm256_add_epi64(high, _mm256_mul_epu32(ab, cd));
}

FINLINE int VECTORCALL hsum_256_epi32(__m256i v)
{
	__m128i sum128 = _mm_add_epi32(_mm256_castsi256_si128(v), _mm256_extracti128_si256(v, 1));
	return hsum_128_epi32avx(sum128);
}

FINLINE int64 VECTORCALL hsum_256_epi64(__m256i v)
{
	return _mm256_cvtsi256_si32(v) + _mm256_extract_epi64(v, 1) + _mm256_extract_epi64(v, 2) + _mm256_extract_epi64(v, 3);
}

FINLINE double VECTORCALL hsum_256_pd(__m256d v)
{
	__m128d sum128 = _mm_add_pd(_mm256_castpd256_pd128(v), _mm256_extractf128_pd(v, 1));
	return hsum_128_pdavx(sum128);
}

// from: Faster Population Counts Using AVX2 Instructions resource paper
FINLINE int64 VECTORCALL popcount256_epi64(__m256i v)
{
	static const __m256i lookup = _mm256_setr_epi8(0, 1, 1, 2, 1, 2, 2, 3, 1, 2,
		2, 3, 2, 3, 3, 4, 0, 1, 1, 2, 1, 2, 2, 3,
		1, 2, 2, 3, 2, 3, 3, 4);
	static const __m256i low_mask = _mm256_set1_epi8(0x0f);
	__m256i lo =  _mm256_and_si256(v, low_mask);
	__m256i hi = _mm256_and_si256(_mm256_srli_epi32(v, 4), low_mask);
	__m256i popcnt1 = _mm256_shuffle_epi8(lookup, lo);
	__m256i popcnt2 = _mm256_shuffle_epi8(lookup, hi);
	__m256i total = _mm256_add_epi8(popcnt1, popcnt2);
	v = _mm256_sad_epu8(total, _mm256_setzero_si256());
	return _mm256_cvtsi256_si32(v) + _mm256_extract_epi64(v, 1) + _mm256_extract_epi64(v, 2) + _mm256_extract_epi64(v, 3);
}

FINLINE __m256i VECTORCALL popcnt256si(__m256i v) // returns 4 64 bit integer that contains pop counts
{
	const __m256i lookup = _mm256_setr_epi8(0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4);
	const __m256i low_mask = _mm256_set1_epi8(0x0f);
	__m256i lo = _mm256_and_si256(v, low_mask);
	__m256i hi = _mm256_and_si256(_mm256_srli_epi32(v, 4), low_mask);
	__m256i popcnt1 = _mm256_shuffle_epi8(lookup, lo);
	__m256i popcnt2 = _mm256_shuffle_epi8(lookup, hi);
	return _mm256_sad_epu8(_mm256_add_epi8(popcnt1, popcnt2), _mm256_setzero_si256());
}

FINLINE float VECTORCALL hsum256_ps_avx(__m256 v) {
	__m128 vlow = _mm256_castps256_ps128(v);
	__m128 vhigh = _mm256_extractf128_ps(v, 1); // high 128
	vlow = _mm_add_ps(vlow, vhigh);     // add the low 128
	return hsum_ps_sse3(vlow);         // and inline the sse3 version, which is optimal for AVX
}
#endif // AX_SUPPORT_AVX2