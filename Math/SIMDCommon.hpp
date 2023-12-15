#pragma once
#include "Math.hpp"

#if defined(__GNUC__) || defined(__clang__)
    
#endif

AX_NAMESPACE 

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

#ifdef AX_SUPPORT_SSE

#define MakeShuffleMask(x,y,z,w)           (x | (y<<2) | (z<<4) | (w<<6))
// vec(0, 1, 2, 3) -> (vec[x], vec[y], vec[z], vec[w])
#define VecSwizzleMask(vec, mask)        _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(vec), mask))
#define VecSwizzle(vec, x, y, z, w)      VecSwizzleMask(vec, MakeShuffleMask(x,y,z,w))
#define VecSwizzle1(vec, x)              VecSwizzleMask(vec, MakeShuffleMask(x,x,x,x))
// special swizzle                   
#define VecSwizzle_0022(vec)             _mm_moveldup_ps(vec)
#define VecSwizzle_1133(vec)             _mm_movehdup_ps(vec)
	
// return (vec1[x], vec1[y], vec2[z], vec2[w])
#define VecShuffle(vec1, vec2, x,y,z,w)  _mm_shuffle_ps(vec1, vec2, MakeShuffleMask(x,y,z,w))
// special shuffle                  
#define VecShuffle_0101(vec1, vec2)      _mm_movelh_ps(vec1, vec2)
#define VecShuffle_2323(vec1, vec2)      _mm_movehl_ps(vec2, vec1)

constexpr uint32_t AX_SELECT_0 = 0x00000000;
constexpr uint32_t AX_SELECT_1 = 0xFFFFFFFF;

#define g_XSelect1000 _mm_castsi128_ps(_mm_setr_epi32(AX_SELECT_1, AX_SELECT_0, AX_SELECT_0, AX_SELECT_0))
#define g_XSelect1100 _mm_castsi128_ps(_mm_setr_epi32(AX_SELECT_1, AX_SELECT_1, AX_SELECT_0, AX_SELECT_0))
#define g_XSelect1110 _mm_castsi128_ps(_mm_setr_epi32(AX_SELECT_1, AX_SELECT_1, AX_SELECT_1, AX_SELECT_0))
#define g_XSelect1011 _mm_castsi128_ps(_mm_setr_epi32(AX_SELECT_1, AX_SELECT_0, AX_SELECT_1, AX_SELECT_1))

#define g_XIdentityR0 _mm_setr_ps(1.0f, 0.0f, 0.0f, 0.0f)
#define g_XIdentityR1 _mm_setr_ps(0.0f, 1.0f, 0.0f, 0.0f)
#define g_XIdentityR2 _mm_setr_ps(0.0f, 0.0f, 1.0f, 0.0f)
#define g_XIdentityR3 _mm_setr_ps(0.0f, 0.0f, 0.0f, 1.0f)

#define g_XMaskXY _mm_castsi128_ps(_mm_setr_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0x00000000, 0x00000000))
#define g_XMask3  _mm_castsi128_ps(_mm_setr_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000))
#define g_XMaskX  _mm_castsi128_ps(_mm_setr_epi32(0xFFFFFFFF, 0x00000000, 0x00000000, 0x00000000))
#define g_XMaskY  _mm_castsi128_ps(_mm_setr_epi32(0x00000000, 0xFFFFFFFF, 0x00000000, 0x00000000))
#define g_XMaskZ  _mm_castsi128_ps(_mm_setr_epi32(0x00000000, 0x00000000, 0xFFFFFFFF, 0x00000000))
#define g_XMaskW  _mm_castsi128_ps(_mm_setr_epi32(0x00000000, 0x00000000, 0x00000000, 0xFFFFFFFF))

#define g_XOne         _mm_setr_ps( 1.0f, 1.0f, 1.0f, 1.0f)
#define g_XNegativeOne _mm_setr_ps( -1.0f, -1.0f, -1.0f, -1.0f)

__forceinline __m128 VECTORCALL SSESelect(const __m128 V1, const __m128 V2, const __m128& Control)
{
    return _mm_blendv_ps(V1, V2, Control);
    // return _mm_or_ps(_mm_andnot_ps(Control, V1), _mm_and_ps(V2, Control));
}

inline __constexpr int _mm_shuffle(int fp3, int fp2, int fp1, int fp0)
{
	return (((fp3) << 6) | ((fp2) << 4) | ((fp1) << 2) | ((fp0)));
}

__forceinline __m128 VECTORCALL SSESplatX(const __m128 V1) { return _mm_permute_ps(V1, _mm_shuffle(0, 0, 0, 0)); }
__forceinline __m128 VECTORCALL SSESplatY(const __m128 V1) { return _mm_permute_ps(V1, _mm_shuffle(1, 1, 1, 1)); }
__forceinline __m128 VECTORCALL SSESplatZ(const __m128 V1) { return _mm_permute_ps(V1, _mm_shuffle(2, 2, 2, 2)); }
__forceinline __m128 VECTORCALL SSESplatW(const __m128 V1) { return _mm_permute_ps(V1, _mm_shuffle(3, 3, 3, 3)); }

__forceinline float VECTORCALL Min3(__m128 ab)
{
	__m128 xy = _mm_min_ps(SSESplatX(ab), SSESplatY(ab));
	return _mm_cvtss_f32(_mm_min_ps(xy, SSESplatZ(ab)));
}

__forceinline float VECTORCALL Max3(__m128 ab)
{
	__m128 xy = _mm_max_ps(SSESplatX(ab), SSESplatY(ab));
	return _mm_cvtss_f32(_mm_max_ps(xy, SSESplatZ(ab)));
}

__forceinline void VECTORCALL SSEStoreVector3(float* f, __m128 vec)
{
	_mm_store_ss(f + 0, vec);
	_mm_store_ss(f + 1, _mm_permute_ps(vec, _mm_shuffle(1,1,1,1)));
	_mm_store_ss(f + 2, _mm_permute_ps(vec, _mm_shuffle(2,2,2,2)));
}

__forceinline float VECTORCALL SSEVectorGetX(__m128 V) {
	return _mm_cvtss_f32(V);
}

__forceinline float VECTORCALL SSEVectorGetY(__m128 V) {
	return _mm_cvtss_f32(_mm_shuffle_ps(V, V, _mm_shuffle(1, 1, 1, 1)));
}

__forceinline float VECTORCALL SSEVectorGetZ(__m128 V) {
	return _mm_cvtss_f32(_mm_shuffle_ps(V, V, _mm_shuffle(2, 2, 2, 2)));
}

__forceinline float VECTORCALL SSEVectorGetW(__m128 V) {
	return _mm_cvtss_f32(_mm_shuffle_ps(V, V, _mm_shuffle(3, 3, 3, 3)));
}

__forceinline __m128 VECTORCALL SSEVectorLength(const __m128 V)
{
	return _mm_sqrt_ps(_mm_dp_ps(V, V, 0x7f));
}

__forceinline __m128 VECTORCALL SSEVectorDistance(const __m128 A, const __m128 B)
{
  __m128 diff = _mm_sub_ps(A, B);
  return _mm_sqrt_ps(_mm_dp_ps(diff, diff, 0x7f));
}

__forceinline __m128 VECTORCALL SSEVectorNormalize(const __m128 V)
{
	return _mm_mul_ps(_mm_rsqrt_ps(_mm_dp_ps(V, V, 0x7f)), V);
}

__forceinline float VECTORCALL SSEVectorLengthf(const __m128 v)
{
	return _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(v, v, 0x71)));
}

__forceinline __m128 VECTORCALL SSEVector3Cross(const __m128 vec0, const __m128 vec1)
{
	__m128 tmp0 = _mm_shuffle_ps(vec0, vec0, _MM_SHUFFLE(3,0,2,1));
	__m128 tmp1 = _mm_shuffle_ps(vec1, vec1, _MM_SHUFFLE(3,1,0,2));
	__m128 tmp2 = _mm_mul_ps(tmp0, vec1);
	__m128 tmp3 = _mm_mul_ps(tmp0, tmp1);
	__m128 tmp4 = _mm_shuffle_ps(tmp2, tmp2, _MM_SHUFFLE(3,0,2,1));
	return _mm_sub_ps(tmp3, tmp4);
}

__forceinline __m128 VECTORCALL SSEVector3Dot(const __m128 V1, const __m128 V2)
{
	__m128 vDot = _mm_mul_ps(V1, V2);
	__m128 vTemp = _mm_permute_ps(vDot, _mm_shuffle(2, 1, 2, 1));
	vDot = _mm_add_ss(vDot, vTemp);
	vTemp = _mm_permute_ps(vTemp, _mm_shuffle(1, 1, 1, 1));
	vDot = _mm_add_ss(vDot, vTemp);
	// Splat x
	return _mm_permute_ps(vDot, _mm_shuffle(0, 0, 0, 0));
}

__forceinline __m128i VECTORCALL Multiply64Bit(__m128i ab, __m128i cd)
{
	__m128i b = _mm_srli_epi64(ab, 32);
	__m128i bc = _mm_mul_epu32(b, cd);
	__m128i d = _mm_srli_epi64(cd, 32);
	__m128i ad = _mm_mul_epu32(ab, d);
	__m128i high = _mm_add_epi64(bc, ad);
	high = _mm_slli_epi64(high, 32);
	return _mm_add_epi64(high, _mm_mul_epu32(ab, cd));
}

__forceinline int VECTORCALL hsum_128_epi32avx(__m128i x)
{
	__m128i hi64 = _mm_unpackhi_epi64(x, x); // 3-operand non-destructive AVX lets us save a byte without needing a movdqa
	__m128i sum64 = _mm_add_epi32(hi64, x);
	__m128i hi32 = _mm_shuffle_epi32(sum64, _mm_shuffle(2, 3, 0, 1));    // Swap the low two elements
	__m128i sum32 = _mm_add_epi32(sum64, hi32);
	return _mm_cvtsi128_si32(sum32);       // movd
}

__forceinline double VECTORCALL hsum_128_pdavx(__m128d x)
{
	__m128d hi64 = _mm_unpackhi_pd(x, x); // 3-operand non-destructive AVX lets us save a byte without needing a movdqa
	__m128d sum64 = _mm_add_pd(hi64, x);
	// todo: fix
	//__m128d hi32 = _mm_shuffle_pd(sum64, sum64, _mm_shuffle_pd(2, 3, 0, 1));    // Swap the low two elements
	//__m128d sum32 = _mm_add_pd(sum64, hi32);
	//return _mm_cvtsd_f64(sum32);       // movd
	return 0.0;
}

__forceinline float VECTORCALL hsum_ps_sse3(__m128 v) {
	__m128 shuf = _mm_movehdup_ps(v);        // broadcast elements 3,1 to 2,0
	__m128 sums = _mm_add_ps(v, shuf);
	shuf = _mm_movehl_ps(shuf, sums); // high half -> low half
	sums = _mm_add_ss(sums, shuf);
	return _mm_cvtss_f32(sums);
}

__forceinline __m128 VECTORCALL _mm_fabs_ps(__m128 x)
{
    __m128 y = _mm_cmple_ps(x, _mm_setzero_ps());
    return SSESelect(x, _mm_sub_ps(_mm_setzero_ps(), x), y);
}

__forceinline __m128 VECTORCALL _mm_copysign_ps(__m128 x, __m128 y)
{
    return _mm_or_ps(_mm_and_ps(x, _mm_set1_ps(0x7fffffff)),
                     _mm_and_ps(y, _mm_set1_ps(0x80000000)));
}

#if defined(__GNUC__) || defined(__clang__)

inline __m128 VECTORCALL _mm_sin_ps(__m128 x)
{ 
  __m128 xx = _mm_mul_ps(x, _mm_mul_ps(x, x));
  __m128 t  = _mm_sub_ps(x, _mm_mul_ps(xx, _mm_set1_ps(0.16666666666f))); 
  
  xx = _mm_mul_ps(x, _mm_mul_ps(x, x));
  t = _mm_add_ps(t, _mm_mul_ps(xx, _mm_set1_ps(0.00833333333f)));
  
  xx = _mm_mul_ps(x, _mm_mul_ps(x, x));
  t  = _mm_sub_ps(t, _mm_mul_ps(xx, _mm_set1_ps(0.00019841269f)));
  
  xx = _mm_mul_ps(x, _mm_mul_ps(x, x));
  t = _mm_add_ps(t, _mm_mul_ps(xx, _mm_set1_ps(2.75573e-06f)));
  return t;
}

inline __m128 VECTORCALL _mm_cos_ps(__m128 x)
{
  __m128 xx = _mm_mul_ps(x, x);
  __m128 t  = _mm_sub_ps(_mm_set1_ps(1.0f), _mm_mul_ps(xx, _mm_set1_ps(0.5f))); 
  
  xx = _mm_mul_ps(x, _mm_mul_ps(x, x));
  t  = _mm_add_ps(t, _mm_mul_ps(xx, _mm_set1_ps(0.04166666666f)));
  
  xx = _mm_mul_ps(x, _mm_mul_ps(x, x));
  t  = _mm_sub_ps(t, _mm_mul_ps(xx, _mm_set1_ps(0.00138888888f)));
  
  xx = _mm_mul_ps(x, _mm_mul_ps(x, x));
  t  = _mm_add_ps(t, _mm_mul_ps(xx, _mm_set1_ps(2.48016e-05f)));
  return t;
}

__forceinline __m128 VECTORCALL _mm_atan_ps(__m128 x)
{
  static const float sa1 =  0.99997726f, sa3 = -0.33262347f, sa5  = 0.19354346f,
                     sa7 = -0.11643287f, sa9 =  0.05265332f, sa11 = -0.01172120f;
    
  const __m128 xx = _mm_mul_ps(x, x);
  // (a9 + x_sq * a11
  __m128 y = _mm_fmadd_ps(xx, _mm_set_ps1(sa11), _mm_set_ps1(sa9));
  y = _mm_fmadd_ps(xx, y, _mm_set_ps1(sa7));
  y = _mm_fmadd_ps(xx, y, _mm_set_ps1(sa5));
  y = _mm_fmadd_ps(xx, y, _mm_set_ps1(sa3));
  y = _mm_fmadd_ps(xx, y, _mm_set_ps1(sa1));
  return _mm_mul_ps(x, y);
}

inline __m128 VECTORCALL _mm_atan2_ps(__m128 y, __m128 x)
{
  __m128 ay = _mm_fabs_ps(y), ax = _mm_fabs_ps(x);
  __m128 invert = _mm_cmpgt_ps(ay, ax);
  __m128 z = SSESelect(_mm_div_ps(ax, ay), _mm_div_ps(ay, ax), invert);
  __m128 th = _mm_atan_ps(z);
  th = SSESelect(th, _mm_sub_ps(_mm_set1_ps(PIDiv2), th), invert);
  th = SSESelect(th, _mm_sub_ps(_mm_set1_ps(PI), th), _mm_cmplt_ps(x, _mm_setzero_ps()));
  return _mm_copysign_ps(th, y);
}

inline __m128 VECTORCALL _mm_sincos_ps(__m128* cv, __m128 x)
{
  __m128 xx = _mm_mul_ps(x, x);
  __m128 t  = _mm_sub_ps(_mm_set1_ps(1.0f), _mm_mul_ps(xx, _mm_set1_ps(0.5f))); 
  xx = _mm_mul_ps(x, x);
  __m128 st = _mm_sub_ps(x, _mm_mul_ps(xx, _mm_set1_ps(0.16666666666f))); 
  
  xx = _mm_mul_ps(x, x);
  t  = _mm_add_ps(t, _mm_mul_ps(xx, _mm_set1_ps(0.04166666666f)));
  xx = _mm_mul_ps(x, x);
  st = _mm_add_ps(st, _mm_mul_ps(xx, _mm_set1_ps(0.00833333333f)));
  
  xx = _mm_mul_ps(x, x);
  t  = _mm_sub_ps(t, _mm_mul_ps(xx, _mm_set1_ps(0.00138888888f)));
  xx = _mm_mul_ps(x, x);
  st = _mm_sub_ps(st, _mm_mul_ps(xx, _mm_set1_ps(0.00019841269f)));
  
  xx = _mm_mul_ps(x, x);
  t  = _mm_add_ps(t, _mm_mul_ps(xx, _mm_set1_ps(2.48016e-05f)));
  xx = _mm_mul_ps(x, x);
  st = _mm_add_ps(st, _mm_mul_ps(xx, _mm_set1_ps(2.75573e-06f)));
  *cv = t;
  return st;
}

#endif //__clang__ || __gnu

#else // AX_SUPPORT_SSE

#define VecSwizzle1(vec, x)         MakeVec4(vec[x], vec[x], vec[x], vec[x])
#define VecSwizzle(vec, x, y, z, w) MakeVec4(vec[x], vec[y], vec[z], vec[w])
#define VecShuffle_0101(a, b)       MakeVec4(a.x, a.y, b.x, b.y)  // 00, 01, 10, 11
#define VecShuffle_2323(a, b)       MakeVec4(a.z, a.w, b.z, b.w)  // 02, 03, 12, 13
// return (vec1[x], vec1[y], vec2[z], vec2[w])
#define VecShuffle(vec1, vec2, x, y, z, w) MakeVec4(vec1[x], vec1[y], vec2[z], vec2[w])

#endif
#ifdef AX_SUPPORT_AVX2

__forceinline __m256i VECTORCALL SSESelect(const __m256i V1, const __m256i V2, const __m256i& Control)
{
    return _mm256_blendv_epi8(V1, V2, Control);
}

// https://stackoverflow.com/questions/17863411/sse-multiplication-of-2-64-bit-integers
__forceinline __m256i VECTORCALL Multiply64Bit(__m256i ab, __m256i cd)
{
	__m256i b = _mm256_srli_epi64(ab, 64);
	__m256i bc = _mm256_mul_epu32(b, cd);
	__m256i d = _mm256_srli_epi64(cd, 64);
	__m256i ad = _mm256_mul_epu32(ab, d);
	__m256i high = _mm256_add_epi64(bc, ad);
	high = _mm256_slli_epi64(high, 64);
	return _mm256_add_epi64(high, _mm256_mul_epu32(ab, cd));
}

__forceinline int VECTORCALL hsum_256_epi32(__m256i v)
{
	__m128i sum128 = _mm_add_epi32(_mm256_castsi256_si128(v), _mm256_extracti128_si256(v, 1));
	return hsum_128_epi32avx(sum128);
}

__forceinline int64 VECTORCALL hsum_256_epi64(__m256i v)
{
	return _mm256_cvtsi256_si32(v) + _mm256_extract_epi64(v, 1) + _mm256_extract_epi64(v, 2) + _mm256_extract_epi64(v, 3);
}

__forceinline double VECTORCALL hsum_256_pd(__m256d v)
{
	__m128d sum128 = _mm_add_pd(_mm256_castpd256_pd128(v), _mm256_extractf128_pd(v, 1));
	return hsum_128_pdavx(sum128);
}

// from: Faster Population Counts Using AVX2 Instructions resource paper
__forceinline int64 VECTORCALL popcount256_epi64(__m256i v)
{
	const __m256i lookup = _mm256_setr_epi8(0, 1, 1, 2, 1, 2, 2, 3, 1, 2,
		2, 3, 2, 3, 3, 4, 0, 1, 1, 2, 1, 2, 2, 3,
		1, 2, 2, 3, 2, 3, 3, 4);
	const __m256i low_mask = _mm256_set1_epi8(0x0f);
	__m256i lo =  _mm256_and_si256(v, low_mask);
	__m256i hi = _mm256_and_si256(_mm256_srli_epi32(v, 4), low_mask);
	__m256i popcnt1 = _mm256_shuffle_epi8(lookup, lo);
	__m256i popcnt2 = _mm256_shuffle_epi8(lookup, hi);
	__m256i total = _mm256_add_epi8(popcnt1, popcnt2);
	v = _mm256_sad_epu8(total, _mm256_setzero_si256());
	return _mm256_cvtsi256_si32(v) + _mm256_extract_epi64(v, 1) + _mm256_extract_epi64(v, 2) + _mm256_extract_epi64(v, 3);
}

__forceinline __m256i VECTORCALL popcnt256si(__m256i v) // returns 4 64 bit integer that contains pop counts
{
	const __m256i lookup = _mm256_setr_epi8(0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4);
	const __m256i low_mask = _mm256_set1_epi8(0x0f);
	__m256i lo = _mm256_and_si256(v, low_mask);
	__m256i hi = _mm256_and_si256(_mm256_srli_epi32(v, 4), low_mask);
	__m256i popcnt1 = _mm256_shuffle_epi8(lookup, lo);
	__m256i popcnt2 = _mm256_shuffle_epi8(lookup, hi);
	return _mm256_sad_epu8(_mm256_add_epi8(popcnt1, popcnt2), _mm256_setzero_si256());
}

__forceinline float VECTORCALL hsum256_ps_avx(__m256 v) {
	__m128 vlow = _mm256_castps256_ps128(v);
	__m128 vhigh = _mm256_extractf128_ps(v, 1); // high 128
	vlow = _mm_add_ps(vlow, vhigh);     // add the low 128
	return hsum_ps_sse3(vlow);         // and inline the sse3 version, which is optimal for AVX
}
#endif // AX_SUPPORT_AVX2

AX_END_NAMESPACE 