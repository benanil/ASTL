
/******************************************************************************************
*  Purpose:                                                                               *
*    Conversion of 16 bit floating point values and 32 bit floating point values          *
*    SSE and AVX used for performance but not required, scalar versions exists            *
*  Good To Know:                                                                          *
*    ...                                                                                  *
*  Author:                                                                                *
*    Anilcan Gulkaya 2024 anilcangulkaya7@gmail.com github @benanil                       *
*******************************************************************************************/

// todo better check for arm neon fp16

#pragma once

#include "SIMDVectorMath.hpp"

/*//////////////////////////////////////////////////////////////////////////*/
/*                             Half                                         */
/*//////////////////////////////////////////////////////////////////////////*/

#ifdef AX_SUPPORT_NEON
typedef float16_t half;
#else
typedef ushort half;
constexpr half OneFP16 = 15360;
constexpr half MinusOneFP16 = 48128;
constexpr half ZeroFP16 = 0;
constexpr half HalfFP16 = 14336; // fp16 0.5
constexpr half Sqrt2FP16 = 15784; // fp16 sqrt(2)

typedef uint half2;
constexpr half2 Half2Up    = OneFP16 << 16u;
constexpr half2 Half2Down  = MinusOneFP16 << 16u;
constexpr half2 Half2Left  = MinusOneFP16;
constexpr half2 Half2Right = OneFP16;
constexpr half2 Half2One   = OneFP16 | (OneFP16 << 16);
constexpr half2 Half2Zero  = 0;

#define MakeHalf2(x, y) ((x) | ((y) << 16))
#define Half2SetX(v, x) (v &= 0xFFFF0000u, v |= x;)
#define Half2SetY(v, y) (v &= 0x0000FFFFu, v |= y;)
#endif

typedef uint half2;

// todo better check for half support
purefn float ConvertHalfToFloat(half x) 
{
#if defined(AX_SUPPORT_AVX2) 
    return _mm_cvtss_f32(_mm_cvtph_ps(_mm_set1_epi16(x))); 
#elif defined(AX_SUPPORT_NEON)
    return vgetq_lane_f32(vcvt_f32_f16(vdup_n_f16(x)), 0);
#else
    uint h = x;
    uint h_e = h & 0x00007c00u;
    uint h_m = h & 0x000003ffu;
    uint h_s = h & 0x00008000u;
    uint h_e_f_bias = h_e + 0x0001c000u;

    uint f_s  = h_s        << 0x00000010u;
    uint f_e  = h_e_f_bias << 0x0000000du;
    uint f_m  = h_m        << 0x0000000du;
    uint f_result = f_s | f_e | f_m;
        
    return BitCast<float>(f_result);
#endif
}

purefn half ConvertFloatToHalf(float Value) 
{
#if defined(AX_SUPPORT_AVX2)
    return _mm_extract_epi16(_mm_cvtps_ph(_mm_set_ss(Value), 0), 0);
#elif defined(AX_SUPPORT_NEON)
    return vget_lane_f16(vcvt_f16_f32(vdupq_n_f32(Value)), 0);
#else
    uint32_t Result; // branch removed version of DirectxMath function
    uint32_t IValue = BitCast<uint32_t>(Value);
    uint32_t Sign = (IValue & 0x80000000u) >> 16U;
    IValue = IValue & 0x7FFFFFFFu;      // Hack off the sign
    // if (IValue > 0x47FFEFFFu) {
    //     return 0x7FFFU | Sign; // The number is too large to be represented as a half.  Saturate to infinity.
    // }
    uint32_t mask = 0u - (IValue < 0x38800000u);
    uint32_t b = IValue + 0xC8000000U;
    uint32_t a = (0x800000u | (IValue & 0x7FFFFFu)) >> (113u - (IValue >> 23u));
    
    IValue = (mask & a) | (~mask & b);
    Result = ((IValue + 0x0FFFu + ((IValue >> 13u) & 1u)) >> 13u) & 0x7FFFu; 
    return (half)(Result | Sign);
#endif
}

inline void ConvertHalf2ToFloat2(float* result, uint32_t h) 
{
#if defined(AX_SUPPORT_AVX2)
    _mm_storel_pi((__m64 *)result, _mm_cvtph_ps(_mm_set1_epi16(h)));
#elif defined(AX_SUPPORT_NEON)
    float16x4_t halfVec = vreinterpret_f16_u32(vdup_n_u32(h));
    vst1_f32(result, vget_low_f32(vcvt_f32_f16(halfVec)));
#else
    uint64_t h2 = (uint64_t)(h & 0x0000FFFFull) | (uint64_t(h & 0xFFFF0000ull) << 16ull);

    uint64_t h_e = h2 & 0x00007c0000007c00ull;
    uint64_t h_m = h2 & 0x000003ff000003ffull;
    uint64_t h_s = h2 & 0x0000800000008000ull;
    uint64_t h_e_f_bias = h_e + 0x0001c0000001c000ull;

    uint64_t f_s  = h_s        << 0x00000010ull;
    uint64_t f_e  = h_e_f_bias << 0x0000000dull;
    uint64_t f_m  = h_m        << 0x0000000dull;
    uint64_t f_result = f_s | f_e | f_m;
        
    result[0] = BitCast<float>((uint32_t)(f_result & 0xFFFFFFFFu));
    result[1] = BitCast<float>((uint32_t)(f_result >> 32ull));
#endif
}

inline uint32_t ConvertFloat2ToHalf2(const float* float2)
{
#if defined(AX_SUPPORT_NEON)
    float32x2_t x = vld1_dup_f32(float2);
    float32x4_t x4 = vcombine_f32(x, x);
    return vget_lane_u32(vreinterpret_u32_f16(vcvt_f16_f32(x4)), 0);
#elif defined(AX_SUPPORT_AVX)
    return _mm_extract_epi32(_mm_cvtps_ph(_mm_set_ss(float2), 0), 0);
#endif
    uint32_t result = 0;
    result  = ConvertFloatToHalf(float2[0]);
    result |= (uint32_t)ConvertFloatToHalf(float2[1]) << 16;
    return result;
}

// input half4 is 4x 16 bit integer for example it can be uint64_t
inline void ConvertHalf4ToFloat4(float* result, const half* half4)
{
#ifdef AX_SUPPORT_AVX2
    _mm_storeu_ps(result, _mm_cvtph_ps(_mm_loadu_si64(half4)));

#elif defined(AX_SUPPORT_NEON)
    vst1q_f32(result, vcvt_f32_f16(vld1_dup_f16(half4)));

#elif defined(AX_SUPPORT_SSE)
    Vector4x32u h4 = VeciLoad64((const uint64_t*)half4);
    h4 = VeciUnpackLow16(h4, VeciZero());   // [half4.xy, half4.xy, half4.zw, half4.zw] 
    
    Vector4x32u h_e = VeciAnd(h4, VeciSet1(0x00007c00));
    Vector4x32u h_m = VeciAnd(h4, VeciSet1(0x000003ff));
    Vector4x32u h_s = VeciAnd(h4, VeciSet1(0x00008000));
    Vector4x32u h_e_f_bias = VeciAdd(h_e, VeciSet1(0x0001c000));
    
    Vector4x32u f_s  = VeciSll32(h_s, 0x00000010);
    Vector4x32u f_e  = VeciSll32(h_e_f_bias, 0x0000000d);
    Vector4x32u f_m  = VeciSll32(h_m, 0x0000000d);
    Vector4x32u f_em = VeciOr(f_e, f_m);

    Vector4x32u i_result = VeciOr(f_s, f_em);
    VecStore(result, VeciToVecf(i_result));
    
#else // no intrinsics
    ConvertHalf2ToFloat2(result, *(uint32_t*)half4);
    ConvertHalf2ToFloat2(result + 2, *((uint32_t*)(half4) + 1));
#endif
}

inline void ConvertFloat4ToHalf4(half* result, const float* float4)
{
#ifdef AX_SUPPORT_AVX2

    *((long long*)result) = _mm_extract_epi64(_mm_cvtps_ph(_mm_loadu_ps(float4), _MM_FROUND_TO_NEAREST_INT), 0);

#elif defined(AX_SUPPORT_NEON)

    *(float16x4_t*)result = vcvt_f16_f32(vld1q_f32(float4));

#elif defined(AX_SUPPORT_SSE)

    Vector4x32u IValue = VeciLoad((const unsigned int*)float4);
    Vector4x32u Sign = VeciSrl32(VeciAnd(IValue, VeciSet1(0x80000000u)), 16);
    IValue = VeciAnd(IValue, VeciSet1(0x7FFFFFFFu));      // Hack off the sign
    
    Vector4x32u mask = VeciCmpLt(IValue, VeciSet1(0x38800000u));
    Vector4x32u b = VeciAdd(IValue, VeciSet1(0xC8000000u));
    Vector4x32u a = VeciOr(VeciSet1(0x800000u), VeciAnd(IValue, VeciSet1(0x7FFFFFu)));
    a = VeciSrl(a, VeciSub(VeciSet1(113u), VeciSrl32(IValue, 23u)));
    
    IValue = VeciBlend(b, a, mask);

    Vector4x32u Result = VeciAdd(IValue, VeciSet1(0x0FFFu));
    Result = VeciAdd(Result, VeciAnd(VeciSrl32(IValue, 13u), VeciSet1(1u)));
    Result = VeciSrl32(Result, 13u);
    Result = VeciAnd(Result, VeciSet1(0x7FFFu));
    Result = VeciOr(Result, Sign);

    #ifdef AX_SUPPORT_SSE
        const int shufleMask = MakeShuffleMask(0, 2, 1, 3);
        __m128i lo = _mm_shufflelo_epi16(Result, shufleMask);
        __m128i hi = _mm_shufflehi_epi16(lo, shufleMask);
        Result = _mm_shuffle_epi32(hi, shufleMask);
        *((long long*)result) = _mm_extract_epi64(Result, 0);
    #else
        // todo test
        // Narrow the 32-bit to 16-bit, effectively extracting the lower 16 bits of each element
        uint16x4_t low16_bits = vmovn_u32(Result);  // Narrow to 16 bits per element
        // Directly cast the `uint16x4_t` to `uint64_t`
        vst1_u64((uint64_t*)result, vreinterpret_u64_u16(low16_bits));
    #endif

#else // no intrinsics
    *(uint32_t*)result = ConvertFloat2ToHalf2(float4);
    *((uint32_t*)result + 1) = ConvertFloat2ToHalf2(float4 + 2);
#endif
}


#ifdef AX_SUPPORT_AVX2

// convert 8 float and half with one instruction
#define ConvertFloat8ToHalf8(result, float8) _mm_storeu_si128((__m128i*)result, _mm256_cvtps_ph(_mm256_loadu_ps(float8), _MM_FROUND_TO_NEAREST_INT))

#define ConvertHalf8ToFloat8(result, half8)  _mm256_storeu_ps(result,  _mm256_cvtph_ps(_mm_loadu_si128((const __m128i*)half8)))

#else

inline void ConvertHalf8ToFloat8(float* float8, const half* half8)
{
    ConvertHalf4ToFloat4(float8    , half8);
    ConvertHalf4ToFloat4(float8 + 4, half8 + 4);
}

inline void ConvertFloat8ToHalf8(half* result, const float* float8)
{
    ConvertFloat4ToHalf4(result    , float8);
    ConvertFloat4ToHalf4(result + 4, float8 + 4);
}

#endif // AX_SUPPORT_AVX2


inline void ConvertHalfToFloatN(float* res, const half* x, const int n) 
{   
    for (int i = 0; i < n; i += 8, x += 8, res += 8)
        ConvertHalf8ToFloat8(res, x);
 
    for (int i = 0; i < (n & 7); i++, res++, x++)
        *res = ConvertHalfToFloat(*x);
}

inline void ConvertFloatToHalfN(half* res, const float* x, const int n) 
{   
    for (int i = 0; i < n; i += 8, x += 8, res += 8)
        ConvertFloat8ToHalf8(res, x);
 
    for (int i = 0; i < (n & 7); i++, res++, x++)
        *res = ConvertFloatToHalf(*x);
}

purefn float3 ConvertHalf3ToFloat3(half* h) {
	float3 res; 
    ConvertHalf2ToFloat2(&res.x, *(uint32_t*)h); 
    res.z = ConvertHalfToFloat(h[2]); 
    return res;
}