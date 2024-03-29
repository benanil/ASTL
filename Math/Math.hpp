// most of the functions are accurate and faster than stl 
// convinient for game programming

#pragma once

#include "../Common.hpp" // includes MIN, MAX, Clamp, Abs and FAbs

AX_NAMESPACE 

// constants
constexpr float PI       = 3.14159265358f;
constexpr float RadToDeg = 180.0f / PI;
constexpr float DegToRad = PI / 180.0f;
constexpr float OneDivPI = 1.0f / PI;
constexpr float PIDiv2   = PI / 2.0f;
constexpr float TwoPI    = PI * 2.0f;
constexpr float Sqrt2    = 1.414213562f;
constexpr float Epsilon  = 0.0001f;
// for integer constants use stdint.h's INT32_MIN, INT64_MIN, INT32_MAX...
// for float constants use float.h's FLT_MAX, FLT_MIN, DBL_MAX, DBL_MIN

/*//////////////////////////////////////////////////////////////////////////*/
/*                           Base Functions                                 */
/*//////////////////////////////////////////////////////////////////////////*/

__forceinline __constexpr float Lerp(float x, float y, float t) {
	return x + (y - x) * t;
}

__forceinline __constexpr double Lerp(double x, double y, double t) {
	return x + (y - x) * t;
}

__forceinline __constexpr float SqrtConstexpr(float a) {
    // from: Jack W. Cerenshaw's math toolkit for real time development book: page 63 Listing 4.3
    // I've removed some of the branches. slightly slower than sqrtf
    // double version is also here: https://gist.github.com/benanil/9d1668c0befb24263e27bd04dfa2e90f#file-mathfunctionswithoutstl-c-L230
    const double A = 0.417319242, B = 0.5901788532;
    union { double fp; struct { unsigned lo, hi; }; } x {};
    if (a <= 0.001) return 0.0f;
    x.fp = (double)a;
    // grab the exponent
    unsigned expo = (x.hi >> 20u) - 0x3fe;
    x.hi &= 0x000fffffu;
    x.hi += 0x3fe00000u;
    // get square root of normalized number
    double root = A + B * x.fp;
    root = 0.5 * (x.fp / root + root); // you can even remove this haha if you do that you might want to reduce this: 0.414213562
    root = 0.5 * (x.fp / root + root);
    // root = 0.5 * (x.fp / root + root); // iterate 3 times probably overkill
    // now rebuild the result
    x.fp = root;
    bool isOdd = expo & 1;
    expo += isOdd;
    x.fp *= 1.0 + (isOdd * 0.414213562);
    expo = (expo + (expo < 0u)) >> 1;
    // put it back
    expo += 0x3feu;
    x.hi &= 0x000fffffu;
    x.hi += expo << 20u;
    return (float)x.fp;
}

__forceinline float Sqrt(float a) {
#ifdef AX_SUPPORT_SSE
	return _mm_cvtss_f32(_mm_sqrt_ps(_mm_set_ps1(a)));
#elif defined(__clang__)
	return __builtin_sqrt(a);
#else
	return SqrtConstexpr(a);
#endif
}

// original doom rsqrt implementation with comments :)
// https://en.wikipedia.org/wiki/Fast_inverse_square_root
// https://rrrola.wz.cz/inv_sqrt.html   <- fast and 2.5x accurate
__forceinline float RSqrt(float x) {
#ifdef AX_SUPPORT_SSE
	// when I compile with godbolt -O1 expands to one instruction vrsqrtss
	// which is approximately equal latency as mulps.
	// RSqrt(float):                             # @RSqrt3(float)
  //       vrsqrtss        xmm0, xmm0, xmm0
  //       ret
	return _mm_cvtss_f32(_mm_rsqrt_ps(_mm_set_ps1(x)));
#elif defined(AX_ARM)
	return vget_lane_f32(vrsqrte_f32(vdup_n_f32(x)), 0);
#else
	float f = x;
	uint32_t i = BitCast<uint32_t>(f);
	i = (uint32_t)(0x5F1FFFF9ul - (i >> 1));
	f = BitCast<float>(i);
	return 0.703952253f * f * (2.38924456f - x * f * f);
#endif
}

__forceinline __constexpr float Fract(float a) { a = Abs(a); return a - int(a); }

// https://github.com/id-Software/DOOM-3/blob/master/neo/idlib/math/Math.h
__forceinline float Exp(float f) {
    const int IEEE_FLT_MANTISSA_BITS  =	23;
    const int IEEE_FLT_EXPONENT_BITS  =	8;
    const int IEEE_FLT_EXPONENT_BIAS  =	127;
    const int IEEE_FLT_SIGN_BIT       =	31;
    int i, s, e, m, exponent;
    float x, x2, y, p, q;
    
    x = f * 1.44269504088896340f;		// multiply with ( 1 / log( 2 ) )
    
    i = BitCast<int>(x);
    s = ( i >> IEEE_FLT_SIGN_BIT );
    e = ( ( i >> IEEE_FLT_MANTISSA_BITS ) & ( ( 1 << IEEE_FLT_EXPONENT_BITS ) - 1 ) ) - IEEE_FLT_EXPONENT_BIAS;
    m = ( i & ( ( 1 << IEEE_FLT_MANTISSA_BITS ) - 1 ) ) | ( 1 << IEEE_FLT_MANTISSA_BITS );
    i = ( ( m >> ( IEEE_FLT_MANTISSA_BITS - e ) ) & ~( e >> 31 ) ) ^ s;
    
    exponent = ( i + IEEE_FLT_EXPONENT_BIAS ) << IEEE_FLT_MANTISSA_BITS;
    y = BitCast<float>(exponent);
    x -= (float) i;
    if ( x >= 0.5f ) {
        x -= 0.5f;
        y *= 1.4142135623730950488f;	// multiply with sqrt( 2 )
    }
    x2 = x * x;
    p = x * ( 7.2152891511493f + x2 * 0.0576900723731f );
    q = 20.8189237930062f + x2;
    x = y * ( q + p ) / ( q - p );
    return x;
}

// https://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/
// https://github.com/ekmett/approximate/blob/master/cbits/fast.c#L81
// should be much more precise with large b
__forceinline __constexpr float Pow(float a, float b) {
    // calculate approximation with fraction of the exponent
    int e = (int) b;
    union {
        double d;
        int x[2];
    } u = { (double)a };
    u.x[1] = (int)((b - e) * (u.x[1] - 1072632447) + 1072632447.0f);
    u.x[0] = 0;
    
    // exponentiation by squaring with the exponent's integer part
    // double r = u.d makes everything much slower, not sure why
    float r = 1.0;
    while (e) {
        if (e & 1) {
            r *= a;
        }
        a *= a;
        e >>= 1;
    }
    
    return (float)(u.d * r);
}

// https://github.com/ekmett/approximate/blob/master/cbits/fast.c#L81 <--you can find double versions
__forceinline __constexpr float Log(float x) {
    return (BitCast<int>(x) - 1064866805) * 8.262958405176314e-8f;
}

// if you want log10 for integer you can look at Algorithms.hpp
__forceinline __constexpr float Log10(float x) { 
    return Log(x) / 2.30258509299f; // ln(x) / ln(10)
} 

// you might look at this link as well: https://tech.ebayinc.com/engineering/fast-approximate-logarithms-part-i-the-basics/
__forceinline __constexpr float Log2(float x) {
    return Log(x) / 0.6931471805599453094f; // ln(x) / ln(2) 
}

// https://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
__forceinline __constexpr unsigned int Log2(unsigned int v)
{
    constexpr int MultiplyDeBruijnBitPosition[32] = 
    {
        0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
        8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31
    };
    
    // first round down to one less than a power of 2
    v |= v >> 1; v |= v >> 2; v |= v >> 4; v |= v >> 8; v |= v >> 16;
    return MultiplyDeBruijnBitPosition[(uint32_t)(v * 0x07C4ACDDU) >> 27];
}

inline __constexpr unsigned int Log10(unsigned int v)
{                                                   // this would work too
    unsigned int const PowersOf10[] = {             // if (n <= 9) return 1;
        1, 10, 100, 1000, 10000, 100000,            // if (n <= 99) return 2;
        1000000, 10000000, 100000000, 1000000000    // if (n <= 999) return 3;
    };                                              // ...
    unsigned int t = (Log2(v) + 1) * 1233 >> 12;    // if (n <= 2147483647) return 10; // 2147483647 = int max
    return t - (v < PowersOf10[t]);                                                      
}                                                                                        

__forceinline __constexpr bool IsZero(float x) {
    return Abs(x) <= 0.0001f; 
}

__forceinline __constexpr bool AlmostEqual(float x, float  y) {
    return Abs(x-y) <= 0.001f;
}

__forceinline __constexpr float Sign(float x) 
{
    int res = BitCast<int>(1.0f);
    res |= BitCast<int>(x) & 0x80000000;
    return BitCast<float>(res);
} 

__forceinline __constexpr float CopySign(float x, float y) {
    int ix = BitCast<int>(x) & 0x7fffffff;
    int iy = BitCast<int>(y) & 0x80000000;
    return BitCast<float>(ix | iy);
}

__forceinline __constexpr bool IsNan(float f) {
    uint32 intValue = BitCast<uint32>(f);
    uint32 exponent = (intValue >> 23) & 0xFF;
    uint32 fraction = intValue & 0x7FFFFF;
    return (exponent == 0xFF) && (fraction != 0);
}

template<typename T>
__forceinline __constexpr T FMod(T x, T y) {
    T quotient = x / y;
    T whole = (T)((int)quotient);  // truncate quotient to integer
    T remainder = x - whole * y;
    remainder += (remainder < (T)0.0) * y;
    return remainder;
}

template<typename T>
__forceinline __constexpr T Floor(T x) {
    T whole = (T)(int)x;  // truncate quotient to integer
    return x - (x-whole);
}

//  ######################################  
//  #####  [TRIGONOMETRIC FUNCTINS]  #####  
//  ######################################  

// https://mazzo.li/posts/vectorized-atan2.html
__forceinline __constexpr float ATan(float x) {
    const float x_sq = x * x;
    const float a1 =  0.99997726f, a3 = -0.33262347f, a5  = 0.19354346f,
                a7 = -0.11643287f, a9 =  0.05265332f, a11 = -0.01172120f;
    return x * (a1 + x_sq * (a3 + x_sq * (a5 + x_sq * (a7 + x_sq * (a9 + x_sq * a11)))));
}

__forceinline __constexpr float ATan2(float y, float x) {
    // from here: https://yal.cc/fast-atan2/  
    float ay = Abs(y), ax = Abs(x);
    int invert = ay > ax;
    float z = invert ? ax/ay : ay/ax;// [0,1]
    float th = ATan(z);              // [0,π/4]
    if (invert)   th = PIDiv2 - th;  // [0,π/2]
    if (x < 0.0f) th = PI - th;      // [0,π]
    return CopySign(th, y);          // [-π,π] 
}

__forceinline float ASin(float z) {
	return ATan2(z, Sqrt(1.0f-(z * z)));
}

// https://en.wikipedia.org/wiki/Sine_and_cosine
// warning: accepts input between -TwoPi and TwoPi  if (Abs(x) > TwoPi) use x = FMod(x + PI, TwoPI) - PI;

__forceinline __constexpr float RepeatPI(float x) {
  return FMod(x + PI, TwoPI) - PI;
}

__forceinline __constexpr float Sin(float x) {
    float xx = x * x * x;                // x^3
    float t = x - (xx * 0.16666666666f); // x3/!3  6 = !3 = 1.6666 = rcp(3)
    t += (xx *= x * x) * 0.00833333333f; // x5/!5  120 = !5 = 0.0083 = rcp(5)
    t -= (xx *= x * x) * 0.00019841269f; // x7/!7  5040 = !7
    t += (xx * x * x)  * 2.75573e-06f;   // 362880 = !9
    return t;
}

// warning: accepts input between -TwoPi and TwoPi  if (Abs(x) > TwoPi) use x = FMod(x + PI, TwoPI) - PI;
__forceinline __constexpr float Cos(float x) {
    float xx = x * x;                    // x^2
    float t = 1.0f - (xx * 0.5f);        // 1-(x2/!2) 0.5f = rcp(!2)
    t += (xx *= x * x) * 0.04166666666f; // t + (x4/!4) 0.04166 = rcp(24) = 24.0f = !4
    t -= (xx *= x * x) * 0.00138888888f; // t - (x6/!6) 720.0f = !6
    t += (xx * x * x)  * 2.48016e-05f;   // t + (x8/!8) 40320 = !8
    return t;
}

__forceinline __constexpr void SinCos(float x, float* s, float* c) {
	*s = Sin(x); *c = Cos(x); 
}

__forceinline __constexpr float Tan(float radians) {
    float rr = radians * radians;
    float a = 9.5168091e-03f;
    a *= rr; a += 2.900525e-03f;
    a *= rr; a += 2.45650893e-02f;
    a *= rr; a += 5.33740603e-02f;
    a *= rr; a += 1.333923995e-01f;
    a *= rr; a += 3.333314036e-01f;
    a *= rr; a += 1.0f;
    return a * radians;
}

// inspired from Casey Muratori's performance aware programming
// this functions makes the code more readable. OpenCL and Cuda has the same functions as well
__forceinline float ATan2PI(float y, float x) { return ATan2(y, x) / PI; }
__forceinline float ASinPI(float z) { return ASin(z) / PI; }
__forceinline float ACos(float x)   { return PIDiv2 - ASin(x); }
__forceinline float ACosPI(float x) { return ACos(x) / PI; }
__forceinline __constexpr float CosPI(float x)  { return Cos(x) / PI; }
__forceinline __constexpr float SinPI(float x)  { return Sin(x) / PI; }

//  ######################################  
//  #####      [HALF FUNCTINS]       #####  
//  ######################################  

typedef ushort half;

__forceinline float ConvertHalfToFloat(half x) {
#if defined(AX_SUPPORT_SSE) 
	return _mm_cvtss_f32(_mm_cvtph_ps(_mm_set1_epi16(x))); 
#elif defined(__ARM_NEON__)
	return _cvtsh_ss(x); 
#else
    uint Mantissa, Exponent, Result;
    Mantissa = (uint)(x & 0x03FF);
    
    if ((x & 0x7C00) != 0)  // The value is normalized
    {
        Exponent = (uint)((x >> 10) & 0x1F);
    }
    else if (Mantissa != 0) // The value is denormalized
    {
        // Normalize the value in the resulting float
        Exponent = 1;
        do {
            Exponent--;
            Mantissa <<= 1;
        } while ((Mantissa & 0x0400) == 0);
    
        Mantissa &= 0x03FF;
    }
    else { // The value is zero
        Exponent = (uint)-112;
    }
    
    Result = ((x & 0x8000) << 16) | ((Exponent + 112) << 23) | (Mantissa << 13);
    return BitCast<float>(Result);
#endif
}
// faster but not always safe version. taken from: https://stackoverflow.com/questions/1659440/32-bit-to-16-bit-floating-point-conversion
// const uint e = (x & 0x7C00) >> 10; // exponent
// const uint m = (x & 0x03FF) << 13; // mantissa
// const uint v = BitCast<uint>((float)m) >> 23; // evil log2 bit hack to count leading zeros in denormalized format
// uint a = (x & 0x8000) << 16 | (e != 0) * ((e + 112) << 23 | m);
// a |= ((e == 0) & (m != 0)) * ((v - 37) << 23 | ((m << (150 - v)) & 0x007FE000));
// return BitCast<float>(a); // sign : normalized : denormalized

// converts maximum 4 float
__forceinline void ConvertHalfToFloat(float* res, const half* x, short n) {
#if defined(AX_SUPPORT_SSE)
	alignas(16) float a[4];
    half b[8]; 
	SmallMemCpy(b, x, n * sizeof(half));
	_mm_store_ps(a, _mm_cvtph_ps(_mm_loadu_si16(b))); // MSVC does not have scalar instructions.
	SmallMemCpy(res, a, n * sizeof(float));
#else
	for (int i = 0; i < n; i++)
		res[i] = ConvertHalfToFloat(x[i]);
#endif
}

__forceinline half ConvertFloatToHalf(float Value) {
#if defined(AX_SUPPORT_SSE)
	return _mm_extract_epi16(_mm_cvtps_ph(_mm_set_ss(Value), 0), 0);
#elif defined(__ARM_NEON__)
	return _cvtss_sh(Value, 0);
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

// converts maximum 4 half
__forceinline void ConvertFloatToHalf(half* res, const float* x, short n) {
#if defined(AX_SUPPORT_SSE)
    alignas(16) half a[8];
    float b[8]; 
    SmallMemCpy(b, x, n * sizeof(float));
    _mm_store_si128((__m128i*)a, _mm_cvtps_ph(_mm_loadu_ps(x), 0)); // MSVC does not have scalar instructions.
    SmallMemCpy(res, a, n * sizeof(half));
#else
    for (int i = 0; i < n; i++)
        res[i] = ConvertFloatToHalf(x[i]);
#endif
}

__forceinline void ConvertFloatToHalf4(half* res, const float* x)
{
#if defined(AX_SUPPORT_SSE) 
    _mm_store_si128((__m128i*)res, _mm_cvtps_ph(_mm_loadu_ps(x), _MM_FROUND_TO_NEAREST_INT)); 
#else
    ConvertFloatToHalf(res, x, 4);
#endif
}

// packs -1,1 range float to short
__forceinline short PackSnorm16(float x) {
    return (short)Clamp(x * (float)INT16_MAX, (float)INT16_MIN, (float)INT16_MAX);
}

__forceinline float UnpackSnorm16(short x) {
    return (float)x / (float)INT16_MAX;
}

// packs 0,1 range float to short
__forceinline short PackUnorm16(float x) {
    return (short)Clamp(x * (float)INT16_MAX, (float)INT16_MIN, (float)INT16_MAX);
}

__forceinline float UnpackUnorm16(short x) {
    return (float)x / (float)UINT16_MAX;
}

//  ######################################  
//  #####      [COLOR FUNCTINS]      #####  
//  ######################################  

__forceinline uint MakeColorRGBAU32(uint8 r, uint8 g, uint8 b, uint8 a) {
    return r | (uint(g) << 8) | (uint(b) << 16) | (uint(a) << 24);
}

__forceinline uint PackColorRGBU32(float r, float g, float b) {
    return (uint)(r * 255.0f) | ((uint)(g * 255.0f) << 8) | ((uint)(b * 255.0f) << 16);
}

__forceinline uint PackColorRGBU32(float* c) {
    return (uint)(*c * 255.0f) | ((uint)(c[1] * 255.0f) << 8) | ((uint)(c[2] * 255.0f) << 16);
}

__forceinline uint PackColorRGBAU32(float* c) {
    return (uint)(*c * 255.0f) | ((uint)(c[1] * 255.0f) << 8) | ((uint)(c[2] * 255.0f) << 16) | ((uint)(c[3] * 255.0f) << 24);
}

__forceinline void UnpackColorRGBf(unsigned color, float* colorf) {
    const float toFloat = 1.0f / 255.0f;
    colorf[0] = float(color >> 0  & 0xFF) * toFloat;
    colorf[1] = float(color >> 8  & 0xFF) * toFloat;
    colorf[2] = float(color >> 16 & 0xFF) * toFloat;
}

__forceinline void UnpackColorRGBAf(unsigned color, float* colorf) {
    const float toFloat = 1.0f / 255.0f;
    colorf[0] = float(color >> 0  & 0xFF) * toFloat;
    colorf[1] = float(color >> 8  & 0xFF) * toFloat;
    colorf[2] = float(color >> 16 & 0xFF) * toFloat;
    colorf[3] = float(color >> 24) * toFloat;
}

template<typename T> __forceinline __constexpr T Min3(T a, T b, T c) {
    T res = a < b ? a : b;
    return res < c ? c : res;
}

template<typename T> __forceinline __constexpr T Max3(T a, T b, T c) {
    T res = a > b ? a : b;
    return res > c ? res : c;
}

AX_END_NAMESPACE 