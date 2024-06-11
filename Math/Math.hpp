
// most of the functions are accurate and faster than stl 
// convinient for game programming, be aware of speed and preciseness tradeoffs because cstd has more accurate functions

// Sections of this file are:
// Essential    : square, sqrt, exp, pow...
// Trigonometry : sin, cos, tan, atan, atan2...
// Half         : IEEE 16bit float, conversion functions
// Color        : packing and unpacking rgba8 color
// Ease         : easeIn, easeOut...

#pragma once

#include "../Common.hpp" // includes MIN, MAX, Clamp, Abs and FAbs

AX_NAMESPACE 

// constants
constexpr float PI       = 3.14159265358f;
constexpr float RadToDeg = 180.0f / PI;
constexpr float DegToRad = PI / 180.0f;
constexpr float OneDivPI = 1.0f / PI;
constexpr float HalfPI   = PI / 2.0f;
constexpr float TwoPI    = PI * 2.0f;
constexpr float Sqrt2    = 1.414213562f;
constexpr float Epsilon  = 0.0001f;
// for integer constants use stdint.h's INT32_MIN, INT64_MIN, INT32_MAX...
// for float constants use float.h's FLT_MAX, FLT_MIN, DBL_MAX, DBL_MIN

/*//////////////////////////////////////////////////////////////////////////*/
/*                              Essential                                   */
/*//////////////////////////////////////////////////////////////////////////*/

inline_constexpr float Sqr(float x) {
    return x * x;
}

inline_constexpr float Lerp(float x, float y, float t) {
	return x + (y - x) * t;
}

inline_constexpr double Lerp(double x, double y, double t) {
	return x + (y - x) * t;
}

inline_constexpr double SqrtConstexpr(double a) {
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
    return x.fp;
}

inline_constexpr float SqrtConstexpr(float a) {
    return (float)SqrtConstexpr((double)a);
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

// https://github.com/id-Software/DOOM-3/blob/master/neo/idlib/math/Math.h
inline_constexpr float Exp(float f) {
    const int IEEE_FLT_MANTISSA_BITS  =	23;
    const int IEEE_FLT_EXPONENT_BITS  =	8;
    const int IEEE_FLT_EXPONENT_BIAS  =	127;
    const int IEEE_FLT_SIGN_BIT       =	31;
    
    float x = f * 1.44269504088896340f;		// multiply with ( 1 / log( 2 ) )
    
    int i = BitCast<int>(x);
    int s = ( i >> IEEE_FLT_SIGN_BIT );
    int e = ( ( i >> IEEE_FLT_MANTISSA_BITS ) & ( ( 1 << IEEE_FLT_EXPONENT_BITS ) - 1 ) ) - IEEE_FLT_EXPONENT_BIAS;
    int m = ( i & ( ( 1 << IEEE_FLT_MANTISSA_BITS ) - 1 ) ) | ( 1 << IEEE_FLT_MANTISSA_BITS );
    i = ( ( m >> ( IEEE_FLT_MANTISSA_BITS - e ) ) & ~( e >> 31 ) ) ^ s;
    
    int exponent = ( i + IEEE_FLT_EXPONENT_BIAS ) << IEEE_FLT_MANTISSA_BITS;
    float y = BitCast<float>(exponent);
    x -= (float) i;
    if ( x >= 0.5f ) {
        x -= 0.5f;
        y *= 1.4142135623730950488f;	// multiply with sqrt( 2 )
    }
    float x2 = x * x;
    float p = x * ( 7.2152891511493f + x2 * 0.0576900723731f );
    float q = 20.8189237930062f + x2;
    x = y * ( q + p ) / ( q - p );
    return x;
}

// https://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/
// https://github.com/ekmett/approximate/blob/master/cbits/fast.c#L81
// should be much more precise with large b
inline_constexpr float Pow(float a, float b) {
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
inline_constexpr float Log(float x) {
    return (BitCast<int>(x) - 1064866805) * 8.262958405176314e-8f;
}

// if you want log10 for integer you can look at Algorithms.hpp
inline_constexpr float Log10(float x) { 
    return Log(x) / 2.30258509299f; // ln(x) / ln(10)
} 

// you might look at this link as well: https://tech.ebayinc.com/engineering/fast-approximate-logarithms-part-i-the-basics/
inline_constexpr float Log2(float x) {
    return Log(x) / 0.6931471805599453094f; // ln(x) / ln(2) 
}

// https://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
inline_constexpr unsigned int Log2(unsigned int v) {
    constexpr int MultiplyDeBruijnBitPosition[32] = 
    {
        0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
        8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31
    };
    
    // first round down to one less than a power of 2
    v |= v >> 1; v |= v >> 2; v |= v >> 4; v |= v >> 8; v |= v >> 16;
    return MultiplyDeBruijnBitPosition[(uint32_t)(v * 0x07C4ACDDU) >> 27];
}

inline_constexpr unsigned int Log10(unsigned int v) {
                                                    // this would work too
    unsigned int const PowersOf10[] = {             // if (n <= 9) return 1;
        1, 10, 100, 1000, 10000, 100000,            // if (n <= 99) return 2;
        1000000, 10000000, 100000000, 1000000000    // if (n <= 999) return 3;
    };                                              // ...
    unsigned int t = (Log2(v) + 1) * 1233 >> 12;    // if (n <= 2147483647) return 10; // 2147483647 = int max
    return t - (v < PowersOf10[t]);                                                      
}                                                                                        

inline_constexpr bool IsZero(float x) {
    return Abs(x) <= 0.0001f; 
}

inline_constexpr bool AlmostEqual(float x, float  y) {
    return Abs(x-y) <= 0.001f;
}

inline_constexpr float Sign(float x) {
    int res = BitCast<int>(1.0f);
    res |= BitCast<int>(x) & 0x80000000;
    return BitCast<float>(res);
} 

inline_constexpr int Sign(int x) {
    return x < 0 ? -1 : 1; // equal to above float version
} 

inline_constexpr float CopySign(float x, float y) {
    int ix = BitCast<int>(x) & 0x7fffffff;
    int iy = BitCast<int>(y) & 0x80000000;
    return BitCast<float>(ix | iy);
}

inline_constexpr bool IsNan(float f) {
    uint32 intValue = BitCast<uint32>(f);
    uint32 exponent = (intValue >> 23) & 0xFF;
    uint32 fraction = intValue & 0x7FFFFF;
    return (exponent == 0xFF) && (fraction != 0);
}

template<typename T>
inline_constexpr T FMod(T x, T y) {
    T quotient = x / y;
    T whole = (T)((int)quotient);  // truncate quotient to integer
    T remainder = x - whole * y;
    remainder += (remainder < (T)0.0) * y;
    return remainder;
}

template<typename T>
inline_constexpr T Floor(T x) {
    T whole = (T)(int)x;  // truncate quotient to integer
    return x - (x-whole);
}

template<typename T>
inline_constexpr T Ceil(T x) {
    T whole = (T)(int)x;  // truncate quotient to integer
    return whole + float(x > whole);
}

inline_constexpr float Fract(float a) {
    return a - int(a); 
}

/*//////////////////////////////////////////////////////////////////////////*/
/*                      Trigonometric Functions                             */
/*//////////////////////////////////////////////////////////////////////////*/

// https://mazzo.li/posts/vectorized-atan2.html
inline_constexpr float ATan(float x) {
    const float x_sq = x * x;
    const float a1 =  0.99997726f, a3 = -0.33262347f, a5  = 0.19354346f,
                a7 = -0.11643287f, a9 =  0.05265332f, a11 = -0.01172120f;
    return x * (a1 + x_sq * (a3 + x_sq * (a5 + x_sq * (a7 + x_sq * (a9 + x_sq * a11)))));
}

// is not working negative angles
inline_constexpr float ATan2(float y, float x) {
    // from here: https://yal.cc/fast-atan2/  
    float ay = Abs(y), ax = Abs(x);
    int invert = ay > ax;
    float z = invert ? ax/ay : ay/ax;// [ 0 , 1]
    float th = ATan(z);              // [ 0 , pi/4]
    if (invert)   th = HalfPI - th;  // [ 0 , pi/2]
    if (x < 0.0f) th = PI - th;      // [ 0 , pi]
    return CopySign(th, y);          // [-pi, pi] 
}

inline_constexpr float ASin(float z) {
	return ATan2(z, SqrtConstexpr(1.0f-(z * z)));
}

// https://en.wikipedia.org/wiki/Sine_and_cosine
// warning: accepts input between -TwoPi and TwoPi  if (Abs(x) > TwoPi) use x = FMod(x + PI, TwoPI) - PI;

inline_constexpr float RepeatPI(float x) {
  return FMod(x + PI, TwoPI) - PI;
}

// Accepts input between -TwoPi and TwoPi, use SinR if value is bigger than this range
// https://x.com/nthnblair/status/1790836310531559701
// https://www.desmos.com/calculator/owkgky7wh4
// https://github.com/LancePutnam/Gamma/blob/master/Gamma/scl.h  :function sinfast
inline_constexpr float Sin(float x) 
{
    union fi32 { float f; int i; } u = {x};
    
    int lz = 0, gtpi = 0;
    lz = u.i >> 31; // get less than zero
    u.i &= 0x7fffffff; // make positive
    gtpi = u.f > PI;
    
    u.f -= gtpi * PI; // subtract pi if greater than pi
    u.f *= 0.63655f; // constant founded using desmos
    u.f *= 2.0f - u.f;
    u.f *= 0.225f * u.f + 0.775f; 
    
    u.i |= gtpi << 31; // negate if was bigger than PI
    u.i ^= lz << 31; // flip sign if was negative
    return u.f; 
}

// Accepts input between -TwoPi and TwoPi, use CosR if value is bigger than this range  
inline_constexpr float Cos(float a)
{
    int lz = 0, greater = 0;
    lz = a < 0.0f;
    a *= -2.0f * lz + 1.0f; // make positive
    greater = a > PI;
    
    a -= PI * greater;
    a *= 0.159f;
    a = 1.0f - 32.0f * a * a * (0.75f - a);
    return greater ? -a : a; // sqrt(1.0f - (Sin(a)*Sin(a))); // for better approximation
}

// calculates sin(x) between [0,pi]
inline_constexpr float Sin0pi(float x) {
    x *= 0.63655f; // constant founded using desmos
    x *= 2.0f - x;
    return x * (0.225f * x + 0.775f); 
}

// calculates cos(x) between [0,pi]
inline_constexpr float Cos0pi(float a) {
    a *= 0.159f;
    return 1.0f - 32.0f * a * a * (0.75f - a);
}

// R suffix allows us to use with greater range than -TwoPI, TwoPI
inline_constexpr float SinR(float x) {
    return Sin(FMod(x + PI, TwoPI) - PI);
}

// R suffix allows us to use with greater range than -TwoPI, TwoPI
inline_constexpr float CosR(float x) {
    return Cos(FMod(x + PI, TwoPI) - PI);
}

inline void SinCos(float x, float* sp, float* cp) 
{
    *sp = Sin(x);
    *cp = Cos(x);
}

inline_constexpr float Tan(float radians) {
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
inline_constexpr float ATan2PI(float y, float x) { return ATan2(y, x) / PI; }
inline_constexpr float ASinPI(float z) { return ASin(z) / PI; }
inline_constexpr float ACos(float x)   { return HalfPI - ASin(x); }
inline_constexpr float ACosPI(float x) { return ACos(x) / PI; }
inline_constexpr float CosPI(float x)  { return Cos(x) / PI; }
inline_constexpr float SinPI(float x)  { return Sin(x) / PI; }

/*//////////////////////////////////////////////////////////////////////////*/
/*                             Half                                         */
/*//////////////////////////////////////////////////////////////////////////*/

typedef ushort half;
constexpr half OneFP16 = 15360;
constexpr half MinusOneFP16 = 48128;
constexpr half ZeroFP16 = 0;
constexpr half HalfFP16 = 14336; // fp16 0.5
constexpr half Sqrt2FP16 = 15784; // fp16 sqrt(2)
#define HALF2XY(x, y) ((x) | ((y) << 16));

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

// maybe write 4 quantity version that is faster 
__forceinline void ConvertHalfToFloat(float* res, const half* x, short n) 
{
	for (int i = 0; i < n; i++)
		res[i] = ConvertHalfToFloat(x[i]);
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

__forceinline void ConvertFloatToHalf4(half* res, const float* x) {
#if defined(AX_SUPPORT_SSE) 
    _mm_store_si128((__m128i*)res, _mm_cvtps_ph(_mm_loadu_ps(x), _MM_FROUND_TO_NEAREST_INT)); 
#else
    ConvertFloatToHalf(res, x, 4);
#endif
}

// packs -1,1 range float to short
inline_constexpr short PackSnorm16(float x) {
    return (short)Clamp(x * (float)INT16_MAX, (float)INT16_MIN, (float)INT16_MAX);
}

inline_constexpr float UnpackSnorm16(short x) {
    return (float)x / (float)INT16_MAX;
}

// packs 0,1 range float to short
inline_constexpr short PackUnorm16(float x) {
    return (short)Clamp(x * (float)INT16_MAX, (float)INT16_MIN, (float)INT16_MAX);
}

inline_constexpr float UnpackUnorm16(short x) {
    return (float)x / (float)UINT16_MAX;
}

/*//////////////////////////////////////////////////////////////////////////*/
/*                            Color                                         */
/*//////////////////////////////////////////////////////////////////////////*/

inline_constexpr uint MakeColorRGBAU32(uint8 r, uint8 g, uint8 b, uint8 a) {
    return r | (uint(g) << 8) | (uint(b) << 16) | (uint(a) << 24);
}

inline_constexpr uint PackColorRGBU32(float r, float g, float b) {
    return (uint)(r * 255.0f) | ((uint)(g * 255.0f) << 8) | ((uint)(b * 255.0f) << 16);
}

inline_constexpr uint PackColorRGBU32(float* c) {
    return (uint)(*c * 255.0f) | ((uint)(c[1] * 255.0f) << 8) | ((uint)(c[2] * 255.0f) << 16);
}

inline_constexpr uint PackColorRGBAU32(float* c) {
    return (uint)(*c * 255.0f) | ((uint)(c[1] * 255.0f) << 8) | ((uint)(c[2] * 255.0f) << 16) | ((uint)(c[3] * 255.0f) << 24);
}

inline_constexpr void UnpackColorRGBf(unsigned color, float* colorf) {
    const float toFloat = 1.0f / 255.0f;
    colorf[0] = float(color >> 0  & 0xFF) * toFloat;
    colorf[1] = float(color >> 8  & 0xFF) * toFloat;
    colorf[2] = float(color >> 16 & 0xFF) * toFloat;
}

inline_constexpr void UnpackColorRGBAf(unsigned color, float* colorf) {
    const float toFloat = 1.0f / 255.0f;
    colorf[0] = float(color >> 0  & 0xFF) * toFloat;
    colorf[1] = float(color >> 8  & 0xFF) * toFloat;
    colorf[2] = float(color >> 16 & 0xFF) * toFloat;
    colorf[3] = float(color >> 24) * toFloat;
}

/*//////////////////////////////////////////////////////////////////////////*/
/*                                 Easing                                   */
/*//////////////////////////////////////////////////////////////////////////*/

// to see visually: https://easings.net/ 
inline_constexpr float EaseIn(float x) {
    return x * x;
}

inline_constexpr float EaseOut(float x) { 
    float r = 1.0f - x;
    return 1.0f - (r * r); 
}

inline_constexpr float EaseInOut(float x) {
    return x < 0.5f ? 2.0f * x * x : 1.0f - Sqr(-2.0f * x + 2.0f) / 2.0f;
}

// integral symbol shaped interpolation, similar to EaseInOut
inline_constexpr float SmoothStep(float x) {
    return x * x * (3.0f - x * 2.0f);
}

inline_constexpr float EaseInSine(float x) {
  return 1.0f - Cos((x * PI) * 0.5f);
}

inline_constexpr float EaseOutSine(float x) {
  return Sin((x * PI) * 0.5f);
}

// Gradually changes a value towards a desired goal over time.
inline float SmoothDamp(float current, float target, float& currentVelocity, float smoothTime, float maxSpeed, float deltaTime)
{
    // Based on Game Programming Gems 4 Chapter 1.10
    smoothTime = MAX(0.0001f, smoothTime);
    float omega = 2.0f / smoothTime;

    float x = omega * deltaTime;
    float exp = 1.0f / (1.0f + x + 0.48f * x * x + 0.235f * x * x * x);
    float change = current - target;
    float originalTo = target;

    // Clamp maximum speed
    float maxChange = maxSpeed * smoothTime;
    change = Clamp(change, -maxChange, maxChange);
    target = current - change;

    float temp = (currentVelocity + omega * change) * deltaTime;
    currentVelocity = (currentVelocity - omega * temp) * exp;
    float output = target + (change + temp) * exp;

    // Prevent overshooting
    if (originalTo - current > 0.0f == output > originalTo)
    {
        output = originalTo;
        currentVelocity = (output - originalTo) / deltaTime;
    }

    return output;
}

inline float Remap(float in, float inMin, float inMax, float outMin, float outMax)
{
    return outMin + (in - inMin) * (outMax - outMin) / (inMax - inMin);
}

inline float Repeat(float t, float length)
{
    return Clamp(t - Floor(t / length) * length, 0.0f, length);
}

__forceinline float Step(float edge, float x)
{
    return float(x > edge);
}

// https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
// position(x0, y0), linePos1(x1, y1), linePos2(x2, y2) 
inline float LineDistance(float x0, float y0, float x1, float y1, float x2, float y2)
{
    float a = ((x2 - x1) * (y0 - y1)) - ((x0 - x1) * (y2 - y1));
    return Abs(a) / Sqrt(Sqr(x2 - x1) + Sqr(y2 - y1));
}

AX_END_NAMESPACE 
