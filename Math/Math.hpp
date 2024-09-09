
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
constexpr float PI        = 3.14159265358f;
constexpr float HalfPI    = PI / 2.0f;
constexpr float QuarterPI = PI / 4.0f;
constexpr float RadToDeg  = 180.0f / PI;
constexpr float DegToRad  = PI / 180.0f;
constexpr float OneDivPI  = 1.0f / PI;
constexpr float TwoPI     = PI * 2.0f;
constexpr float Sqrt2     = 1.414213562f;
constexpr float Epsilon   = 0.0001f;
// for integer constants use stdint.h's INT32_MIN, INT64_MIN, INT32_MAX...
// for float constants use float.h's FLT_MAX, FLT_MIN, DBL_MAX, DBL_MIN

/*//////////////////////////////////////////////////////////////////////////*/
/*                              Essential                                   */
/*//////////////////////////////////////////////////////////////////////////*/

pureconst float Sqr(float x) {
    return x * x;
}

pureconst float Lerp(float x, float y, float t) {
    return x + (y - x) * t;
}

pureconst double Lerp(double x, double y, double t) {
    return x + (y - x) * t;
}

pureconst double SqrtConstexpr(double a) {
    // from: Jack W. Cerenshaw's math toolkit for real time development book: page 63 Listing 4.3
    // I've removed some of the branches. slightly slower than sqrtf
    // double version is also here: https://gist.github.com/benanil/9d1668c0befb24263e27bd04dfa2e90f#file-mathfunctionswithoutstl-c-L230
    const double A = 0.417319242, B = 0.5901788532;
    union { double fp; struct { unsigned lo, hi; }; } x {};
    if (a <= 0.001) return 0.0;
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

pureconst float SqrtConstexpr(float x) {
    float y = BitCast<float>(532545536u + (BitCast<uint32_t>(x) >> 1u));
    y = 0.5f * (y + x / y);
    y = 0.5f * (y + x / y);
    return y;
}

// cbrt, std::cbrt
pureconst float CubeRootConstexpr(float x)
{
    float y = BitCast<float>(709973695u + BitCast<uint32_t>(x)/3u);
    y = y * (2.0f / 3.0f) + (1.0f / 3.0f) * x / (y * y);
    y = y * (2.0f / 3.0f) + (1.0f / 3.0f) * x / (y * y);
    return y;
}

purefn float Sqrt(float a) {
#ifdef AX_SUPPORT_SSE
    return _mm_cvtss_f32(_mm_sqrt_ps(_mm_set_ps1(a)));
#elif defined(__clang__)
    return __builtin_sqrt(a);
#else
    return SqrtConstexpr(a);
#endif
}

// https://en.wikipedia.org/wiki/Fast_inverse_square_root
// https://rrrola.wz.cz/inv_sqrt.html   <- fast and 2.5x accurate
purefn float RSqrt(float x) {
#ifdef AX_SUPPORT_SSE
    // when I compile with godbolt -O1 expands to one instruction vrsqrtss
    // which is approximately equal latency as mulps.
    // RSqrt(float):                             # @RSqrt3(float)
    //       vrsqrtss        xmm0, xmm0, xmm0
    //       ret
    return _mm_cvtss_f32(_mm_rsqrt_ps(_mm_set_ps1(x)));
#elif defined(AX_SUPPORT_NEON)
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
pureconst float Exp(float f) {
    const int IEEE_FLT_MANTISSA_BITS  =	23;
    const int IEEE_FLT_EXPONENT_BITS  =	8;
    const int IEEE_FLT_EXPONENT_BIAS  =	127;
    const int IEEE_FLT_SIGN_BIT       =	31;
    
    float x = f * 1.44269504088896340f; // multiply with ( 1 / log( 2 ) )
    
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
pureconst float Pow(float a, float b) {
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
pureconst float Log(float x) {
    return (BitCast<int>(x) - 1064866805) * 8.262958405176314e-8f;
}

// if you want log10 for integer you can look at Algorithms.hpp
pureconst float Log10(float x) { 
    return Log(x) / 2.30258509299f; // ln(x) / ln(10)
} 

// you might look at this link as well: https://tech.ebayinc.com/engineering/fast-approximate-logarithms-part-i-the-basics/
pureconst float Log2(float x) {
    return Log(x) / 0.6931471805599453094f; // ln(x) / ln(2) 
}

// https://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
pureconst unsigned int Log2(unsigned int v) {
    constexpr int MultiplyDeBruijnBitPosition[32] = 
    {
        0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
        8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31
    };
    
    // first round down to one less than a power of 2
    v |= v >> 1; v |= v >> 2; v |= v >> 4; v |= v >> 8; v |= v >> 16;
    return MultiplyDeBruijnBitPosition[(uint32_t)(v * 0x07C4ACDDU) >> 27];
}

pureconst unsigned int Log10(unsigned int v) {
                                                    // this would work too
    unsigned int const PowersOf10[] = {             // if (n <= 9) return 1;
        1, 10, 100, 1000, 10000, 100000,            // if (n <= 99) return 2;
        1000000, 10000000, 100000000, 1000000000    // if (n <= 999) return 3;
    };                                              // ...
    unsigned int t = (Log2(v) + 1) * 1233 >> 12;    // if (n <= 2147483647) return 10; // 2147483647 = int max
    return t - (v < PowersOf10[t]);                                                      
}                                                                                        

pureconst bool IsZero(float x) {
    return Abs(x) <= 0.0001f; 
}

pureconst bool AlmostEqual(float x, float  y) {
    return Abs(x-y) <= 0.001f;
}

pureconst float Sign(float x) {
    int res = BitCast<int>(1.0f);
    res |= BitCast<int>(x) & 0x80000000;
    return BitCast<float>(res);
} 

pureconst int Sign(int x) {
    return x < 0 ? -1 : 1; // equal to above float version
} 

pureconst float CopySign(float x, float y) {
    int ix = BitCast<int>(x) & 0x7fffffff;
    int iy = BitCast<int>(y) & 0x80000000;
    return BitCast<float>(ix | iy);
}

pureconst bool IsNan(float f) {
    uint32 intValue = BitCast<uint32>(f);
    uint32 exponent = (intValue >> 23) & 0xFF;
    uint32 fraction = intValue & 0x7FFFFF;
    return (exponent == 0xFF) && (fraction != 0);
}

template<typename T>
pureconst T FMod(T x, T y) {
    T quotient = x / y;
    T whole = (T)((int)quotient);  // truncate quotient to integer
    T remainder = x - whole * y;
    remainder += (remainder < (T)0.0) * y;
    return remainder;
}

template<typename T>
pureconst T Floor(T x) {
    T whole = (T)(int)x;  // truncate quotient to integer
    return x - (x-whole);
}

template<typename T>
pureconst T Ceil(T x) {
    T whole = (T)(int)x;  // truncate quotient to integer
    return whole + float(x > whole);
}

pureconst float Fract(float a) {
    return a - int(a); 
}

/*//////////////////////////////////////////////////////////////////////////*/
/*                      Trigonometric Functions                             */
/*//////////////////////////////////////////////////////////////////////////*/

// https://mazzo.li/posts/vectorized-atan2.html
pureconst float ATan(float x) {
    const float x_sq = x * x;
    const float a1 =  0.99997726f, a3 = -0.33262347f, a5  = 0.19354346f,
                a7 = -0.11643287f, a9 =  0.05265332f, a11 = -0.01172120f;
    return x * (a1 + x_sq * (a3 + x_sq * (a5 + x_sq * (a7 + x_sq * (a9 + x_sq * a11)))));
}

// Warning! if y and x is zero this will return HalfPI instead of 0.0f unlike cstdlib
pureconst float ATan2(float y, float x) {
    // https://gist.github.com/volkansalma/2972237
    float abs_y = Abs(y) + 1e-10f;      // kludge to prevent 0/0 condition
    float r = (x - CopySign(abs_y, x)) / (abs_y + Abs(x));
    float angle = HalfPI - CopySign(PI / 4.f, x);
    angle += (0.1963f * r * r - 0.9817f) * r;
    return CopySign(angle, y);
}

purefn float ASin(float z) {
    return ATan2(z, Sqrt(1.0f-(z * z)));
}

// Valid input range -1..1 output is -pi..pi
purefn float ACos(float x)   
{
    // Lagarde 2014, "Inverse trigonometric functions GPU optimization for AMD GCN architecture"
    // This is the approximation of degree 1, with a max absolute error of 9.0x10^-3
    float y = Abs(x);
    float p = -0.1565827f * y + 1.570796f;
    p *= Sqrt(1.0f - y);
    return x >= 0.0f ? p : PI - p;
}

purefn float ACosPositive(float x)
{
    float p = -0.1565827f * x + 1.570796f;
    return p * Sqrt(1.0f - x);
}

// https://en.wikipedia.org/wiki/Sine_and_cosine
// warning: accepts input between -TwoPi and TwoPi  if (Abs(x) > TwoPi) use x = FMod(x + PI, TwoPI) - PI;

pureconst float RepeatPI(float x) {
    return FMod(x + PI, TwoPI) - PI;
}

// Accepts input between -TwoPi and TwoPi, use SinR if value is bigger than this range
// https://x.com/nthnblair/status/1790836310531559701
// https://www.desmos.com/calculator/owkgky7wh4
// https://github.com/LancePutnam/Gamma/blob/master/Gamma/scl.h  :function sinfast
pureconst float Sin(float x) 
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
pureconst float Cos(float a)
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
pureconst float Sin0pi(float x) {
    x *= 0.63655f; // constant founded using desmos
    x *= 2.0f - x;
    return x * (0.225f * x + 0.775f); 
}

// calculates cos(x) between [0,pi]
pureconst float Cos0pi(float a) {
    a *= 0.159f;
    return 1.0f - 32.0f * a * a * (0.75f - a);
}

// R suffix allows us to use with greater range than -TwoPI, TwoPI
pureconst float SinR(float x) {
    return Sin(FMod(x + PI, TwoPI) - PI);
}

// R suffix allows us to use with greater range than -TwoPI, TwoPI
pureconst float CosR(float x) {
    return Cos(FMod(x + PI, TwoPI) - PI);
}

forceinline void SinCos(float x, float* sp, float* cp) 
{
    *sp = Sin(x);
    *cp = Cos(x);
}

// https://github.com/id-Software/DOOM-3/blob/master/neo/idlib/math/Math.h
// tanf equivalent
pureconst float Tan(float a) {
    float s = 0.0f;
    bool reciprocal = false;

    if (( a < 0.0f ) || (a >= PI)) {
        a -= Floor(a / PI) * PI;
    }
    
    if ( a < HalfPI ) {
        bool greater = a > QuarterPI;
        reciprocal = greater;
        if (greater) a = HalfPI - a;
    } else {
        bool greater = a > HalfPI + QuarterPI;
        reciprocal = !greater;
        a = greater ? a - PI : HalfPI - a;
    }

    s = a * a;
    s = a * ((((((9.5168091e-03f * s + 2.900525e-03f ) * s + 2.45650893e-02f) * s + 5.33740603e-02f) * s + 1.333923995e-01f) * s + 3.333314036e-01f) * s + 1.0f);
    return reciprocal ? 1.0f / s : s;
}

// inspired from Casey Muratori's performance aware programming
// this functions makes the code more readable. OpenCL and Cuda has the same functions as well
pureconst float ATan2PI(float y, float x) { return ATan2(y, x) / PI; }
purefn    float ASinPI(float z) { return ASin(z) / PI; }
purefn    float ACosPI(float x) { return ACos(x) / PI; }
pureconst float CosPI(float x)  { return Cos(x) / PI; }
pureconst float SinPI(float x)  { return Sin(x) / PI; }

// packs -1,1 range float to short
pureconst short PackSnorm16(float x) {
    return (short)Clamp(x * (float)INT16_MAX, (float)INT16_MIN, (float)INT16_MAX);
}

pureconst float UnpackSnorm16(short x) {
    return (float)x / (float)INT16_MAX;
}

// packs 0,1 range float to short
pureconst short PackUnorm16(float x) {
    return (short)Clamp(x * (float)INT16_MAX, (float)INT16_MIN, (float)INT16_MAX);
}

pureconst float UnpackUnorm16(short x) {
    return (float)x / (float)UINT16_MAX;
}

/*//////////////////////////////////////////////////////////////////////////*/
/*                                 Easing                                   */
/*//////////////////////////////////////////////////////////////////////////*/

// to see visually: https://easings.net/ 
pureconst float EaseIn(float x) {
    return x * x;
}

pureconst float EaseOut(float x) { 
    float r = 1.0f - x;
    return 1.0f - (r * r); 
}

pureconst float EaseInOut(float x) {
    return x < 0.5f ? 2.0f * x * x : 1.0f - Sqr(-2.0f * x + 2.0f) / 2.0f;
}

// integral symbol shaped interpolation, similar to EaseInOut
pureconst float SmoothStep(float x) {
    return x * x * (3.0f - x * 2.0f);
}

pureconst float EaseInSine(float x) {
  return 1.0f - Cos((x * PI) * 0.5f);
}

pureconst float EaseOutSine(float x) {
  return Sin((x * PI) * 0.5f);
}

// Gradually changes a value towards a desired goal over time.
purefn float SmoothDamp(float current, float target, float& currentVelocity, float smoothTime, float maxSpeed, float deltaTime)
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

purefn float Remap(float in, float inMin, float inMax, float outMin, float outMax)
{
    return outMin + (in - inMin) * (outMax - outMin) / (inMax - inMin);
}

purefn float Repeat(float t, float length)
{
    return Clamp(t - Floor(t / length) * length, 0.0f, length);
}

purefn float Step(float edge, float x)
{
    return float(x > edge);
}

// https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
// position(x0, y0), linePos1(x1, y1), linePos2(x2, y2) 
purefn float LineDistance(float x0, float y0, float x1, float y1, float x2, float y2)
{
    float a = ((x2 - x1) * (y0 - y1)) - ((x0 - x1) * (y2 - y1));
    return Abs(a) * RSqrt(Sqr(x2 - x1) + Sqr(y2 - y1));
}

AX_END_NAMESPACE 
