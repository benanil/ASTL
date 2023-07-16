// most of the functions are accurate and faster than stl 
// convinient for game programming

// todo IsNan function

#pragma once

#include "../Common.hpp" // includes Min, Max, Clamp, Abs and FAbs

// constants
constexpr float PI = 3.14159265358f;
constexpr float RadToDeg = 180.0f / PI;
constexpr float DegToRad = PI / 180.0f;
constexpr float OneDivPI = 1.0f / PI;
constexpr float PIDiv2 = PI / 2.0f;
constexpr float TwoPI = PI * 2.0f;
constexpr float Sqrt2 = 1.414213562f;
// for integer constants use stdint.h's INT32_MIN, INT64_MIN, INT32_MAX...
// for float constants use float.h's FLT_MAX, FLT_MIN, DBL_MAX, DBL_MIN

//  ######################################  
//  #####  [BASE FUNCTINS]  #####  
//  ######################################  

// for constant sqrt look at here: https://gist.github.com/alexshtf/eb5128b3e3e143187794
// or here for constexpr: https://gist.github.com/benanil/9d1668c0befb24263e27bd04dfa2e90f#file-mathfunctionswithoutstl-c-L230
FINLINE float Sqrt(float x) 
{
#ifdef _MSC_VER
	return _mm_cvtss_f32(_mm_sqrt_ps(_mm_set_ps1(x)));  
#else
	return __builtin_sqrt(x);
#endif
}

// original doom rsqrt implementation with comments :) maybe use _mm_rsqrt_ps instead ? and __builtin_sqrt if sse not supported
// used constant(0x5f375a86f) is from Chriss Lomont's Fast Inverse Square Root paper
FINLINE constexpr float RSqrt(float x) 
{
	float x2 = x * 0.5f; // doom rsqrt
	float y  = x;
	int i    = BitCast<int>(y);        // evil floating point bit level hacking
	i = 0x5f375a86f - ( i >> 1 );      // what the fuck? 
	y = BitCast<float>(i);              
	y = y * ( 1.5f - ( x2 * y * y ) ); // 1st iteration
	y = y * ( 1.5f - ( x2 * y * y ) ); // 2th iteration
	return y;
}

FINLINE constexpr float Fract(float a) { a = Abs(a); return a - int(a); }
									
FINLINE constexpr double Exp(double a)
{
	union { double d; long long x; } u{}, v{};
	u.x = (long long)(3248660424278399LL * a + 0x3fdf127e83d16f12LL);
	v.x = (long long)(0x3fdf127e83d16f12LL - 3248660424278399LL * a);
	return u.d / v.d;
}

// https://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/
// https://github.com/ekmett/approximate/blob/master/cbits/fast.c#L81
// should be much more precise with large b
FINLINE constexpr float Pow(float a, float b) 
{
	// calculate approximation with fraction of the exponent
	int e = (int) b;
	union {
		double d;
		int x[2];
	} u = { (double)a };
	u.x[1] = (int)((b - e) * (u.x[1] - 1072632447) + 1072632447);
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
FINLINE constexpr float Log(float x)
{
	return (BitCast<int>(x) - 1064866805) * 8.262958405176314e-8f;
}

// if you want log10 for integer you can look at Algorithms.hpp
FINLINE constexpr float Log10(float x)
{ 
	return Log(x) / 2.30258509299f; // ln(x) / ln(10)
} 

// you might look at this link as well: https://tech.ebayinc.com/engineering/fast-approximate-logarithms-part-i-the-basics/
FINLINE constexpr float Log2(float x)
{
	return Log(x) / 0.6931471805599453094f; // ln(x) / ln(2) 
}

FINLINE constexpr bool IsZero(float x) noexcept { return Abs(x) <= 0.0001f; }
FINLINE constexpr bool AlmostEqual(float x, float  y) noexcept { return Abs(x-y) <= 0.001f; }
FINLINE constexpr float Sign(float x, float y) { return y < 0.0f ? -x : x; } //!< maybe use templates?

FINLINE constexpr float CopySign(float x, float y) 
{
	int ix = BitCast<int>(x) & 0x7fffffff;
	int iy = BitCast<int>(y) & 0x80000000;
	return BitCast<float>(ix | iy);
}

FINLINE constexpr bool IsNan(float f)
{
	uint32 intValue = BitCast<uint32>(f);
	uint32 exponent = (intValue >> 23) & 0xFF;
	uint32 fraction = intValue & 0x7FFFFF;
	return (exponent == 0xFF) & (fraction != 0);
}

FINLINE constexpr float FMod(float x, float y) {
	float quotient = x / y;
	float whole = (float)((int)quotient);  // truncate quotient to integer
	float remainder = x - whole * y;
	remainder += (remainder < 0.0f) * y;
	return remainder;
}

FINLINE constexpr float Floor(float x) {
	float whole = (float)(int)x;  // truncate quotient to integer
	return x - (x-whole);
}

template<typename T> FINLINE constexpr 
bool IsPowerOfTwo(T x) noexcept { return (x != 0) && ((x & (x - 1)) == 0); }

//  ######################################  
//  #####  [TRIGONOMETRIC FUNCTINS]  #####  
//  ######################################  

// https://mazzo.li/posts/vectorized-atan2.html
FINLINE constexpr float ATan(float x) {
	const float x_sq = x * x;
	constexpr float a1 =  0.99997726f, a3 = -0.33262347f, a5  = 0.19354346f,
	                a7 = -0.11643287f, a9 =  0.05265332f, a11 = -0.01172120f;
	return x * (a1 + x_sq * (a3 + x_sq * (a5 + x_sq * (a7 + x_sq * (a9 + x_sq * a11)))));
}

FINLINE constexpr float ATan2(float y, float x) {
	// from here: https://yal.cc/fast-atan2/  
	float ay = Abs(y), ax = Abs(x);
	int invert = ay > ax;
	float z = invert ? ax/ay : ay/ax;// [0,1]
	float th = ATan(z);              // [0,π/4]
	if(invert) th = PIDiv2 - th;     // [0,π/2]
	if(x < 0)  th = PI - th;         // [0,π]
	return CopySign(th, y);          // [-π,π] 
}

FINLINE float ASin(float z) 
{
	return ATan2(z, Sqrt(1.0f-(z * z)));
}

// https://en.wikipedia.org/wiki/Sine_and_cosine
// warning: accepts input between -TwoPi and TwoPi  if (Abs(x) > TwoPi) use x = FMod(x + PI, TwoPI) - PI;
FINLINE constexpr float Sin(float x) 
{
	float xx = x * x * x;                // x^3
	float t = x - (xx * 0.16666666666f); // x3/!3  6 = !3 = 1.6666 = rcp(3)
	t += (xx *= x * x) * 0.00833333333f; // x5/!5  120 = !5 = 0.0083 = rcp(5)
	t -= (xx *= x * x) * 0.00019841269f; // x7/!7  5040 = !7
	t += (xx * x * x)  * 2.75573e-06f;   // 362880 = !9
	return t;
}

// warning: accepts input between -TwoPi and TwoPi  if (Abs(x) > TwoPi) use x = FMod(x + PI, TwoPI) - PI;
FINLINE constexpr float Cos(float x) 
{
	float xx = x * x;                    // x^2
	float t = 1.0f - (xx * 0.5f);        // 1-(x2/!2) 0.5f = rcp(!2)
	t += (xx *= x * x) * 0.04166666666f; // t + (x4/!4) 0.04166 = rcp(24) = 24.0f = !4
	t -= (xx *= x * x) * 0.00138888888f; // t - (x6/!6) 720.0f = !6
	t += (xx * x * x)  * 2.48016e-05f;   // t + (x8/!8) 40320 = !8
	return t;
}

FINLINE constexpr void SinCos(float x, float* s, float* c)
{
	*s = Sin(x); *c = Cos(x); 
}

FINLINE constexpr float Tan(float radians) 
{
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
FINLINE float ATan2PI(float y, float x) { return ATan2(y, x) / PI; }
FINLINE float ASinPI(float z) { return ASin(z) / PI; }
FINLINE float ACos(float x)   { return PIDiv2 - ASin(x); }
FINLINE float ACosPI(float x) { return ACos(x) / PI; }
FINLINE float CosPI(float x)  { return Cos(x) / PI; }
FINLINE float SinPI(float x)  { return Sin(x) / PI; }

//  ######################################  
//  #####      [HALF FUNCTINS]       #####  
//  ######################################  

typedef ushort half;

// taken from stack overflow
FINLINE float ConvertHalfToFloat(half x)
{
	const uint e = (x & 0x7C00) >> 10; // exponent
	const uint m = (x & 0x03FF) << 13; // mantissa
	const uint v = BitCast<uint>((float)m) >> 23; // evil log2 bit hack to count leading zeros in denormalized format
	uint a = (x & 0x8000) << 16 | (e != 0) * ((e + 112) << 23 | m);
	a |= ((e == 0) & (m != 0)) * ((v - 37) << 23 | ((m << (150 - v)) & 0x007FE000));
	return BitCast<float>(a); // sign : normalized : denormalized
}

FINLINE half ConvertFloatToHalf(float Value)
{
	const uint b = BitCast<uint>(Value) + 0x00001000; // round-to-nearest-even: add last bit after truncated mantissa
	const uint e = (b & 0x7F800000) >> 23; // exponent
	const uint m = b & 0x007FFFFF; // mantissa; in line below: 0x007FF000 = 0x00800000-0x00001000 = decimal indicator flag - initial rounding
	uint a = (b & 0x80000000) >> 16 | (e > 112) * ((((e - 112) << 10) & 0x7C00) | m >> 13);
	return a | ((e < 113) & (e > 101))*((((0x007FF000 + m) >> (125-e)) + 1) >> 1) | (e > 143) * 0x7FFF; // sign : normalized : denormalized : saturate
}

// taken from XNA math
FINLINE float ConvertHalfToFloatSafe(half Value)
{
	uint Mantissa, Exponent, Result;
	Mantissa = (uint)(Value & 0x03FF);

	if ((Value & 0x7C00) != 0)  // The value is normalized
	{
		Exponent = (uint)((Value >> 10) & 0x1F);
	}
	else if (Mantissa != 0)     // The value is denormalized
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

	Result = ((Value & 0x8000) << 16) | ((Exponent + 112) << 23) | (Mantissa << 13);
	return *(float*)&Result;
}

FINLINE half ConvertFloatToHalfSafe(float Value)
{
	uint Result;
	uint IValue = ((uint*)(&Value))[0];
	uint Sign = (IValue & 0x80000000U) >> 16U;
	IValue = IValue & 0x7FFFFFFFU;      // Hack off the sign

	if (IValue > 0x47FFEFFFU) {
		// The number is too large to be represented as a half.  Saturate to infinity.
		Result = 0x7FFFU;
	}
	else if (IValue < 0x38800000U) {
		// The number is too small to be represented as a normalized half.
		// Convert it to a denormalized value.
		uint Shift = 113U - (IValue >> 23U);
		IValue = (0x800000U | (IValue & 0x7FFFFFU)) >> Shift;
	}
	else {
		// Rebias the exponent to represent the value as a normalized half.
		IValue += 0xC8000000U;
	}

	Result = ((IValue + 0x0FFFU + ((IValue >> 13U) & 1U)) >> 13U) & 0x7FFFU; 
	return (half)(Result | Sign);
}

//  ######################################  
//  #####      [COLOR FUNCTINS]      #####  
//  ######################################  

FINLINE uint PackColorRGBU32(float r, float g, float b) {
	return (uint)(r * 255.0f) | ((uint)(g * 255.0f) << 8) | ((uint)(b * 255.0f) << 16);
}

FINLINE uint PackColorRGBU32(float* c) {
	return (uint)(*c * 255.0f) | ((uint)(c[1] * 255.0f) << 8) | ((uint)(c[2] * 255.0f) << 16);
}

FINLINE uint PackColorRGBAU32(float* c) {
	return (uint)(*c * 255.0f) | ((uint)(c[1] * 255.0f) << 8) | ((uint)(c[2] * 255.0f) << 16) | ((uint)(c[3] * 255.0f) << 24);
}

FINLINE void UnpackColorRGBf(unsigned color, float* colorf)
{
	constexpr float toFloat = 1.0f / 255.0f;
	colorf[0] = float(color >> 0  & 0xFF) * toFloat;
	colorf[1] = float(color >> 8  & 0xFF) * toFloat;
	colorf[2] = float(color >> 16 & 0xFF) * toFloat;
}

FINLINE void UnpackColorRGBAf(unsigned color, float* colorf) {
	constexpr float toFloat = 1.0f / 255.0f;
	colorf[0] = float(color >> 0  & 0xFF) * toFloat;
	colorf[1] = float(color >> 8  & 0xFF) * toFloat;
	colorf[2] = float(color >> 16 & 0xFF) * toFloat;
	colorf[3] = float(color >> 24) * toFloat;
}
