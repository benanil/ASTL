#pragma once
#include "../Common.hpp"
#include <math.h>

// constants
constexpr float PI = 3.14159265358f;
constexpr float PI_2 = 1.5707963267f;  // pi/2
constexpr float RadToDeg = 180.0f / PI;
constexpr float DegToRad = PI / 180.0f;
constexpr float OneDivPI = 1.0f / PI;
constexpr float PIDiv2 = PI / 2.0f;
constexpr float TwoPI = PI * 2.0f;

//  ######################################  
//  #####  [BASE FUNCTINS]  #####  
//  ######################################  

FINLINE _NODISCARD float Pow(float x, float y) { return powf(x, y);   }
FINLINE _NODISCARD float Sqrt(float x)         { return sqrtf(x);  }
FINLINE _NODISCARD float Log(float x)          { return logf(x);   }
FINLINE _NODISCARD float Log10(float x)        { return log10f(x); } // if you want log10 for integer you can look at Algorithms.hpp
FINLINE _NODISCARD float Log2(float x)         { return log2f(x);  }

FINLINE constexpr _NODISCARD float FAbs(float x) { return x < 0.0f ? -x : x; }
FINLINE constexpr _NODISCARD double FAbs(double x) { return x < 0.0 ? -x : x; }
FINLINE constexpr _NODISCARD bool IsZero(float x) noexcept { return FAbs(x) <= 0.0001f; }
FINLINE constexpr _NODISCARD bool AlmostEqual(float x, float  y) noexcept { return FAbs(x-y) <= 0.001f; }

FINLINE constexpr _NODISCARD float FMod(float x, float y) {
	float quotient = x / y;
	float whole = (float)((int)quotient);  // truncate quotient to integer
	float remainder = x - whole * y;
	remainder += (remainder < 0.0f) * y;
	return remainder;
}

FINLINE constexpr _NODISCARD float Floor(float x) {
	float whole = (float)(int)x;  // truncate quotient to integer
	return x - (x-whole);
}

template<typename T> FINLINE constexpr 
bool IsPowerOfTwo(T x) noexcept { return !(x&1) & (x != 0); }

//  ######################################  
//  #####  [TRIGONOMETRIC FUNCTINS]  #####  
//  ######################################  

// https://mazzo.li/posts/vectorized-atan2.html
FINLINE constexpr _NODISCARD float ATan(float x) {
	const float x_sq = x * x;
	constexpr float a1 =  0.99997726f, a3 = -0.33262347f, a5  = 0.19354346f,
	                a7 = -0.11643287f, a9 =  0.05265332f, a11 = -0.01172120f;
	return x * (a1 + x_sq * (a3 + x_sq * (a5 + x_sq * (a7 + x_sq * (a9 + x_sq * a11)))));
}

// from here: https://yal.cc/fast-atan2/  
FINLINE _NODISCARD float ATan2(float y, float x) {
	float ay = FAbs(y), ax = FAbs(x);
	int invert = ay > ax;
	float z = invert ? ax/ay : ay/ax;// [0,1]
	float th = ATan(z);              // [0,π/4]
	if(invert) th = PI_2 - th;       // [0,π/2]
	if(x < 0)  th = PI - th;         // [0,π]
	// with removing this function we can make this function constexpr
	// and currently not compatible with other than MSVC
	return __copysignf(th, y);       // [-π,π] 
}

FINLINE _NODISCARD float ASin(float z) 
{
	return ATan2(z, Sqrt(1.0f-(z * z)));
	// this way is probably much faster, 
	// vectorizable and can be constexpr 
	// but it doesn't work correctly when value is closer to -1
	float z3 = z * z * z;
	float t = z + (0.500000f * z3 * 0.333333333f); // z3 / 3.0f 0.333 = rcp(3)
	t += (0.500000000f) * ((z3 *= z * z) * 0.200000000f);  // z5 / 5.0f 0.2 = rcp(5)
	t += (0.375000000f) * ((z3 *= z * z) * 0.142857143f); // z7 / 9.0f
	t += (0.234375000f) * ((z3 *= z * z) * 0.111111111f);
	t += (0.102539062f) * ((z3 *= z * z) * 0.090909091f);
	t += (0.028839111f) * ((z3 *= z * z) * 0.076923077f);
	t += (0.004956722f) * ((z3 * z * z)  * 0.066666667f);
	return t; //0.000029            //* 0.058824
}

FINLINE constexpr _NODISCARD float Sin(float x) 
{
	x = FMod(x + PI, TwoPI) - PI;
	float xx = x * x * x;                // x^3
	float t = x - (xx * 0.16666666666f); // x3/!3  6 = !3 = 1.6666 = rcp(3)
	t += (xx *= x * x) * 0.00833333333f; // x5/!5  120 = !5 = 0.0083 = rcp(5)
	t -= (xx *= x * x) * 0.00019841269f; // x7/!7  5040 = !7
	t += (xx * x * x) / 362880.0f;       // 362880 = !9
	return t;
}

FINLINE constexpr _NODISCARD float Cos(float x) 
{
	x = FMod(x + PI, TwoPI) - PI;
	float xx = x * x;                    // x^2
	float t = 1.0f - (xx * 0.5f);        // 1-(x2/!2) 0.5f = rcp(!2)
	t += (xx *= x * x) * 0.04166666666f; // t + (x4/!4) 0.04166 = rcp(24) = 24.0f = !4
	t -= (xx *= x * x) * 0.00138888888f; // t - (x6/!6) 720.0f = !6
	t += (xx * x * x) / 40320.0f;        // t + (x8/!8) 40320 = !8
	return t;
}

FINLINE constexpr _NODISCARD void SinCos(float x, float* s, float* c)
{
	*s = Sin(x); *c = Cos(x); 
}

// from chat gpt
FINLINE constexpr _NODISCARD float Tan(float x, int depth = 10) 
{
	if (depth == 0.0f) return x;
	float h = Tan(x * 0.5f, depth - 1);
	return (2.0f * h) / (1.0f - (h * h));
}

// inspired from Casey Muratori's performance aware programming
// this functions makes the code more readable. OpenCL and Cuda has the same functions as well
FINLINE _NODISCARD float ATan2PI(float y, float x) { return ATan2(y, x) / PI; }
FINLINE _NODISCARD float ASinPI(float z) { return ASin(z) / PI; }
FINLINE _NODISCARD float ACos(float x)   { return PIDiv2 - ASin(x); }
FINLINE _NODISCARD float ACosPI(float x) { return ACos(x) / PI; }
FINLINE _NODISCARD float CosPI(float x)  { return Cos(x) / PI; }
FINLINE _NODISCARD float SinPI(float x)  { return Sin(x) / PI; }

FINLINE _NODISCARD constexpr float Tan10(float x) { return Tan(x); }

//  ######################################  
//  #####      [MISC FUNCTINS]       #####  
//  ######################################  

template<typename RealT>
FINLINE _NODISCARD RealT Lerp(const RealT from, const RealT to,
                              const RealT t) noexcept {
	return from + (to - from) * t;
}

// original doom rsqrt implementation with comments :)
FINLINE _NODISCARD constexpr float RSqrt(float number) 
{
	float x2 = number * 0.5F;
	float y  = number;
	long i   = *(long*) &y;              // evil floating point bit level hacking
	i = 0x5f3759df - ( i >> 1 );         // what the fuck? 
	y = * ( float * ) &i;
	y = y * ( 1.5f - ( x2 * y * y ) );   // 1st iteration
	return y;
}

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
