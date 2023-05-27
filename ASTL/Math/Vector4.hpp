#pragma once
#include "SIMDCommon.hpp"
#include "Vector.hpp"

// todo support non SIMD Vector4
#ifdef AX_SUPPORT_SSE2

AX_ALIGNED(16) struct Vector4
{
	union
	{
		struct { float x, y, z, w; };
		float arr[4];
		__m128 vec;
	};

	static constexpr int NumElements = 4;
	static constexpr uint64 ElementSize = sizeof(float);
	using ElemType = float;

	float& operator[] (int index) { return arr[index]; }
	float operator[] (int index) const { return arr[index]; }

	constexpr Vector4() : x(0), y(0), z(0), w(0) {}
	constexpr Vector4(__m128 _vec) : vec(_vec) {}
	constexpr Vector4(float scale) : x(scale), y(scale), z(scale), w(scale) {}
	constexpr Vector4(float _x, float _y, float _z, float _w) : x(_x), y(_y), z(_z), w(_w) {}

	constexpr Vector4(const Vector3f& a, float f) : x(a.x), y(a.y), z(a.z), w(f) {}

	operator __m128() const { return vec; }

	FINLINE static __m128 VECTORCALL Normalize(const __m128 V)
	{
		return _mm_mul_ps(_mm_rsqrt_ps(_mm_dp_ps(V, V, 0xff)), V);
	}

	FINLINE static __m128 VECTORCALL Dot(const __m128 V1, const __m128 V2)
	{
		return _mm_dp_ps(V1, V2, 0xff);
	}

	float Length() { return hsum_ps_sse3(_mm_dp_ps(vec, vec, 0xff)); }

	Vector4& Normalized() { vec = _mm_mul_ps(_mm_rsqrt_ps(_mm_dp_ps(vec, vec, 0xff)), vec); return *this; }

	Vector3f xyz() const { Vector3f v;  SSEStoreVector3(&v.x, vec); return v; }

	Vector4 VECTORCALL operator + (const Vector4 b) const { return _mm_add_ps(vec, b.vec); }
	Vector4 VECTORCALL operator - (const Vector4 b) const { return _mm_sub_ps(vec, b.vec); }
	Vector4 VECTORCALL operator * (const Vector4 b) const { return _mm_mul_ps(vec, b.vec); }
	Vector4 VECTORCALL operator / (const Vector4 b) const { return _mm_div_ps(vec, b.vec); }

	Vector4& VECTORCALL operator += (const Vector4 b) { vec = _mm_add_ps(vec, b.vec); return *this; }
	Vector4& VECTORCALL operator -= (const Vector4 b) { vec = _mm_sub_ps(vec, b.vec); return *this; }
	Vector4& VECTORCALL operator *= (const Vector4 b) { vec = _mm_mul_ps(vec, b.vec); return *this; }
	Vector4& VECTORCALL operator /= (const Vector4 b) { vec = _mm_div_ps(vec, b.vec); return *this; }

	Vector4 operator  *  (const float b) const { return _mm_mul_ps(vec, _mm_set_ps1(b)); }
	Vector4 operator  /  (const float b) const { return _mm_div_ps(vec, _mm_set_ps1(b)); }
	Vector4& operator *= (const float b) { vec = _mm_mul_ps(vec, _mm_set_ps1(b)); return *this; }
	Vector4& operator /= (const float b) { vec = _mm_div_ps(vec, _mm_set_ps1(b)); return *this; }
};

AX_ALIGNED(16) struct RaySSE
{
	__m128 origin;
	__m128 direction;
};

#endif // AX_SUPPORT_SSE2

#ifdef AX_SUPPORT_AVX2

AX_ALIGNED(32) struct Vector4d
{
	union
	{
		struct { double x, y, z, w; };
		double arr[4];
		__m256d vec;
	};

	static constexpr int NumElements = 4;
	static constexpr uint64 ElementSize = sizeof(double);
	using ElemType = double;

	constexpr Vector4d() : x(0), y(0), z(0), w(0) {}
	constexpr Vector4d(double scale) : x(scale), y(scale), z(scale), w(scale) {}
	constexpr Vector4d(__m256d _vec) : vec(_vec) {}
	constexpr Vector4d(double _x, double _y, double _z, double _w) : x(_x), y(_y), z(_z), w(_w) {}

	operator __m256d() const { return vec; }

	FINLINE static __m256d VECTORCALL Normalize(__m256d V)
	{
		return _mm256_div_pd(V, _mm256_sqrt_pd(Dot(V, V))); // v / sqrt(dot(v))
	}

	FINLINE static __m256d VECTORCALL Dot(const __m256d V1, const __m256d V2)
	{
		__m256d vDot = _mm256_mul_pd(V1, V2);
		__m256d vTemp = _mm256_permute4x64_pd(vDot, _MM_SHUFFLE(2, 1, 2, 1));
		vDot = _mm256_add_pd(vDot, vTemp);
		vTemp = _mm256_permute4x64_pd(vTemp, _MM_SHUFFLE(1, 1, 1, 1));
		vDot = _mm256_add_pd(vDot, vTemp);
		return _mm256_permute4x64_pd(vDot, _MM_SHUFFLE(0, 0, 0, 0));
	}

	Vector4d& Normalized() { vec = Normalize(vec); return *this; }

	Vector4d VECTORCALL operator + (const Vector4d b) const { return _mm256_add_pd(vec, b.vec); }
	Vector4d VECTORCALL operator - (const Vector4d b) const { return _mm256_sub_pd(vec, b.vec); }
	Vector4d VECTORCALL operator * (const Vector4d b) const { return _mm256_mul_pd(vec, b.vec); }
	Vector4d VECTORCALL operator / (const Vector4d b) const { return _mm256_div_pd(vec, b.vec); }

	Vector4d& VECTORCALL operator += (const Vector4d b) { vec = _mm256_add_pd(vec, b.vec); return *this; }
	Vector4d& VECTORCALL operator -= (const Vector4d b) { vec = _mm256_sub_pd(vec, b.vec); return *this; }
	Vector4d& VECTORCALL operator *= (const Vector4d b) { vec = _mm256_mul_pd(vec, b.vec); return *this; }
	Vector4d& VECTORCALL operator /= (const Vector4d b) { vec = _mm256_div_pd(vec, b.vec); return *this; }

	Vector4d operator  *  (const double b) const { return _mm256_mul_pd(vec, _mm256_set1_pd(b)); }
	Vector4d operator  /  (const double b) const { return _mm256_div_pd(vec, _mm256_set1_pd(b)); }
	Vector4d& operator *= (const double b) { vec = _mm256_mul_pd(vec, _mm256_set1_pd(b)); return *this; }
	Vector4d& operator /= (const double b) { vec = _mm256_div_pd(vec, _mm256_set1_pd(b)); return *this; }
};

#endif // AX_SUPPORT_SSE2

typedef Vector4 float4;

