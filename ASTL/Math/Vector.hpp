#pragma once
#include "../Random.hpp" // includes math.h

template<typename T>
struct Vector2
{
	union
	{
		struct { T x, y; };
		T arr[2];
	};
	
	static constexpr int NumElements = 2;
	static constexpr uint64 ElementSize = sizeof(T);
	using ElemType = T;

	Vector2()                   : x(0), y(0) { }
	constexpr Vector2(T s)      : x(s), y(s) { }
	constexpr Vector2(T a, T b) : x(a), y(b) { }
	
	float Length()		  const { return Sqrt(LengthSquared()); }
	float LengthSquared() const { return x * x + y * y; }

	T& operator[] (int index) { return arr[index]; }
	T  operator[] (int index) const { return arr[index]; }

	static FINLINE float Distance(Vector2 a, Vector2 b) {
		float diffx = (float)(a.x - b.x);
		float diffy = (float)(a.y - b.y);
		return Sqrt(diffx * diffx + diffy * diffy);
	}

	static FINLINE float DistanceSq(Vector2 a, Vector2 b) {
		float diffx = (float)(a.x - b.x);
		float diffy = (float)(a.y - b.y);
		return diffx * diffx + diffy * diffy;
	}

	static FINLINE Vector2 Rotate(Vector2 vec, float angle)
	{
		float s = Sin(angle), c = Cos(angle);
		return Vector2(vec.x * c - vec.y * s, vec.x * s + vec.y * c);
	}

	void Normalized() const { *this /= Length(); }
	
	Vector2 Normalize(Vector2 other) { return other.Normalize(); }
	
	Vector2 operator - () { return Vector2(-x, -y); }
	Vector2 operator + (Vector2 other) const { return Vector2(x + other.x, y + other.y); }
	Vector2 operator * (Vector2 other) const { return Vector2(x * other.x, y * other.y); }
	Vector2 operator / (Vector2 other) const { return Vector2(x / other.x, y / other.y); }
	Vector2 operator - (Vector2 other) const { return Vector2(x - other.x, y - other.y); }

	Vector2 operator + (T other) const { return Vector2(x + other, y + other); }
	Vector2 operator * (T other) const { return Vector2(x * other, y * other); }
	Vector2 operator / (T other) const { return Vector2(x / other, y / other); }
	Vector2 operator - (T other) const { return Vector2(x - other, y - other); }
	
	Vector2 operator += (Vector2 o) { x += o.x; y += o.y; return *this; }
	Vector2 operator *= (Vector2 o) { x *= o.x; y *= o.y; return *this; }
	Vector2 operator /= (Vector2 o) { x /= o.x; y /= o.y; return *this; }
	Vector2 operator -= (Vector2 o) { x -= o.x; y -= o.y; return *this; }

	Vector2 operator += (T o) { x += o; y += o; return *this; }
	Vector2 operator *= (T o) { x *= o; y *= o; return *this; }
	Vector2 operator /= (T o) { x /= o; y /= o; return *this; }
	Vector2 operator -= (T o) { x -= o; y -= o; return *this; }

	bool operator == (Vector2 b) const { return x == b.x && y == b.y; }
	bool operator != (Vector2 b) const { return x != b.x || y != b.y; }
	bool operator <  (Vector2 b) const { return x < b.x && y < b.y; }

	static constexpr Vector2 Zero()     { return Vector2( 0,  0); } 
	static constexpr Vector2 One()      { return Vector2( 1,  1); } 
	static constexpr Vector2 Up()       { return Vector2( 0,  1); } 
	static constexpr Vector2 Left()     { return Vector2(-1,  0); } 
	static constexpr Vector2 Down()     { return Vector2( 0, -1); } 
	static constexpr Vector2 Right()    { return Vector2( 1,  0); } 
};

template<typename T>
struct Vector3
{
	union
	{
		struct { T x, y, z; };
		T arr[3];
	};

	static constexpr int NumElements = 3;
	using ElemType = T;

	Vector3() : x(0), y(0), z(0) { }
	constexpr Vector3(T scale) : x(scale), y(scale), z(scale) { }
	constexpr Vector3(T a, T b, T c) : x(a), y(b), z(c) { }

	T& operator[] (int index) { return arr[index]; }
	T  operator[] (int index) const { return arr[index]; }

	float Length() const { return sqrtf(LengthSquared()); }
	constexpr float LengthSquared() const { return x * x + y * y + z * z; }

	Vector3& Normalized() { *this /= Length(); return *this; }
	void NormalizeSelf() { *this /= Length(); }

	static float Distance(const Vector3& a, const Vector3& b)
	{
		Vector3 diff = a - b; diff *= diff;
		return Sqrt(diff.x + diff.y + diff.z);
	}

	// distance squared for comparing distances, faster than distance 
	static float DistanceSq(const Vector3& a, const Vector3& b)
	{
		Vector3 diff = a - b; diff *= diff;
		return diff.x + diff.y + diff.z;
	}

	static float Length(const Vector3& vec) { return vec.Length(); }

	static float Dot(const Vector3& a, const Vector3& b)
	{
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}

	static Vector3 Lerp(const Vector3& a, const Vector3& b, float t)
	{
		return Vector3(
			a.x + (b.x - a.x) * t,
			a.y + (b.y - a.y) * t,
			a.z + (b.z - a.z) * t
		);
	}

	static Vector3 Cross(const Vector3& a, const Vector3& b)
	{
		return Vector3(a.y * b.z - b.y * a.z,
			           a.z * b.x - b.z * a.x,
			           a.x * b.y - b.x * a.y);
	}

	static Vector3 Reflect(const Vector3& in, const Vector3& normal)
	{
		return in - normal * Vector3::Dot(normal, in) * 2.0f;
	}
	// for more accuracy you can use sqrt instead of rsqrt: a / sqrt(dot(a,a)) 
	static Vector3 Normalize(const Vector3& a) {
		return a / Sqrt(Vector3::Dot(a, a));
	}

	Vector3 operator - () { return Vector3(-x, -y, -z); }
	Vector3 operator + (const Vector3& other) const { return Vector3(x + other.x, y + other.y, z + other.z); }
	Vector3 operator * (const Vector3& other) const { return Vector3(x * other.x, y * other.y, z * other.z); }
	Vector3 operator / (const Vector3& other) const { return Vector3(x / other.x, y / other.y, z / other.z); }
	Vector3 operator - (const Vector3& other) const { return Vector3(x - other.x, y - other.y, z - other.z); }

	Vector3 operator + (T other) const { return Vector3(x + other, y + other, z + other); }
	Vector3 operator * (T other) const { return Vector3(x * other, y * other, z * other); }
	Vector3 operator / (T other) const { return Vector3(x / other, y / other, z / other); }
	Vector3 operator - (T other) const { return Vector3(x - other, y - other, z - other); }

	Vector3 operator += (const Vector3& o) { x += o.x; y += o.y; z += o.z; return *this; }
	Vector3 operator *= (const Vector3& o) { x *= o.x; y *= o.y; z *= o.z; return *this; }
	Vector3 operator /= (const Vector3& o) { x /= o.x; y /= o.y; z /= o.z; return *this; }
	Vector3 operator -= (const Vector3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }

	Vector3 operator += (T o) { x += o; y += o; z += o; return *this; }
	Vector3 operator *= (T o) { x *= o; y *= o; z *= o; return *this; }
	Vector3 operator /= (T o) { x /= o; y /= o; z /= o; return *this; }
	Vector3 operator -= (T o) { x -= o; y -= o; z -= o; return *this; }

	Vector3 xxx() const { return Vector3(x); }
	Vector3 yyy() const { return Vector3(y); }
	Vector3 zzz() const { return Vector3(z); }

	static constexpr Vector3 Zero()    { return Vector3(0.0, 0.0, 0.0); }
	static constexpr Vector3 One()     { return Vector3(1.0, 1.0, 1.0); }
	static constexpr Vector3 Up()      { return Vector3(0.0, 1.0, 0.0); }
	static constexpr Vector3 Left()    { return Vector3(-1.0, 0.0, 0.0); }
	static constexpr Vector3 Down()    { return Vector3(0.0, -1.0, 0.0); }
	static constexpr Vector3 Right()   { return Vector3(1.0, 0.0, 0.0); }
	static constexpr Vector3 Forward() { return Vector3(0.0, 0.0, 1.0); }
	static constexpr Vector3 Backward(){ return Vector3(0.0, 0.0, -1.0); }
};

using Vector2d = Vector2<double>;
using Vector2f = Vector2<float>;
using Vector2i = Vector2<int>;
using Vector2s = Vector2<short>;
using Vector2c = Vector2<char>;

using Vector3d = Vector3<double>;
using Vector3f = Vector3<float>;
using Vector3i = Vector3<int>;
using Vector3s = Vector3<short>;
using Vector3c = Vector3<char>;

typedef Vector3f float3;
typedef Vector2f float2;

// recommended to use simd instructions instead. this functions are slow in hot loops
template<typename T> FINLINE Vector3<T> Min(const Vector3<T>& a, const Vector3<T>& b) { return Vector3<T>(Min(a.x, b.x), Min(a.y, b.y), Min(a.z, b.z)); }
template<typename T> FINLINE Vector3<T> Max(const Vector3<T>& a, const Vector3<T>& b) { return Vector3<T>(Max(a.x, b.x), Max(a.y, b.y), Max(a.z, b.z)); }
template<typename T> FINLINE Vector2<T> Min(const Vector2<T>& a, const Vector2<T>& b) { return Vector2<T>(Min(a.x, b.x), Min(a.y, b.y)); }
template<typename T> FINLINE Vector2<T> Max(const Vector2<T>& a, const Vector2<T>& b) { return Vector2<T>(Max(a.x, b.x), Max(a.y, b.y)); }

template<typename T> FINLINE T Max3(const Vector3<T>& a) { return Max(Max(a.x, a.y), a.z); }
template<typename T> FINLINE T Min3(const Vector3<T>& a) { return Min(Min(a.x, a.y), a.z); }

FINLINE Vector2f ToVector2f(const Vector2i& vec) { return Vector2f((float)vec.x, (float)vec.y); }
FINLINE Vector2i ToVector2i(const Vector2f& vec) { return Vector2i((int)vec.x, (int)vec.y);  }

FINLINE uint64 VecToHash(Vector2c vec) { return uint64(vec.x*3 ^vec.y) | (uint64(vec.y) << 8ull ); }
FINLINE uint64 VecToHash(Vector2s vec) { return (uint64)WangHash(uint64(vec.x) | (uint64(vec.y) << 16ull)); }
FINLINE uint64 VecToHash(Vector2i vec) { return MurmurHash(uint64(vec.x) | (uint64(vec.y) << 32ull)); }
FINLINE uint64 VecToHash(Vector3c vec) { return (uint64)WangHash(uint64(vec.x) | (uint64(vec.y) << 8ull) | (uint64(vec.z) << 16ull)); }

FINLINE uint64 VecToHash(Vector3s vec){
	return WangHash(uint64(vec.x) | (uint64(vec.y) << 16ull) | (uint64(vec.z) << 24ull));
}

FINLINE uint64 VecToHash(Vector3i vec) {
	return MurmurHash(uint64(vec.x) | (uint64(vec.y) << 32ull)) + WangHash(vec.z);
}

struct Ray
{
	Vector3f origin;
	Vector3f direction;
	Ray() {}
	Ray(Vector3f o, Vector3f d) : origin(o), direction(d) {}
};

//   ###   VECTOR4   ###

#include "SIMDCommon.hpp"

#ifdef AX_SUPPORT_SSE

AX_ALIGNED(16) struct Vector4f
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

	constexpr Vector4f() : x(0), y(0), z(0), w(0) {}
	constexpr Vector4f(__m128 _vec) : vec(_vec) {}
	constexpr Vector4f(float scale) : x(scale), y(scale), z(scale), w(scale) {}
	constexpr Vector4f(float _x, float _y, float _z, float _w) : x(_x), y(_y), z(_z), w(_w) {}
	constexpr Vector4f(const Vector3f& a, float f) : x(a.x), y(a.y), z(a.z), w(f) {}

	// root of evil
	operator __m128() const { return vec; }

	FINLINE static __m128 VECTORCALL Normalize(const __m128 V)
	{
		return _mm_mul_ps(_mm_rsqrt_ps(_mm_dp_ps(V, V, 0xff)), V);
	}

	FINLINE static __m128 VECTORCALL Dot(const __m128 V1, const __m128 V2)
	{
		return _mm_dp_ps(V1, V2, 0xff);
	}

	float Length() const { return Sqrt(_mm_cvtss_f32(_mm_dp_ps(vec, vec, 0xff))); }

	Vector4f& Normalized() { vec = _mm_mul_ps(_mm_rsqrt_ps(_mm_dp_ps(vec, vec, 0xff)), vec); return *this; }

	Vector3f xyz() const { Vector3f v;  SSEStoreVector3(&v.x, vec); return v; }

	Vector4f VECTORCALL operator + (const Vector4f b) const { return _mm_add_ps(vec, b.vec); }
	Vector4f VECTORCALL operator - (const Vector4f b) const { return _mm_sub_ps(vec, b.vec); }
	Vector4f VECTORCALL operator * (const Vector4f b) const { return _mm_mul_ps(vec, b.vec); }
	Vector4f VECTORCALL operator / (const Vector4f b) const { return _mm_div_ps(vec, b.vec); }

	Vector4f& VECTORCALL operator += (const Vector4f b) { vec = _mm_add_ps(vec, b.vec); return *this; }
	Vector4f& VECTORCALL operator -= (const Vector4f b) { vec = _mm_sub_ps(vec, b.vec); return *this; }
	Vector4f& VECTORCALL operator *= (const Vector4f b) { vec = _mm_mul_ps(vec, b.vec); return *this; }
	Vector4f& VECTORCALL operator /= (const Vector4f b) { vec = _mm_div_ps(vec, b.vec); return *this; }

	Vector4f operator  *  (const float b) const { return _mm_mul_ps(vec, _mm_set_ps1(b)); }
	Vector4f operator  /  (const float b) const { return _mm_div_ps(vec, _mm_set_ps1(b)); }
	Vector4f& operator *= (const float b) { vec = _mm_mul_ps(vec, _mm_set_ps1(b)); return *this; }
	Vector4f& operator /= (const float b) { vec = _mm_div_ps(vec, _mm_set_ps1(b)); return *this; }
};

AX_ALIGNED(16) struct RaySSE
{
	__m128 origin;
	__m128 direction;
};

#endif // AX_SUPPORT_SSE

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

	double Length() const { return _mm256_cvtsd_f64(_mm256_sqrt_pd(Dot(this->vec, this->vec))); }

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

#endif // AX_SUPPORT_AVX2

#if !defined(AX_SUPPORT_SSE ) || !defined(AX_SUPPORT_AVX2)
template<typename T>
struct Vector4
{
	union
	{
		struct { float x, y, z, w; };
		float arr[4];
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

	FINLINE static Vector4 Normalize(const Vector4& V)
	{
		return V * RSqrt(Dot(V, V));
	}

	FINLINE static Vector4 VECTORCALL Dot(const Vector4& V1, const Vector4& V2)
	{
		return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
	}

	float Length() const { return Sqrt(Dot(this)); }

	Vector4& Normalized() { Vector4 norm = Normalize(this); x = norm.x; y = norm.y; z = norm.z; w = norm.w; return *this; }

	Vector3f xyz() const { return Vector3<T>(x, y, z); }

	Vector4 operator - () { return Vector4(-x, -y, -z, -w); 
	Vector4 operator + (const Vector4& other) const { return Vector4(x + other.x, y + other.y, z + other.z, w + other.w); }
	Vector4 operator * (const Vector4& other) const { return Vector4(x * other.x, y * other.y, z * other.z, w * other.w); }
	Vector4 operator / (const Vector4& other) const { return Vector4(x / other.x, y / other.y, z / other.z, w / other.w); }
	Vector4 operator - (const Vector4& other) const { return Vector4(x - other.x, y - other.y, z - other.z, w - other.w); }

	Vector4 operator + (T other) const { return Vector4(x + other, y + other, z + other, w + other); }
	Vector4 operator * (T other) const { return Vector4(x * other, y * other, z * other, w * other); }
	Vector4 operator / (T other) const { return Vector4(x / other, y / other, z / other, w / other); }
	Vector4 operator - (T other) const { return Vector4(x - other, y - other, z - other, w - other); }

	Vector4 operator += (const Vector4& o) { x += o.x; y += o.y; z += o.z; w += o.w; return *this; }
	Vector4 operator *= (const Vector4& o) { x *= o.x; y *= o.y; z *= o.z; w *= o.w; return *this; }
	Vector4 operator /= (const Vector4& o) { x /= o.x; y /= o.y; z /= o.z; w /= o.w; return *this; }
	Vector4 operator -= (const Vector4& o) { x -= o.x; y -= o.y; z -= o.z; w -= o.w; return *this; }

	Vector4 operator += (T o) { x += o; y += o; z += o; w += o; return *this; }
	Vector4 operator *= (T o) { x *= o; y *= o; z *= o; w *= o; return *this; }
	Vector4 operator /= (T o) { x /= o; y /= o; z /= o; w /= o; return *this; }
	Vector4 operator -= (T o) { x -= o; y -= o; z -= o; w -= o; return *this; }
	};

#if !AX_SUPPORT_SSE
	using Vector4f = Vector4<float>;
#endif
#if !AX_SUPPORT_AVX2
	using Vector4d = Vector4<double>;
#endif // !AX_SUPPORT_SSE

#else
	typedef Vector4f float4;
#endif