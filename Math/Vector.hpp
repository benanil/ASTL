#pragma once

#include "Math.hpp"
#include "SIMDCommon.hpp"
#include "../Random.hpp" 

AX_NAMESPACE 

template<typename T>
struct Vector2
{
	union
	{
		struct { T x, y; };
		T arr[2];
	};
	
	static const int NumElements = 2;
	static const uint64 ElementSize = sizeof(T);
	using ElemType = T;

	float Length()		  const { return Sqrt(LengthSquared()); }
	float LengthSquared() const { return x * x + y * y; }

	T& operator[] (int index) { return arr[index]; }
	T  operator[] (int index) const { return arr[index]; }

	static __forceinline float Distance(Vector2 a, Vector2 b) {
		float diffx = (float)(a.x - b.x);
		float diffy = (float)(a.y - b.y);
		return Sqrt(diffx * diffx + diffy * diffy);
	}

	static __forceinline float DistanceSq(Vector2 a, Vector2 b) {
		float diffx = (float)(a.x - b.x);
		float diffy = (float)(a.y - b.y);
		return diffx * diffx + diffy * diffy;
	}

	static __forceinline Vector2 Rotate(Vector2 vec, float angle)
	{
		float s = Sin(angle), c = Cos(angle);
		return Vector2(vec.x * c - vec.y * s, vec.x * s + vec.y * c);
	}

	void Normalized() const { *this /= Length(); }
	
	Vector2 Normalize(Vector2 other) { return other.Normalize(); }
	
	Vector2 operator - () { return { -x, -y }; }
	Vector2 operator + (Vector2 other) const { return {x + other.x, y + other.y}; }
	Vector2 operator * (Vector2 other) const { return {x * other.x, y * other.y}; }
	Vector2 operator / (Vector2 other) const { return {x / other.x, y / other.y}; }
	Vector2 operator - (Vector2 other) const { return {x - other.x, y - other.y}; }

	Vector2 operator + (T other) const { return {x + other, y + other}; }
	Vector2 operator * (T other) const { return {x * other, y * other}; }
	Vector2 operator / (T other) const { return {x / other, y / other}; }
	Vector2 operator - (T other) const { return {x - other, y - other}; }
	
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

	static __constexpr Vector2 Zero()     { return {T( 0), T( 0)}; } 
	static __constexpr Vector2 One()      { return {T( 1), T( 1)}; } 
	static __constexpr Vector2 Up()       { return {T( 0), T( 1)}; } 
	static __constexpr Vector2 Left()     { return {T(-1), T( 0)}; } 
	static __constexpr Vector2 Down()     { return {T( 0), T(-1)}; } 
	static __constexpr Vector2 Right()    { return {T( 1), T( 0)}; } 
};

template<typename T> __forceinline Vector2<T> MakeVec2(T scale = 0) { Vector2<T> v; v.x = v.y = scale; return v; }
template<typename T> __forceinline Vector2<T> MakeVec2(T a, T b) { Vector2<T> v; v.x = a; v.y = b;  return v; }

template<typename T>
struct Vector3
{
	union
	{
		struct { T x, y, z; };
		T arr[3];
	};

	static const int NumElements = 3;
	using ElemType = T;

	T& operator[] (int index) { return arr[index]; }
	T  operator[] (int index) const { return arr[index]; }

	float Length() const { return Sqrt(LengthSquared()); }
	__constexpr float LengthSquared() const { return x * x + y * y + z * z; }

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
		Vector3 v;
		v.x = a.x + (b.x - a.x) * t;
		v.y = a.y + (b.y - a.y) * t;
		v.z = a.z + (b.z - a.z) * t;
		return v;
	}

	static Vector3 Cross(const Vector3& a, const Vector3& b)
	{
		Vector3 v;
		v.x = a.y * b.z - b.y * a.z;
		v.y = a.z * b.x - b.z * a.x;
		v.z = a.x * b.y - b.x * a.y;
		return v;
	}

	static Vector3 Reflect(const Vector3& in, const Vector3& normal)
	{
		return in - normal * Vector3::Dot(normal, in) * 2.0f;
	}
	// for more accuracy you can use sqrt instead of rsqrt: a / sqrt(dot(a,a)) 
	static Vector3 Normalize(const Vector3& a) {
		return a / Sqrt(Vector3::Dot(a, a));
	}

	Vector3 operator- () { return { -x, -y, -z }; }
	Vector3 operator + (const Vector3& other) const { return { x + other.x, y + other.y, z + other.z }; }
	Vector3 operator * (const Vector3& other) const { return { x * other.x, y * other.y, z * other.z }; }
	Vector3 operator / (const Vector3& other) const { return { x / other.x, y / other.y, z / other.z }; }
	Vector3 operator - (const Vector3& other) const { return { x - other.x, y - other.y, z - other.z }; }

	Vector3 operator + (T other) const { return { x + other, y + other, z + other }; }
	Vector3 operator * (T other) const { return { x * other, y * other, z * other }; }
	Vector3 operator / (T other) const { return { x / other, y / other, z / other }; }
	Vector3 operator - (T other) const { return { x - other, y - other, z - other }; }

	Vector3 operator += (const Vector3& o) { x += o.x; y += o.y; z += o.z; return *this; }
	Vector3 operator *= (const Vector3& o) { x *= o.x; y *= o.y; z *= o.z; return *this; }
	Vector3 operator /= (const Vector3& o) { x /= o.x; y /= o.y; z /= o.z; return *this; }
	Vector3 operator -= (const Vector3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }

	Vector3 operator += (T o) { x += o; y += o; z += o; return *this; }
	Vector3 operator *= (T o) { x *= o; y *= o; z *= o; return *this; }
	Vector3 operator /= (T o) { x /= o; y /= o; z /= o; return *this; }
	Vector3 operator -= (T o) { x -= o; y -= o; z -= o; return *this; }

	Vector3 xxx() const { return {x, x, x}; }
	Vector3 yyy() const { return {y, y, y}; }
	Vector3 zzz() const { return {z, z, z}; }

	static const Vector3 Zero()    { return { 0.0,  0.0,  0.0}; }
	static const Vector3 One()     { return { 1.0,  1.0,  1.0}; }
	static const Vector3 Up()      { return { 0.0,  1.0,  0.0}; }
	static const Vector3 Left()    { return {-1.0,  0.0,  0.0}; }
	static const Vector3 Down()    { return { 0.0, -1.0,  0.0}; }
	static const Vector3 Right()   { return { 1.0,  0.0,  0.0}; }
	static const Vector3 Forward() { return { 0.0,  0.0,  1.0}; }
	static const Vector3 Backward(){ return { 0.0,  0.0, -1.0}; }
};

template<typename T> __forceinline Vector3<T> MakeVec3()              { Vector3<T> v; v.x = v.y = v.z = T(0);    return v; } 
template<typename T> __forceinline Vector3<T> MakeVec3(T scale)       { Vector3<T> v; v.x = v.y = v.z = scale;   return v; } 
template<typename T> __forceinline Vector3<T> MakeVec3(T a, T b, T c) { Vector3<T> v; v.x = a; v.y = b; v.z = c; return v; }

using Vector2d = Vector2<double>;
using Vector2f = Vector2<float>;
using Vector2i = Vector2<int>;
using Vector2s = Vector2<short>;
using Vector2h = Vector2<half>;
using Vector2c = Vector2<char>;

using Vector3d = Vector3<double>;
using Vector3f = Vector3<float>;
using Vector3i = Vector3<int>;
using Vector3s = Vector3<short>;
using Vector3h = Vector3<half>;
using Vector3c = Vector3<char>;

typedef Vector3f float3;
typedef Vector2f float2;
typedef Vector3s half3;
typedef Vector2s half2;

__forceinline float2 ConvertToFloat2(const half* h) {
	float2 res; ConvertHalfToFloat(&res.x, h, 2); return res;
}

__forceinline float3 ConvertToFloat3(const half* h) {
	float3 res; ConvertHalfToFloat(&res.x, h, 3); return res;
}

__forceinline half2 ConvertToHalf2(const float* h) {
	half2 res; 
	res.x = ConvertFloatToHalf(*h++);
	res.y = ConvertFloatToHalf(*h);
	return res;
}

__forceinline half3 ConvertToHalf3(const float* h) {
	half3 res; 
	res.x = ConvertFloatToHalf(*h++);
	res.y = ConvertFloatToHalf(*h++);
	res.z = ConvertFloatToHalf(*h);
	return res;
}

// recommended to use simd instructions instead. this functions are slow in hot loops
template<typename T> __forceinline Vector3<T> Min(const Vector3<T>& a, const Vector3<T>& b) { return { MIN(a.x, b.x), MIN(a.y, b.y), MIN(a.z, b.z) }; }
template<typename T> __forceinline Vector3<T> Max(const Vector3<T>& a, const Vector3<T>& b) { return { MAX(a.x, b.x), MAX(a.y, b.y), MAX(a.z, b.z) }; }
template<typename T> __forceinline Vector2<T> Min(const Vector2<T>& a, const Vector2<T>& b) { return { MIN(a.x, b.x), MIN(a.y, b.y)}; }
template<typename T> __forceinline Vector2<T> Max(const Vector2<T>& a, const Vector2<T>& b) { return { MAX(a.x, b.x), MAX(a.y, b.y)}; }

template<typename T> __forceinline T Max3(const Vector3<T>& a) { return MAX(MAX(a.x, a.y), a.z); }
template<typename T> __forceinline T Min3(const Vector3<T>& a) { return MIN(MIN(a.x, a.y), a.z); }

__forceinline Vector2f ToVector2f(const Vector2i& vec) { return {(float)vec.x, (float)vec.y }; }
__forceinline Vector2i ToVector2i(const Vector2f& vec) { return { (int)vec.x , (int)vec.y   }; }

__forceinline uint64 VecToHash(Vector2c vec) { return uint64(vec.x*3 ^vec.y) | (uint64(vec.y) << 8ull ); }
__forceinline uint64 VecToHash(Vector2s vec) { return (uint64)WangHash(uint64(vec.x) | (uint64(vec.y) << 16ull)); }
__forceinline uint64 VecToHash(Vector2i vec) { return MurmurHash(uint64(vec.x) | (uint64(vec.y) << 32ull)); }
__forceinline uint64 VecToHash(Vector3c vec) { return (uint64)WangHash(uint64(vec.x) | (uint64(vec.y) << 8ull) | (uint64(vec.z) << 16ull)); }

__forceinline uint64 VecToHash(Vector3s vec){
	return WangHash(uint64(vec.x) | (uint64(vec.y) << 16ull) | (uint64(vec.z) << 24ull));
}

__forceinline uint64 VecToHash(Vector3i vec) {
	return MurmurHash(uint64(vec.x) | (uint64(vec.y) << 32ull)) + WangHash(vec.z);
}

//   ###   VECTOR4   ###

#ifdef AX_SUPPORT_SSE
struct alignas(16) Vector4f
{
	union
	{
		struct { float x, y, z, w; };
		float arr[4];
		__m128 vec;
	};

	static const int NumElements = 4;
	static const uint64 ElementSize = sizeof(float);
	using ElemType = float;

	float& operator[] (int index) { return arr[index]; }
	float  operator[] (int index) const { return arr[index]; }

	__forceinline static __m128 VECTORCALL Normalize(const __m128 V)
	{
		return _mm_mul_ps(_mm_rsqrt_ps(_mm_dp_ps(V, V, 0xff)), V);
	}

	__forceinline static __m128 VECTORCALL Dot(const __m128 V1, const __m128 V2)
	{
		return _mm_dp_ps(V1, V2, 0xff);
	}

	float Length() const { return Sqrt(_mm_cvtss_f32(_mm_dp_ps(vec, vec, 0xff))); }

	Vector4f& Normalized() { vec = _mm_mul_ps(_mm_rsqrt_ps(_mm_dp_ps(vec, vec, 0xff)), vec); return *this; }

	Vector3f xyz() const { Vector3f v;  SSEStoreVector3(&v.x, vec); return v; }

	Vector4f VECTORCALL operator + (const Vector4f b) const { Vector4f v; v.vec = _mm_add_ps(vec, b.vec); return v; }
	Vector4f VECTORCALL operator - (const Vector4f b) const { Vector4f v; v.vec = _mm_sub_ps(vec, b.vec); return v; }
	Vector4f VECTORCALL operator * (const Vector4f b) const { Vector4f v; v.vec = _mm_mul_ps(vec, b.vec); return v; }
	Vector4f VECTORCALL operator / (const Vector4f b) const { Vector4f v; v.vec = _mm_div_ps(vec, b.vec); return v; }

	Vector4f& VECTORCALL operator += (const Vector4f b) { vec = _mm_add_ps(vec, b.vec); return *this; }
	Vector4f& VECTORCALL operator -= (const Vector4f b) { vec = _mm_sub_ps(vec, b.vec); return *this; }
	Vector4f& VECTORCALL operator *= (const Vector4f b) { vec = _mm_mul_ps(vec, b.vec); return *this; }
	Vector4f& VECTORCALL operator /= (const Vector4f b) { vec = _mm_div_ps(vec, b.vec); return *this; }

	Vector4f operator  *  (const float b) const { Vector4f v; v.vec = _mm_mul_ps(vec, _mm_set_ps1(b)); return v; }
	Vector4f operator  /  (const float b) const { Vector4f v; v.vec = _mm_div_ps(vec, _mm_set_ps1(b)); return v; }
	Vector4f& operator *= (const float b) { vec = _mm_mul_ps(vec, _mm_set_ps1(b)); return *this; }
	Vector4f& operator /= (const float b) { vec = _mm_div_ps(vec, _mm_set_ps1(b)); return *this; }
};

__forceinline Vector4f MakeVec4(float scale = 0.0f)                     { Vector4f v; v.vec = _mm_set1_ps(scale);            return v; }
__forceinline Vector4f MakeVec4(float _x, float _y, float _z, float _w) { Vector4f v; v.vec = _mm_setr_ps(_x, _y, _z, _w);   return v; }
__forceinline Vector4f MakeVec4(const Vector3<float>& a, float f)       { Vector4f v; v.vec = _mm_setr_ps(a.x, a.y, a.z, f); return v; }

struct RaySSE
{
	__m128 origin;
	__m128 direction;
};

#endif // AX_SUPPORT_SSE

#ifdef AX_SUPPORT_AVX2
struct alignas(32) Vector4d
{
	union
	{
		struct { double x, y, z, w; };
		double arr[4];
		__m256d vec;
	};

	static const int NumElements = 4;
	static const uint64 ElementSize = sizeof(double);
	using ElemType = double;

	double Length() const { return _mm256_cvtsd_f64(_mm256_sqrt_pd(Dot(this->vec, this->vec))); }

	__forceinline static __m256d VECTORCALL Normalize(__m256d V)
	{
		return _mm256_div_pd(V, _mm256_sqrt_pd(Dot(V, V))); // v / sqrt(dot(v))
	}

	__forceinline static __m256d VECTORCALL Dot(const __m256d V1, const __m256d V2)
	{
		__m256d vDot = _mm256_mul_pd(V1, V2);
		__m256d vTemp = _mm256_permute4x64_pd(vDot, _MM_SHUFFLE(2, 1, 2, 1));
		vDot = _mm256_add_pd(vDot, vTemp);
		vTemp = _mm256_permute4x64_pd(vTemp, _MM_SHUFFLE(1, 1, 1, 1));
		vDot = _mm256_add_pd(vDot, vTemp);
		return _mm256_permute4x64_pd(vDot, _MM_SHUFFLE(0, 0, 0, 0));
	}

	Vector4d& Normalized() { vec = Normalize(vec); return *this; }

	Vector4d VECTORCALL operator + (const Vector4d b) const { Vector4d v; v.vec = _mm256_add_pd(vec, b.vec); return v; }
	Vector4d VECTORCALL operator - (const Vector4d b) const { Vector4d v; v.vec = _mm256_sub_pd(vec, b.vec); return v; }
	Vector4d VECTORCALL operator * (const Vector4d b) const { Vector4d v; v.vec = _mm256_mul_pd(vec, b.vec); return v; }
	Vector4d VECTORCALL operator / (const Vector4d b) const { Vector4d v; v.vec = _mm256_div_pd(vec, b.vec); return v; }

	Vector4d& VECTORCALL operator += (const Vector4d b) { vec = _mm256_add_pd(vec, b.vec); return *this; }
	Vector4d& VECTORCALL operator -= (const Vector4d b) { vec = _mm256_sub_pd(vec, b.vec); return *this; }
	Vector4d& VECTORCALL operator *= (const Vector4d b) { vec = _mm256_mul_pd(vec, b.vec); return *this; }
	Vector4d& VECTORCALL operator /= (const Vector4d b) { vec = _mm256_div_pd(vec, b.vec); return *this; }

	Vector4d operator  *  (const double b) const { Vector4d v; v.vec = _mm256_mul_pd(vec, _mm256_set1_pd(b)); return v; }
	Vector4d operator  /  (const double b) const { Vector4d v; v.vec = _mm256_div_pd(vec, _mm256_set1_pd(b)); return v; }
	Vector4d& operator *= (const double b) { vec = _mm256_mul_pd(vec, _mm256_set1_pd(b)); return *this; }
	Vector4d& operator /= (const double b) { vec = _mm256_div_pd(vec, _mm256_set1_pd(b)); return *this; }
};

__forceinline Vector4d MakeVec4(double scale = 0.0)                         { Vector4d v; v.vec = _mm256_set1_pd(scale);            return v; }
__forceinline Vector4d MakeVec4(double _x, double _y, double _z, double _w) { Vector4d v; v.vec = _mm256_setr_pd(_x, _y, _z, _w);   return v; }
__forceinline Vector4d MakeVec4(const Vector3<double>& a, double f)         { Vector4d v; v.vec = _mm256_setr_pd(a.x, a.y, a.z, f); return v; }

#endif // AX_SUPPORT_AVX2

template<typename T>
struct Vector4
{
	union
	{
		struct { T x, y, z, w; };
		T arr[4];
	};

	static const int NumElements = 4;
	using ElemType = T;

	T& operator[] (int index) { return arr[index]; }
	T  operator[] (int index) const { return arr[index]; }

	float Length() const { return Sqrt(LengthSquared()); }
	__constexpr float LengthSquared() const { return x * x + y * y + z * z + w * w; }

	Vector4& Normalized() { *this /= Length(); return *this; }
	void NormalizeSelf() { *this /= Length(); }

	static float Length(const Vector4& vec) { return vec.Length(); }

	static float Dot(const Vector4& a, const Vector4& b)
	{
		return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
	}

	static Vector4 Lerp(const Vector4& a, const Vector4& b, float t)
	{
		Vector4 v;
		v.x = a.x + (b.x - a.x) * t;
		v.y = a.y + (b.y - a.y) * t;
		v.z = a.z + (b.z - a.z) * t;
		v.w = a.w + (b.w - a.w) * t;
		return v;
	}

	// for more accuracy you can use sqrt instead of rsqrt: a / sqrt(dot(a,a)) 
	static Vector4 Normalize(const Vector4& a) {
		return a / Sqrt(Vector4::Dot(a, a));
	}

	Vector4 operator- () { return { -x, -y, -z }; }
	Vector4 operator + (const Vector4& other) const { return { x + other.x, y + other.y, z + other.z, w + other.w }; }
	Vector4 operator * (const Vector4& other) const { return { x * other.x, y * other.y, z * other.z, w * other.w }; }
	Vector4 operator / (const Vector4& other) const { return { x / other.x, y / other.y, z / other.z, w / other.w }; }
	Vector4 operator - (const Vector4& other) const { return { x - other.x, y - other.y, z - other.z, w - other.w }; }

	Vector4 operator + (T other) const { return { x + other, y + other, z + other, w + other }; }
	Vector4 operator * (T other) const { return { x * other, y * other, z * other, w * other }; }
	Vector4 operator / (T other) const { return { x / other, y / other, z / other, w / other }; }
	Vector4 operator - (T other) const { return { x - other, y - other, z - other, w - other }; }

	Vector4 operator += (const Vector4& o) { x += o.x; y += o.y; z += o.z; w += o.w; return *this; }
	Vector4 operator *= (const Vector4& o) { x *= o.x; y *= o.y; z *= o.z; w *= o.w; return *this; }
	Vector4 operator /= (const Vector4& o) { x /= o.x; y /= o.y; z /= o.z; w /= o.w; return *this; }
	Vector4 operator -= (const Vector4& o) { x -= o.x; y -= o.y; z -= o.z; w -= o.w; return *this; }

	Vector4 operator += (T o) { x += o; y += o; z += o; w += o; return *this; }
	Vector4 operator *= (T o) { x *= o; y *= o; z *= o; w *= o; return *this; }
	Vector4 operator /= (T o) { x /= o; y /= o; z /= o; w /= o; return *this; }
	Vector4 operator -= (T o) { x -= o; y -= o; z -= o; w -= o; return *this; }
};
	
template<typename T> __forceinline Vector4<T> MakeVec4(T scale = 0)              { Vector4<T> v; v.x = v.y = v.z = v.w = scale;              return v;}
template<typename T> __forceinline Vector4<T> MakeVec4(T _x, T _y, T _z, T _w)   { Vector4<T> v; v.x = _x; v.y = _y; v.z = _z; v.w = _w;     return v;}
template<typename T> __forceinline Vector4<T> MakeVec4(const Vector3<T>& a, T f) { Vector4<T> v; v.x = a.x; v.y = a.y; v.z = a.z; v.w = a.w; return v; }

#ifndef AX_SUPPORT_SSE
	using Vector4f = Vector4<float>;
#endif
#ifndef AX_SUPPORT_AVX2
	using Vector4d = Vector4<double>;
#endif // !AX_SUPPORT_SSE

struct Ray
{
	Vector3f origin;
	Vector3f direction;
	Ray() {}
	Ray(Vector3f o, Vector3f d) : origin(o), direction(d) {}
};

AX_END_NAMESPACE 