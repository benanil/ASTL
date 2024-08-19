
/*****************************************************************
*   Purpose:                                                     *
*      Vector2x and Vector3x types and functions.                *
*      where x is: float, double, int, short, half               *
*      vectors doesn't have constructors for performance,        *
*      call MakeVecX to create vector.                           *
*   Be Aware:                                                    *
*      Vector4x is not exist, because we want to have-           *
*      better use of hardware, use vec_t or veci_t instead.      *
*      Declared in SIMDVectorMath.hpp                            *
*   Author : Anilcan Gulkaya 2023 anilcangulkaya7@gmail.com      *
*****************************************************************/

#pragma once

#include "Math.hpp"

AX_NAMESPACE

template<typename T>
struct Vector2
{
    union
    {
        struct { T x, y; };
        struct { T arr[2]; };
    };
    
    static const int NumElements = 2;
    static const uint64 ElementSize = sizeof(T);
    using ElemType = T;
    
    float Length()        const { return Sqrt(LengthSquared()); }
    float LengthSafe()    const { float sq = LengthSquared(); return sq < 0.01f ? 0.0f : Sqrt(sq); }
    float LengthSquared() const { return x * x + y * y; }
    
    static float Length(Vector2 v)        { return Sqrt(v.LengthSquared()); }
    static float LengthSquared(Vector2 v) { return v.x * v.x + v.y * v.y; }
    
    T& operator[] (int index) { return arr[index]; }
    T  operator[] (int index) const { return arr[index]; }
    
    static Vector2 Lerp(const Vector2& a, const Vector2& b, float t)
    {
        Vector2 v;
        v.x = a.x + (b.x - a.x) * t;
        v.y = a.y + (b.y - a.y) * t;
        return v;
    }
    
    purefn static float Distance(Vector2 a, Vector2 b) {
        float diffx = (float)(a.x - b.x);
        float diffy = (float)(a.y - b.y);
        return Sqrt(diffx * diffx + diffy * diffy);
    }
    
    purefn static float DistanceSq(Vector2 a, Vector2 b) {
        float diffx = (float)(a.x - b.x);
        float diffy = (float)(a.y - b.y);
        return diffx * diffx + diffy * diffy;
    }
    
    purefn static Vector2 Rotate(Vector2 vec, float angle)
    {
        float s = Sin(angle), c = Cos(angle);
        return Vector2(vec.x * c - vec.y * s, vec.x * s + vec.y * c);
    }
    
    Vector2& Normalized() { *this /= Length(); return *this; }
    void NormalizeSelf()  { *this /= Length(); }
    static Vector2 Normalize(Vector2 x) { return x / x.Length(); }
    static Vector2 NormalizeEst(const Vector2& a) {
        return a * RSqrt(a.x * a.x + a.y * a.y);
    }
    
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

template<typename T> purefn Vector2<T> MakeVec2(T scale = 0) { Vector2<T> v; v.x = v.y = scale; return v; }
template<typename T> purefn Vector2<T> MakeVec2(T a, T b) { Vector2<T> v; v.x = a; v.y = b;  return v; }

template<typename T>
struct Vector3
{
    union
    {
        struct { T x, y, z; };
        struct { T arr[3]; };
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
    
    purefn static float Dot(const Vector3& a, const Vector3& b)
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
    
    static Vector3 Normalize(const Vector3& a) {
        return a / Sqrt(Vector3::Dot(a, a));
    }
    
    static Vector3 NormalizeEst(const Vector3& a) {
        return a * RSqrt(Vector3::Dot(a, a));
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

template<typename T> purefn Vector3<T> MakeVec3()              { Vector3<T> v; v.x = v.y = v.z = T(0);    return v; }
template<typename T> purefn Vector3<T> MakeVec3(T scale)       { Vector3<T> v; v.x = v.y = v.z = scale;   return v; }
template<typename T> purefn Vector3<T> MakeVec3(T a, T b, T c) { Vector3<T> v; v.x = a; v.y = b; v.z = c; return v; }
template<typename T> purefn Vector3<T> MakeVec3(T* p) { Vector3<T> v; v.x = p[0]; v.y = p[1]; v.z = p[2]; return v; }

using Vector2d = Vector2<double>;
using Vector2f = Vector2<float>;
using Vector2i = Vector2<int>;
using Vector2s = Vector2<short>;

using Vector3d = Vector3<double>;
using Vector3f = Vector3<float>;
using Vector3i = Vector3<int>;
using Vector3s = Vector3<short>;

typedef Vector3f float3;
typedef Vector2f float2;

struct Rect2D
{
    Vector2f pos;
    Vector2f size;
};

purefn Vector2f Normalize(Vector2f v) {
	return v / Sqrt(v.x * v.x + v.y * v.y);
}

purefn Vector3f Normalize(Vector3f v) {
	return v / Sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

purefn float Dot(Vector3f a, Vector3f b) 
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

// macro min and macro max used for more inlining (even in debug mode)
// recommended to use simd instructions instead. this functions are slow in hot loops
template<typename T> purefn Vector3<T> Min(const Vector3<T>& a, const Vector3<T>& b) { return { MacroMin(a.x, b.x), MacroMin(a.y, b.y), MacroMin(a.z, b.z) }; }
template<typename T> purefn Vector3<T> Max(const Vector3<T>& a, const Vector3<T>& b) { return { MacroMax(a.x, b.x), MacroMax(a.y, b.y), MacroMax(a.z, b.z) }; }
template<typename T> purefn Vector2<T> Min(const Vector2<T>& a, const Vector2<T>& b) { return { MacroMin(a.x, b.x), MacroMin(a.y, b.y)}; }
template<typename T> purefn Vector2<T> Max(const Vector2<T>& a, const Vector2<T>& b) { return { MacroMax(a.x, b.x), MacroMax(a.y, b.y)}; }

template<typename T> purefn T Max3(const Vector3<T>& a) { return MAX(MacroMax(a.x, a.y), a.z); }
template<typename T> purefn T Min3(const Vector3<T>& a) { return MIN(MacroMin(a.x, a.y), a.z); }

purefn Vector2s ToVector2s(const Vector2f& vec) { return {(short)vec.x, (short)vec.y }; }
purefn Vector2f ToVector2f(const Vector2s& vec) { return {(float)vec.x, (float)vec.y }; }
purefn Vector2f ToVector2f(const Vector2i& vec) { return {(float)vec.x, (float)vec.y }; }
purefn Vector2i ToVector2i(const Vector2f& vec) { return { (int)vec.x , (int)vec.y   }; }

pureconst bool PointBoxIntersection(Vector2f min, Vector2f max, Vector2f point)
{
    return point.x <= max.x && point.y <= max.y &&
           point.x >= min.x && point.y >= min.y;
}

pureconst uint32 GreaterThan(Vector2f a, Vector2f b)
{
    return uint32(a.x > b.x) | (uint32(a.y > b.y) << 1);
}

pureconst uint32 LessThan(Vector2f a, Vector2f b)
{
    return uint32(a.x < b.x) | (uint32(a.y < b.y) << 1);
}

pureconst uint32 GreaterThan(Vector3f a, Vector3f b)
{
    return uint32(a.x > b.x) | (uint32(a.y > b.y) << 1) | (uint32(a.y > b.y) << 2);
}

pureconst uint32 LessThan(Vector3f a, Vector3f b)
{
    return uint32(a.x < b.x) | (uint32(a.y < b.y) << 1) | (uint32(a.y < b.y) << 2);
}

pureconst uint32 All2(uint32 msk) { return msk == 0x11u; }
pureconst uint32 All3(uint32 msk) { return msk == 0x111u; }

pureconst uint32 Any2(uint32 msk) { return msk > 0; }
pureconst uint32 Any3(uint32 msk) { return msk > 0; }

// purefn uint64 VecToHash(Vector2s vec) { return (uint64)WangHash(uint64(vec.x) | (uint64(vec.y) << 16ull)); }
// purefn uint64 VecToHash(Vector2i vec) { return MurmurHash(uint64(vec.x) | (uint64(vec.y) << 32ull)); }
// 
// purefn uint64 VecToHash(Vector3s vec){
// 	return WangHash(uint64(vec.x) | (uint64(vec.y) << 16ull) | (uint64(vec.z) << 24ull));
// }
// 
// purefn uint64 VecToHash(Vector3i vec) {
// 	return MurmurHash(uint64(vec.x) | (uint64(vec.y) << 32ull)) + WangHash(vec.z);
// }

AX_END_NAMESPACE 