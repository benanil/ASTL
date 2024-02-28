#pragma once

#include "Vector.hpp"
#include "SIMDVectorMath.hpp"

AX_NAMESPACE

typedef vec_t Quaternion;
struct xyzw { float x, y, z, w; };

#define QIdentity()  VecSetR(0.0f, 0.0f, 0.0f, 1.0f)
#define QNorm(q)     VecNorm(q)
#define MakeQuat(_x, _y, _z, _w)  VecSetR(_x, _y, _z, _w)

inline vec_t VECTORCALL QMul(vec_t Q1, vec_t Q2) 
{
    static const vec_t ControlWZYX = { 1.0f,-1.0f, 1.0f,-1.0f };
    static const vec_t ControlZWXY = { 1.0f, 1.0f,-1.0f,-1.0f };
    static const vec_t ControlYXWZ = { -1.0f, 1.0f, 1.0f,-1.0f };
    vec_t Q2X = Q2, Q2Y = Q2, Q2Z = Q2, vResult = Q2;
    vResult = VecSplatW(vResult);
    Q2X     = VecSplatX(Q2X);
    Q2Y     = VecSplatY(Q2Y);
    Q2Z     = VecSplatZ(Q2Z);
    vResult = VecMul(vResult, Q1);

    vec_t Q1Shuffle = Q1;
    Q1Shuffle = VecShuffleR(Q1Shuffle, Q1Shuffle, 0, 1, 2, 3);
    Q2X       = VecMul(Q2X, Q1Shuffle);
    Q1Shuffle = VecShuffleR(Q1Shuffle, Q1Shuffle, 2, 3, 0, 1);
    Q2X       = VecMul(Q2X, ControlWZYX);
    Q2Y       = VecMul(Q2Y, Q1Shuffle);
    Q1Shuffle = VecShuffleR(Q1Shuffle, Q1Shuffle, 0, 1, 2, 3);
    Q2Y       = VecMul(Q2Y, ControlZWXY);
    Q2Z       = VecMul(Q2Z, Q1Shuffle);
    vResult   = VecAdd(vResult, Q2X);
    Q2Z       = VecMul(Q2Z, ControlYXWZ);
    Q2Y       = VecAdd(Q2Y, Q2Z);
    vResult   = VecAdd(vResult, Q2Y);
    return vResult;
}

inline vec_t VECTORCALL QMulVec3(vec_t vec, vec_t quat)
{
    vec_t temp0 = Vec3Cross(quat, vec);
    vec_t temp1 = VecMul(vec, VecSplatZ(quat));
    temp0 = VecAdd(temp0, temp1);
    temp1 = VecMul(Vec3Cross(quat, temp0), VecSet1(2.0f));
    return VecAdd(vec, temp1);
}

inline Vector3f VECTORCALL QMulVec3(Vector3f vec, Quaternion quat)
{
    Vector3f res;
    Vec3Store(&res.x, QMulVec3(VecSetR(vec.x, vec.y, vec.z, 1.0f), quat));
    return res;
}

inline Quaternion VECTORCALL QSlerp(Quaternion Q0, Quaternion Q1, float t)
{
    const vec_t T = VecSet1(t);
    // Result = Q0 * sin((1.0 - t) * Omega) / sin(Omega) + Q1 * sin(t * Omega) / sin(Omega)
    static const vec_t OneMINusEpsilon = VecSetR(1.0f - 0.00001f, 1.0f - 0.00001f, 1.0f - 0.00001f, 1.0f - 0.00001f);
    static const veci_t SignMask2 = VecFromInt(0x80000000, 0x00000000, 0x00000000, 0x00000000);
    
    vec_t CosOmega = VecDot(Q0, Q1);
    vec_t Control = VecCmpLt(CosOmega, VecZero());
    vec_t Sign = VecSelect(VecOne(), VecNegativeOne(), Control);
    CosOmega = VecMul(CosOmega, Sign);
    Control = VecCmpLt(CosOmega, OneMINusEpsilon);
    
    vec_t SinOmega = VecMul(CosOmega, CosOmega);
    SinOmega = VecSub(VecOne(), SinOmega);
    SinOmega = VecSqrt(SinOmega);
    
    vec_t Omega = VecAtan2(SinOmega, CosOmega);
    vec_t V01 = VecSwizzle(T, 2, 3, 0, 1);
    V01 = VecMask(V01, VecMaskXY);
    #if !defined(AX_SUPPORT_SSE) && !defined(AX_ARM)
    uint uf = BitCast<uint>(V01.x) ^ 0x80000000;
    V01.x = BitCast<float>(uf);
    #else
    V01 = VecXor(V01, SignMask2);
    #endif
    V01 = VecAdd(VecIdentityR0, V01);
    
    vec_t S0 = VecMul(V01, Omega);
    S0 = VecSin(S0);
    S0 = VecDiv(S0, SinOmega);
    S0 = VecSelect(V01, S0, Control);
    
    vec_t S1 = VecSplatY(S0);
    S0 = VecSplatX(S0);
    S1 = VecMul(S1, Sign);
    
    vec_t Result = VecMul(Q0, S0);
    S1 = VecMul(S1, Q1); // _mm_fmadd_ps(S1, Q1.vec, Result) ?
    return VecAdd(Result, S1);
}

__forceinline Quaternion static QFromEuler(float x, float y, float z)
{
    x *= 0.5f; y *= 0.5f; z *= 0.5f;
    float c[4], s[4];
    vec_t cv;
    vec_t sv = VecSinCos(&cv, VecSetR(x, y, z, 1.0f));
    VecStore(c, cv);
    VecStore(s, sv);
    
    Quaternion q = VecSet(
    s[0] * c[1] * c[2] - c[0] * s[1] * s[2],
    c[0] * s[1] * c[2] + s[0] * c[1] * s[2],
    c[0] * c[1] * s[2] - s[0] * s[1] * c[2],
    c[0] * c[1] * c[2] + s[0] * s[1] * s[2]);
    return q;
}

inline Vector3f QToEulerAngles(Quaternion qu)
{
    xyzw q;
    VecStore(&q.x, qu);
    Vector3f eulerAngles; // using cstd for trigonometric functions recommended
    eulerAngles.x = ATan2(2.0f * (q.y * q.z + q.w * q.x), q.w * q.w - q.x * q.x - q.y * q.y + q.z * q.z);
    eulerAngles.y = ASin(Clamp(-2.0f * (q.x * q.z - q.w * q.y), -1.0f, 1.0f));
    eulerAngles.z = ATan2(2.0f * (q.x * q.y + q.w * q.z), q.w * q.w + q.x * q.x - q.y * q.y - q.z * q.z);
    return eulerAngles;
}

__forceinline Quaternion  VECTORCALL QFromEuler(Vector3f euler)
{
    return QFromEuler(euler.x, euler.y, euler.z);
}

 inline Quaternion FromLookRotation(Vector3f direction, const Vector3f& up)
{
    xyzw result;
    Vector3f matrix[3] {
        Vector3f::Cross(up, direction), up, direction 
    };
    QuaternionFromMatrix<3>(&result.x, &matrix[0].x);
    return VecLoad(&result.x);
}

__forceinline Quaternion VECTORCALL QConjugate(Quaternion vec)
{
    return VecMul(vec, VecSetR(-1.0f, -1.0f, -1.0f, 1.0f));
}

__forceinline static Quaternion QInverse(Quaternion q)
{
    const float lengthSq = VecDotf(q, q);
    if (lengthSq == 1.0f)
    {
        q = QConjugate(q);
        return q;
    }
    else if (lengthSq >= 0.001f)
    {
        q = VecMul(QConjugate(q), VecSet1(1.0f / lengthSq));
        return q;
    }
    else
    {
    	return QIdentity();
    }
}

inline Vector3f VECTORCALL QGetForward(Quaternion vec) {
    Vector3f res;
    Vec3Store(&res.x, QMulVec3(VecSetR( 0.0f, 0.0f, -1.0f, 0.0f), QConjugate(vec)));
    return res; 
}

inline Vector3f VECTORCALL QGetRight(Quaternion vec) {
    Vector3f res;
    Vec3Store(&res.x, QMulVec3(VecSetR( 1.0f, 0.0f, 0.0f, 0.0f), QConjugate(vec)));
    return res; 
}

inline Vector3f VECTORCALL QGetLeft(Quaternion vec) {
    Vector3f res;
    Vec3Store(&res.x, QMulVec3(VecSetR(-1.0f, 0.0f, 0.0f, 0.0f), QConjugate(vec)));
    return res; 
}

inline Vector3f VECTORCALL QGetUp(Quaternion vec) {
    Vector3f res;
    Vec3Store(&res.x, QMulVec3(VecSetR( 0.0f, 1.0f, 0.0f, 0.0f), QConjugate(vec)));
    return res; 
}

AX_END_NAMESPACE 