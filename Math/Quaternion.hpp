#pragma once

#include "Vector.hpp"
#include "SIMDVectorMath.hpp"

AX_NAMESPACE

typedef Vector4x32f Quaternion;
struct xyzw 
{ 
    float x, y, z, w; 
};

#define QIdentity()  VecSetR(0.0f, 0.0f, 0.0f, 1.0f)
#define QNorm(q)     VecNorm(q)
#define QNormEst(q)  VecNormEst(q)
#define MakeQuat(_x, _y, _z, _w)  VecSetR(_x, _y, _z, _w)

inline Vector4x32f VECTORCALL QMul(Vector4x32f Q1, Vector4x32f Q2) 
{
    const Vector4x32f ControlWZYX = { 1.0f,-1.0f, 1.0f,-1.0f };
    const Vector4x32f ControlZWXY = { 1.0f, 1.0f,-1.0f,-1.0f };
    const Vector4x32f ControlYXWZ = { -1.0f, 1.0f, 1.0f,-1.0f };
    
    Vector4x32f vResult = VecSplatW(Q2);
    Vector4x32f Q2X     = VecSplatX(Q2);
    Vector4x32f Q2Y     = VecSplatY(Q2);
    Vector4x32f Q2Z     = VecSplatZ(Q2);
    vResult = VecMul(vResult, Q1);

    Vector4x32f Q1Shuffle = Q1;
    Q1Shuffle = VecRev(Q1Shuffle);
    Q2X       = VecMul(Q2X, Q1Shuffle);
    Q1Shuffle = VecShuffleR(Q1Shuffle, Q1Shuffle, 2, 3, 0, 1);
    Q2X       = VecMul(Q2X, ControlWZYX);
    Q2Y       = VecMul(Q2Y, Q1Shuffle);
    Q1Shuffle = VecRev(Q1Shuffle);
    Q2Y       = VecMul(Q2Y, ControlZWXY);
    Q2Z       = VecMul(Q2Z, Q1Shuffle);
    vResult   = VecAdd(vResult, Q2X);
    Q2Z       = VecMul(Q2Z, ControlYXWZ);
    Q2Y       = VecAdd(Q2Y, Q2Z);
    vResult   = VecAdd(vResult, Q2Y);
    return vResult;
}

// Angle should be between -twopi, twopi
inline Vector4x32f QFromAxisAngle(Vector3f axis, float angle)
{
    float SinV = Sin(0.5f * angle);
    float CosV = Cos(0.5f * angle);
    Vector4x32f q = VecSetR(axis.x * SinV, axis.y * SinV, axis.z * SinV, CosV);
    return QNorm(q);
}

// below 3 function are same as QFromAxisAngle but with single axis, 
// faster because no normalization and less multipication
inline Vector4x32f QFromXAngle(float angle) {
    return VecSetR(Sin(0.5f * angle), 0.0f, 0.0f, Cos(0.5f * angle));
}

inline Vector4x32f QFromYAngle(float angle) {
    return VecSetR(0.0f, Sin(0.5f * angle), 0.0f, Cos(0.5f * angle));
}

inline Vector4x32f QFromZAngle(float angle) {
    return VecSetR(0.0f, 0.0f, Sin(0.5f * angle), Cos(0.5f * angle));
}

inline Vector4x32f VECTORCALL QMulVec3(Vector4x32f vec, Vector4x32f quat)
{
    Vector4x32f temp0 = Vec3Cross(quat, vec);
    Vector4x32f temp1 = VecMul(vec, VecSplatZ(quat));
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

inline Quaternion VECTORCALL QSlerp(Quaternion q0, Quaternion q1, float t)
{
    const Vector4x32f one = VecSet1(1.0f);
    // from paper: "A Fast and Accurate Estimate for SLERP" by David Eberly
    // but I have used fused instructions and I've made optimizations on sign part for ARM cpu's

    // Common code for computing the scalar coefficients of SLERP
    auto CalculateCoefficient = [one] (Vector4x32f vT, Vector4x32f xm1)
    {
        constexpr float const mu = 1.85298109240830f;
        // Precomputed constants
        const Vector4x32f u0123 = VecSetR( 1.f / ( 1 * 3 ), 1.f / ( 2 * 5 ), 1.f / ( 3 * 7 ), 1.f / ( 4 * 9 ) );
        const Vector4x32f u4567 = VecSetR( 1.f / ( 5 * 11 ), 1.f / ( 6 * 13 ), 1.f / ( 7 * 15 ), mu / ( 8 * 17 ) );
        const Vector4x32f v0123 = VecSetR( 1.f / 3, 2.f / 5, 3.f / 7, 4.f / 9 );
        const Vector4x32f v4567 = VecSetR( 5.f / 11, 6.f / 13, 7.f / 15, mu * 8 / 17 );

        Vector4x32f vTSquared = VecMul(vT, vT);
        Vector4x32f b4567 = VecFmsub(u4567, vTSquared, v4567);
        b4567 = VecMul(b4567, xm1);

        Vector4x32f c = VecAdd(VecSplatW(b4567), one);
        c = VecFmaddLane(c, b4567, one, 2); // multiply by lane is faster with ARM cpu's
        c = VecFmaddLane(c, b4567, one, 1);
        c = VecFmaddLane(c, b4567, one, 0);

        Vector4x32f b0123 = VecFmsub(u0123, vTSquared, v0123);
        b0123 = VecMul(b0123, xm1);
        c = VecFmaddLane(c, b0123, one, 3);
        c = VecFmaddLane(c, b0123, one, 2);
        c = VecFmaddLane(c, b0123, one, 1);
        c = VecFmaddLane(c, b0123, one, 0);
        c = VecMul(c, vT);
        return c;
    };

    Vector4x32f x = VecDot(q0, q1); // cos ( theta ) in all components
    Vector4x32f control = VecCmpLt(x, VecZero());
    Vector4x32f sign = VecSelect(VecOne(), VecNegativeOne(), control);
    q1 = VecMul(sign, q1); // do mul instead of xor

    Vector4x32f xm1 = VecFmsub(x, sign, one);
    Vector4x32f cT = CalculateCoefficient(VecSet1(t), xm1);
    Vector4x32f cD = CalculateCoefficient(VecSet1(1.0f - t), xm1);
    cT = VecMul(cT, q1);
    return VecFmadd(cD, q0, cT);
}

// faster but less precise, more error prone version of slerp
purefn Quaternion VECTORCALL QNLerp(Quaternion a, Quaternion b, float t)
{
    Vector4x32i lz = VecCmpLt(VecDot(a, b), VecZero());
    a = VecSelect(a, VecNeg(a), lz);
    a = VecLerp(a, b, t);
    return VecNormEst(a);
}

purefn Quaternion QFromEuler(float x, float y, float z)
{
    x *= 0.5f; y *= 0.5f; z *= 0.5f;
    float c[4], s[4];
    Vector4x32f cv;
    Vector4x32f sv = VecSinCos(&cv, VecSetR(x, y, z, 1.0f));
    VecStore(c, cv);
    VecStore(s, sv);
    
    Quaternion q = VecSetR(
    s[0] * c[1] * c[2] - c[0] * s[1] * s[2],
    c[0] * s[1] * c[2] + s[0] * c[1] * s[2],
    c[0] * c[1] * s[2] - s[0] * s[1] * c[2],
    c[0] * c[1] * c[2] + s[0] * s[1] * s[2]);
    return q;
}

purefn Quaternion VECTORCALL QFromEuler(Vector3f euler)
{
    return QFromEuler(euler.x, euler.y, euler.z);
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

template<int numCol = 4> // number of columns of matrix, 3 or 4
inline void QuaternionFromMatrix(float* Orientation, const float* m) {
    int i, j, k = 0;
    float root, trace = m[0*numCol+0] + m[1 * numCol + 1] + m[2 * numCol + 2];
    
    if (trace > 0.0f)
    {
        root = RSqrt(trace + 1.0f);
        Orientation[3] = 0.5f / root;
        root = 0.5f * root;
        Orientation[0] = root * (m[1 * numCol + 2] - m[2 * numCol + 1]);
        Orientation[1] = root * (m[2 * numCol + 0] - m[0 * numCol + 2]);
        Orientation[2] = root * (m[0 * numCol + 1] - m[1 * numCol + 0]);
    }
    else
    {
        static const int Next[3] = { 1, 2, 0 };
        i = 0;
        i += m[1 * numCol + 1] > m[0 * numCol + 0]; // if (M.m[1][1] > M.m[0][0]) i = 1
        if (m[2 * numCol + 2] > m[i * numCol + i]) i = 2;
        j = Next[i];
        k = Next[j];
        
        root = RSqrt(m[i * numCol + i] - m[j * numCol + j] - m[k * numCol + k] + 1.0f);
        
        Orientation[i] = 0.5f / root;
        root = 0.5f * root;
        Orientation[j] = root * (m[i * numCol + j] + m[j * numCol + i]);
        Orientation[k] = root * (m[i * numCol + k] + m[k * numCol + i]);
        Orientation[3] = root * (m[j * numCol + k] - m[k*numCol+j]);
    } 
}

template<int numCol = 4> // number of columns of matrix, 3 or 4
void MatrixFromQuaternion(float* mat, Quaternion quat)
{
    xyzw q;
    VecStore(&q.x, quat);
    const float num9 = q.x * q.x, num8 = q.y * q.y,
                num7 = q.z * q.z, num6 = q.x * q.y,
                num5 = q.z * q.w, num4 = q.z * q.x,
                num3 = q.y * q.w, num2 = q.y * q.z,
                num  = q.x * q.w;

    mat[numCol * 0 + 0] = 1.0f - (2.0f * (num8 + num7));
    mat[numCol * 0 + 1] = 2.0f * (num6 + num5);
    mat[numCol * 0 + 2] = 2.0f * (num4 - num3);
    
    mat[numCol * 1 + 0] = 2.0f * (num6 - num5);
    mat[numCol * 1 + 1] = 1.0f - (2.0f * (num7 + num9));
    mat[numCol * 1 + 2] = 2.0f * (num2 + num);
    
    mat[numCol * 2 + 0] = 2.0f * (num4 + num3);
    mat[numCol * 2 + 1] = 2.0f * (num2 - num);
    mat[numCol * 2 + 2] = 1.0f - (2.0f * (num8 + num9));
    
    if_constexpr(numCol == 4)
        mat[numCol * 3 + 3] = 1.0f;
}

inline Quaternion QFromLookRotation(Vector3f direction, const Vector3f& up)
{
    Vector3f matrix[3] = {
        Vector3f::Cross(up, direction), up, direction 
    };
    xyzw result;
    QuaternionFromMatrix<3>(&result.x, &matrix[0].x);
    return VecLoad(&result.x);
}

purefn Quaternion VECTORCALL QConjugate(Quaternion vec)
{
    return VecMul(vec, VecSetR(-1.0f, -1.0f, -1.0f, 1.0f));
}

inline Quaternion QInverse(Quaternion q)
{
    const float lengthSq = VecDotf(q, q);
    if (AlmostEqual(lengthSq, 1.0f))
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
    Vec3Store(&res.x, QMulVec3(VecSetR( 0.0f, 0.0f, 1.0f, 0.0f), QConjugate(vec)));
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