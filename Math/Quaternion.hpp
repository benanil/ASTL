#pragma once

#include "Vector.hpp"
#include "SIMDVectorMath.hpp"

AX_NAMESPACE

struct alignas(16) Quaternion
{
    union
    {
        struct { float x, y, z, w; };
        struct { float arr[4]; };
        struct { vec_t vec; };
    };
    
    const float  operator [] (int index) const { return arr[index]; }
    float& operator [] (int index) { return arr[index]; }
    
    __forceinline static Quaternion Identity() 
    {
        return {0.0f, 0.0f, 0.0f, 1.0f};
    }
    
    inline static Quaternion Normalize(Quaternion quat)
    {
        quat.vec = VecNorm(quat.vec);
        return quat;
    }
    
    inline static vec_t VECTORCALL Mul(const vec_t Q1, const vec_t Q2) 
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
    
    inline static Vector3f VECTORCALL MulVec3(Vector3f vec, Quaternion quat)
    {
        Vector3f res;
        Vec3Store(&res.x, MulVec3(VecSetR(vec.x, vec.y, vec.z, 1.0f), quat.vec));
        return res;
    }
    
    inline static vec_t VECTORCALL MulVec3(vec_t vec, vec_t quat)
    {
        vec_t temp0 = Vec3Cross(quat, vec);
        vec_t temp1 = VecMul(vec, VecSplatZ(quat));
        temp0 = VecAdd(temp0, temp1);
        temp1 = VecMul(Vec3Cross(quat, temp0), VecSet1(2.0f));
        return VecAdd(vec, temp1);
    }
    
    inline Quaternion VECTORCALL Slerp(const Quaternion Q0, const Quaternion Q1, float t)
    {
        const vec_t T = VecSet1(t);
		// Result = Q0 * sin((1.0 - t) * Omega) / sin(Omega) + Q1 * sin(t * Omega) / sin(Omega)
		static const vec_t OneMINusEpsilon = VecSetR(1.0f - 0.00001f, 1.0f - 0.00001f, 1.0f - 0.00001f, 1.0f - 0.00001f);
		static const veci_t SignMask2 = VecFromInt(0x80000000, 0x00000000, 0x00000000, 0x00000000);

		vec_t CosOmega = VecDot(Q0.vec, Q1.vec);

		const vec_t Zero = VecZero();
		vec_t Control = VecCmpLt(CosOmega, Zero);
		vec_t Sign = VecSelect(VecOne(), VecNegativeOne(), Control);
		CosOmega = VecMul(CosOmega, Sign);
		Control = VecCmpLt(CosOmega, OneMINusEpsilon);

		vec_t SinOmega = VecMul(CosOmega, CosOmega);
		SinOmega = VecSub(VecOne(), SinOmega);
		SinOmega = VecSqrt(SinOmega);

		vec_t Omega = VecAtan2(SinOmega, CosOmega);
		vec_t V01 = VecSwizzle(T, 2, 3, 0, 1);
		V01 = VecAnd(V01, VecMaskXY);
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

		vec_t Result = VecMul(Q0.vec, S0);
		S1 = VecMul(S1, Q1.vec); // _mm_fmadd_ps(S1, Q1.vec, Result) ?
		Quaternion res;
		res.vec = VecAdd(Result, S1);
		return res;
    }
    
    __forceinline Quaternion static FromEuler(float x, float y, float z) noexcept
    {
        x *= 0.5f; y *= 0.5f; z *= 0.5f;
        float c[4], s[4];
        vec_t cv;
        vec_t sv = VecSinCos(&cv, VecSetR(x, y, z, 1.0f));
        VecStore(c, cv);
        VecStore(s, sv);
        
        Quaternion q;
        q.w = c[0] * c[1] * c[2] + s[0] * s[1] * s[2];
        q.x = s[0] * c[1] * c[2] - c[0] * s[1] * s[2];
        q.y = c[0] * s[1] * c[2] + s[0] * c[1] * s[2];
        q.z = c[0] * c[1] * s[2] - s[0] * s[1] * c[2];
        return q;
    }
    
    inline Vector3f static ToEulerAngles(const Quaternion& q) noexcept {
        Vector3f eulerAngles; // using std is recommended
        eulerAngles.x = ATan2(2.0f * (q.y * q.z + q.w * q.x), q.w * q.w - q.x * q.x - q.y * q.y + q.z * q.z);
        eulerAngles.y = ASin(Clamp(-2.0f * (q.x * q.z - q.w * q.y), -1.0f, 1.0f));
        eulerAngles.z = ATan2(2.0f * (q.x * q.y + q.w * q.z), q.w * q.w + q.x * q.x - q.y * q.y - q.z * q.z);
        return eulerAngles;
    }
    
    __forceinline Quaternion static VECTORCALL FromEuler(Vector3f euler) noexcept
    {
    	return FromEuler(euler.x, euler.y, euler.z);
    }
    
    static inline Quaternion FromLookRotation(Vector3f direction, const Vector3f& up)
    {
        Quaternion result;
        Vector3f matrix[3] {
            Vector3f::Cross(up, direction), up, direction  
        };
        QuaternionFromMatrix<3>(&result.x, &matrix[0].x);
        return result;
    }
    
    __forceinline static Quaternion Inverse(Quaternion q)
    {
        const float lengthSq = VecDotf(q.vec, q.vec);
        if (lengthSq == 1.0f)
        {
            q.vec = Conjugate(q.vec);
            return q;
        }
        else if (lengthSq >= 0.001f)
        {
            q.vec = VecMul(Conjugate(q.vec), VecSet1(1.0f / lengthSq));
            return q;
        }
        else
        {
        	return Identity();
        }
    }
    
    __forceinline Quaternion Inversed() const
    {
        return Inverse(*this);
    }
    
    __forceinline static vec_t VECTORCALL Conjugate(const vec_t vec)
    {
        static const vec_t NegativeOne3 = VecSetR(-1.0f,-1.0f,-1.0f, 1.0f);
        return VecMul(vec, NegativeOne3);
    }
    
    __forceinline static Quaternion VECTORCALL Conjugate(const Quaternion vec)
    {
        static const Quaternion NegativeOne3 = {-1.0f, -1.0f, -1.0f, 1.0f};
        return vec * NegativeOne3;
    }
    
    Vector3f GetForward() const {
        Vector3f res;
        Vec3Store(&res.x, MulVec3(VecSetR( 0.0f, 0.0f, -1.0f, 0.0f), Conjugate(vec)));
        return res; 
    }
    
    Vector3f GetRight() const {
        Vector3f res;
        Vec3Store(&res.x, MulVec3(VecSetR( 1.0f, 0.0f, 0.0f, 0.0f), Conjugate(vec)));
        return res; 
    }
    
    Vector3f GetLeft() const {
        Vector3f res;
        Vec3Store(&res.x, MulVec3(VecSetR(-1.0f, 0.0f, 0.0f, 0.0f), Conjugate(vec))); 
        return res; 
    }
    
    Vector3f GetUp() const {
        Vector3f res;
        Vec3Store(&res.x, MulVec3(VecSetR( 0.0f, 1.0f, 0.0f, 0.0f), Conjugate(vec)));
        return res; 
    }
    
    __forceinline Quaternion operator *  (float b)             const { Quaternion q; q.vec = VecMul(this->vec, VecSet1(b)); return q; }
    __forceinline Quaternion operator *= (float b)                   { this->vec = Mul(this->vec, VecSet1(b)); return *this; }
    __forceinline Quaternion operator *  (const Quaternion& b) const { Quaternion q; q.vec = Mul(this->vec, b.vec); return q; }
    __forceinline Quaternion operator *= (const Quaternion& b)       { this->vec = Mul(this->vec, b.vec); return *this; }
    __forceinline Quaternion operator +  (const Quaternion& b) const { Quaternion q; q.vec = VecAdd(vec, b.vec); return q; }
    __forceinline Quaternion operator += (const Quaternion& b)       { vec = VecAdd(vec, b.vec); return *this; }
};

__forceinline Quaternion MakeQuat(float scale = 0.0f)
{
    Quaternion v;
    v.vec = VecSet1(scale);
    return v;
}

__forceinline Quaternion MakeQuat(float _x, float _y, float _z, float _w) 
{
    Quaternion v;
    v.vec = VecSetR(_x, _y, _z, _w);
    return v; 
}

__forceinline Quaternion MakeQuat(const float* _x)
{
    Quaternion v;
    v.vec = VecLoad(_x);
    return v;
}


AX_END_NAMESPACE 