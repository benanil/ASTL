#pragma once

#include "Vector.hpp"

#ifdef AX_SUPPORT_SSE

AX_NAMESPACE 

struct alignas(16) Quaternion
{
	union
	{
		struct { float x, y, z, w; };
		float arr[4];
		__m128 vec;
	};

	const float  operator [] (int index) const { return arr[index]; }
	float& operator [] (int index) { return arr[index]; }

	__forceinline static Quaternion Identity() 
	{
		Quaternion q;
		q.vec = _mm_setr_ps(0.0f, 0.0f, 0.0f, 1.0f); 
		return q;
	}

	inline static __m128 VECTORCALL Mul(const __m128 Q1, const __m128 Q2) noexcept
	{
		static const __m128 ControlWZYX = { 1.0f,-1.0f, 1.0f,-1.0f };
		static const __m128 ControlZWXY = { 1.0f, 1.0f,-1.0f,-1.0f };
		static const __m128 ControlYXWZ = { -1.0f, 1.0f, 1.0f,-1.0f };
		__m128 Q2X = Q2, Q2Y = Q2, Q2Z = Q2, vResult = Q2;
		vResult = _mm_shuffle_ps(vResult, vResult, _mm_shuffle(3, 3, 3, 3));
		Q2X = _mm_shuffle_ps(Q2X, Q2X, _mm_shuffle(0, 0, 0, 0));
		Q2Y = _mm_shuffle_ps(Q2Y, Q2Y, _mm_shuffle(1, 1, 1, 1));
		Q2Z = _mm_shuffle_ps(Q2Z, Q2Z, _mm_shuffle(2, 2, 2, 2));
		vResult = _mm_mul_ps(vResult, Q1);
		__m128 Q1Shuffle = Q1;
		Q1Shuffle = _mm_shuffle_ps(Q1Shuffle, Q1Shuffle, _mm_shuffle(0, 1, 2, 3));
		Q2X = _mm_mul_ps(Q2X, Q1Shuffle);
		Q1Shuffle = _mm_shuffle_ps(Q1Shuffle, Q1Shuffle, _mm_shuffle(2, 3, 0, 1));
		Q2X = _mm_mul_ps(Q2X, ControlWZYX);
		Q2Y = _mm_mul_ps(Q2Y, Q1Shuffle);
		Q1Shuffle = _mm_shuffle_ps(Q1Shuffle, Q1Shuffle, _mm_shuffle(0, 1, 2, 3));
		Q2Y = _mm_mul_ps(Q2Y, ControlZWXY);
		Q2Z = _mm_mul_ps(Q2Z, Q1Shuffle);
		vResult = _mm_add_ps(vResult, Q2X);
		Q2Z = _mm_mul_ps(Q2Z, ControlYXWZ);
		Q2Y = _mm_add_ps(Q2Y, Q2Z);
		vResult = _mm_add_ps(vResult, Q2Y);
		return vResult;
	}

	inline static __m128 VECTORCALL MulVec3(__m128 vec, __m128 quat)
	{
		__m128 temp = SSEVector3Cross(quat, vec);
		__m128 temp1 = _mm_mul_ps(vec, SSESplatZ(quat)) ;
		temp = _mm_add_ps(temp, temp1);
		temp1 = _mm_mul_ps(SSEVector3Cross(quat, temp), _mm_set1_ps(2.0f));
		return _mm_add_ps(vec, temp1);
	}

	__forceinline static __m128 VECTORCALL Dot(const Quaternion V1, const Quaternion V2) noexcept
	{
		__m128 vTemp2 = V2.vec;
		__m128 vTemp = _mm_mul_ps(V1.vec, vTemp2);
		vTemp2 = _mm_shuffle_ps(vTemp2, vTemp, _mm_shuffle(1, 0, 0, 0)); // Copy X to the Z position and Y to the W position
		vTemp2 = _mm_add_ps(vTemp2, vTemp);          // Add Z = X+Z; W = Y+W;
		vTemp = _mm_shuffle_ps(vTemp, vTemp2, _mm_shuffle(0, 3, 0, 0));  // Copy W to the Z position
		vTemp = _mm_add_ps(vTemp, vTemp2);           // Add Z and W together
		return _mm_shuffle_ps(vTemp, vTemp, _mm_shuffle(2, 2, 2, 2));    // Splat Z and return
	}

	inline Quaternion VECTORCALL Slerp(const Quaternion Q0, const Quaternion Q1, float t) noexcept
	{
		const __m128 T = _mm_set_ps1(t);
		// Result = Q0 * sin((1.0 - t) * Omega) / sin(Omega) + Q1 * sin(t * Omega) / sin(Omega)
		static const __m128 OneMINusEpsilon = _mm_setr_ps(1.0f - 0.00001f, 1.0f - 0.00001f, 1.0f - 0.00001f, 1.0f - 0.00001f);
		static const __m128 SignMask2 = _mm_castsi128_ps(_mm_setr_epi32(0x80000000, 0x00000000, 0x00000000, 0x00000000));

		__m128 CosOmega = _mm_dp_ps(Q0.vec, Q1.vec, 0xff);

		const __m128 Zero = _mm_setzero_ps();
		__m128 Control = _mm_cmplt_ps(CosOmega, Zero);
		__m128 Sign = SSESelect(g_XMOne, g_XMNegativeOne, Control);

		CosOmega = _mm_mul_ps(CosOmega, Sign);

		Control = _mm_cmplt_ps(CosOmega, OneMINusEpsilon);

		__m128 SinOmega = _mm_mul_ps(CosOmega, CosOmega);
		SinOmega = _mm_sub_ps(g_XMOne, SinOmega);
		SinOmega = _mm_sqrt_ps(SinOmega);

		__m128 Omega = _mm_atan2_ps(SinOmega, CosOmega);

		__m128 V01 = _mm_permute_ps(T, _mm_shuffle(2, 3, 0, 1));
		V01 = _mm_and_ps(V01, g_XMMaskXY);
		V01 = _mm_xor_ps(V01, SignMask2);
		V01 = _mm_add_ps(g_XMIdentityR0, V01);

		__m128 S0 = _mm_mul_ps(V01, Omega);
		S0 = _mm_sin_ps(S0);
		S0 = _mm_div_ps(S0, SinOmega);

		S0 = SSESelect(V01, S0, Control);

		__m128 S1 = SSESplatY(S0);
		S0 = SSESplatX(S0);

		S1 = _mm_mul_ps(S1, Sign);
		__m128 Result = _mm_mul_ps(Q0.vec, S0);
		S1 = _mm_mul_ps(S1, Q1.vec); // _mm_fmadd_ps(S1, Q1.vec, Result) ?
		Quaternion res;
		res.vec = _mm_add_ps(Result, S1);
		return res;
	}

	__forceinline Quaternion static FromEuler(float x, float y, float z) noexcept
	{
		x *= 0.5f; y *= 0.5f; z *= 0.5f;
		float c[4], s[4];
		__m128 cv;
		__m128 sv = _mm_sincos_ps(&cv, _mm_setr_ps(x, y, z, 1.0f));
		_mm_store_ps(c, cv);
		_mm_store_ps(s, sv);

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

	__forceinline static __m128 VECTORCALL Conjugate(const __m128 vec)
	{
		static const __m128 NegativeOne3 = _mm_setr_ps(-1.0f,-1.0f,-1.0f, 1.0f);
		return _mm_mul_ps(vec, NegativeOne3);
	}

	Vector3f GetForward() const {
		Vector3f res;
		SSEStoreVector3(&res.x,MulVec3(_mm_setr_ps( 0.0f, 0.0f, -1.0f, 0.0f), Conjugate(vec)));
		return res; 
	}

	Vector3f GetRight() const {
		Vector3f res;
		SSEStoreVector3(&res.x, MulVec3(_mm_setr_ps( 1.0f, 0.0f, 0.0f, 0.0f), Conjugate(vec)));
		return res; 
	}

	Vector3f GetLeft() const {
		Vector3f res;
		SSEStoreVector3(&res.x, MulVec3(_mm_setr_ps(-1.0f, 0.0f, 0.0f, 0.0f), Conjugate(vec))); 
		return res; 
	}

	Vector3f GetUp() const {
		Vector3f res;
		SSEStoreVector3(&res.x, MulVec3(_mm_setr_ps( 0.0f, 1.0f, 0.0f, 0.0f), Conjugate(vec)));
		return res; 
	}

	__forceinline Quaternion operator *  (const Quaternion& b) { Quaternion q; q.vec = Mul(this->vec, b.vec); return q; }
	__forceinline Quaternion operator *= (const Quaternion& b) { this->vec = Mul(this->vec, b.vec); return *this; }
};

__forceinline Quaternion MakeQuat(float scale = 0.0f)                     
{
	Quaternion v;
	v.vec = _mm_set1_ps(scale);
	return v;
}

__forceinline Quaternion MakeQuat(float _x, float _y, float _z, float _w) 
{
	Quaternion v;
	v.vec = _mm_setr_ps(_x, _y, _z, _w);   
	return v; 
}

#else
// todo
#endif // AX_SUPPORT_SSE

AX_END_NAMESPACE 