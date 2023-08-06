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

	inline static Vector4f VECTORCALL Mul(Vector4f Q1, Vector4f Q2)
	{
		return BitCast<Vector4f>(Mul(Q1.vec, Q2.vec));
	}

	inline static __m128 VECTORCALL Mul(const __m128 Q1, const __m128 Q2) 
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

	inline static Quaternion VECTORCALL MulVec3(Vector3f vec, Quaternion quat)
	{
		return MulVec3(_mm_setr_ps(vec.x, vec.y, vec.z, 1.0f), quat.vec);
	}

	inline static __m128 VECTORCALL MulVec3(__m128 vec, __m128 quat)
	{
		__m128 temp = SSEVector3Cross(quat, vec);
		__m128 temp1 = _mm_mul_ps(vec, SSESplatZ(quat)) ;
		temp = _mm_add_ps(temp, temp1);
		temp1 = _mm_mul_ps(SSEVector3Cross(quat, temp), _mm_set1_ps(2.0f));
		return _mm_add_ps(vec, temp1);
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

	__forceinline static Quaternion Inverse(Quaternion q)
	{
		const float lengthSq = _mm_cvtss_f32(_mm_dp_ps(q.vec, q.vec, 0xff));
		if (lengthSq == 1.0f)
		{
			q.vec = Conjugate(q.vec);
			return q;
		}
		else if (lengthSq >= 0.001f)
		{
			q.vec = _mm_mul_ps(Conjugate(q.vec), _mm_set1_ps(1.0f / lengthSq));
			return q;
		}
		else
		{
			return Identity();
		}
	}

	__forceinline static Quaternion Inversed() const
	{
		return Inverse(*this);
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
	__forceinline Quaternion operator +  (const Quaternion& b) { Quaternion q; q.vec = _mm_add_ps(vec, b.vec); return q; }
	__forceinline Quaternion operator += (const Quaternion& b) { vec = _mm_add_ps(vec, b.vec); return *this; }
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

struct alignas(16) Quaternion
{
	union
	{
		struct { float x, y, z, w; };
		float arr[4];
		Vector4f vec;
	};

	const float  operator [] (int index) const { return arr[index]; }
	float& operator [] (int index) { return arr[index]; }

	__forceinline static Quaternion Identity() 
	{
		Quaternion q;
		q.x = q.y = q.z = 0.0f; q.w = 1.0f;
		return q;
	}

	__forceinline static Quaternion FromAngleAxis(float angle, const Vector3f& axis)
	{
		const float hlf = angle * 0.5f;
		const float s   = Sin(hlf);
		const float c   = Cos(hlf);
		return { axis.x * s, axis.y * s, axis.z * s, c };
	}

	static inline Quaternion FromLookRotation(Vector3f direction, const Vector3f& up)
	{
		Quaternion result;
		const Vector3f forward = direction.Normalized();

		Vector3f v = Vector3f::Cross(forward, up);
		if (v.LengthSquared() >= 0.001f)
		{
			v.NormalizeSelf();
			const Vector3f up    = Vector3f::Cross(v, forward);
			const Vector3f right = Vector3f::Cross(up, forward);
			result.FromAxes(right, up, forward);
		}
		else
		{
			result = Quaternion::FromToRotation(Vector3::Forward, forward);
		}

		return result;
	}

	static inline Quaternion FromToRotation(const Quaternion& start, const Quaternion& end) { return start.Inversed() * end; }

	inline static Quaternion VECTORCALL Mul(const Quaternion q1, const Quaternion q2) noexcept
	{
		static const Vector4f ControlWZYX = { 1.0f,-1.0f, 1.0f,-1.0f };
		static const Vector4f ControlZWXY = { 1.0f, 1.0f,-1.0f,-1.0f };
		static const Vector4f ControlYXWZ = { -1.0f, 1.0f, 1.0f,-1.0f };
		Quaternion Q1 = q1.vec, Q2 = q2.vec;
		Vector4f Q2X = Q2, Q2Y = Q2, Q2Z = Q2, vResult = Q2;
		vResult = VecSwizzle1(vResult, 3);
		Q2X = VecSwizzle1(Q2X, 0);
		Q2Y = VecSwizzle1(Q2Y, 1);
		Q2Z = VecSwizzle1(Q2Z, 2);
		vResult *= Q1;
		Vector4f Q1Shuffle = Q1;
		Q1Shuffle = VecShuffle(Q1Shuffle, Q1Shuffle, 0, 1, 2, 3);
		Q2X *= Q1Shuffle;
		Q1Shuffle = VecShuffle(Q1Shuffle, Q1Shuffle, 2, 3, 0, 1);
		Q2X *= ControlWZYX;
		Q2Y *= Q1Shuffle;
		Q1Shuffle = VecShuffle(Q1Shuffle, Q1Shuffle, 0, 1, 2, 3);
		Q2Y *= ControlZWXY);
		Q2Z *= Q1Shuffle);
		vResult += Q2X;
		Q2Z += ControlYXWZ;
		Q2Y += Q2Z;
		vResult = vResult + Q2Y;
		return vResult;
	}

	inline static Vector4f VECTORCALL MulVec3(Vector3f vec, Vector4f quat)
	{
		Vector3f temp = Vector3f::Cross(quat, vec);
		Vector3f temp1 = vec * quat.z;
		temp += temp1;
		temp1 = Vector3f::Cross(quat, temp) * MakeVec3(2.0f);
		return vec + temp1;
	}

	inline Quaternion VECTORCALL Slerp(const Quaternion Q0, const Quaternion Q1, float t) noexcept
	{
		const Vector4f T = _mm_set_ps1(t);
		// Result = Q0 * sin((1.0 - t) * Omega) / sin(Omega) + Q1 * sin(t * Omega) / sin(Omega)
		static const Vector4f OneMINusEpsilon = { 1.0f - 0.00001f, 1.0f - 0.00001f, 1.0f - 0.00001f, 1.0f - 0.00001f };
		static const Vector4f SignMask2 = _mm_castsi128_ps(_mm_setr_epi32(0x80000000, 0x00000000, 0x00000000, 0x00000000));

		float CosOmega = Vector4f::Dot(Q0.vec, Q1.vec);

		float sign = CosOmega < 0.0f ? -1.0f : 1.0f;

		CosOmega = CosOmega * Sign;

		bool Control = CosOmega < OneMINusEpsilon;

		float SinOmega = CosOmega * CosOmega;
		SinOmega = Sqrt(1.0f - SinOmega);

		float Omega = Atan2(SinOmega, CosOmega);

		Vector4f V01 = MakeVec4(T.z, T.w, T.x, T.y); 
		V01 = _mm_and_ps(V01, g_XMMaskXY);
		V01 = _mm_xor_ps(V01, SignMask2);
		V01 = _mm_add_ps(g_XMIdentityR0, V01);

		Vector4f S0 = _mm_mul_ps(V01, Omega);
		S0 = _mm_sin_ps(S0);
		S0 = _mm_div_ps(S0, SinOmega);

		S0 = SSESelect(V01, S0, Control);

		Vector4f S1 = SSESplatY(S0);
		S0 = SSESplatX(S0);

		S1 = _mm_mul_ps(S1, Sign);
		Vector4f Result = _mm_mul_ps(Q0.vec, S0);
		S1 = _mm_mul_ps(S1, Q1.vec); // _mm_fmadd_ps(S1, Q1.vec, Result) ?
		Quaternion res;
		res.vec = _mm_add_ps(Result, S1);
		return res;
	}

	__forceinline static Quaternion FromEuler(float x, float y, float z) noexcept
	{
		x *= 0.5f; y *= 0.5f; z *= 0.5f;
		Vector3f s { Sin(x), Sin(y), Sin(z) };
		Vector3f c { Cos(x), Cos(y), Cos(z) };

		Quaternion q;
		q.w = c[0] * c[1] * c[2] + s[0] * s[1] * s[2];
		q.x = s[0] * c[1] * c[2] - c[0] * s[1] * s[2];
		q.y = c[0] * s[1] * c[2] + s[0] * c[1] * s[2];
		q.z = c[0] * c[1] * s[2] - s[0] * s[1] * c[2];
		return q;
	}

	inline static Vector3f ToEulerAngles(const Quaternion& q) noexcept {
		Vector3f eulerAngles; // using std is recommended
		eulerAngles.x = ATan2(2.0f * (q.y * q.z + q.w * q.x), q.w * q.w - q.x * q.x - q.y * q.y + q.z * q.z);
		eulerAngles.y = ASin(Clamp(-2.0f * (q.x * q.z - q.w * q.y), -1.0f, 1.0f));
		eulerAngles.z = ATan2(2.0f * (q.x * q.y + q.w * q.z), q.w * q.w + q.x * q.x - q.y * q.y - q.z * q.z);
		return eulerAngles;
	}
	
	__forceinline static Quaternion VECTORCALL FromEuler(Vector3f euler) noexcept
	{
		return FromEuler(euler.x, euler.y, euler.z);
	}

	__forceinline static Quaternion Inverse(Quaternion q)
	{
		const float lengthSq = q.vec.LengthSq();
		if (lengthSq == 1.0f)
		{
			q.vec = Conjugate(q.vec);
			return q;
		}
		else if (lengthSq >= 0.001f)
		{
			q.vec = Conjugate(q.vec) * (1.0f / lengthSq);
			return q;
		}
		else
		{
			return Identity();
		}
	}

	__forceinline static Quaternion Inversed() const
	{
		return Inverse(*this);
	}

	__forceinline static Vector4f VECTORCALL Conjugate(const Vector4f vec)
	{
		static const Vector4f NegativeOne3 = { -1.0f, -1.0f, -1.0f, 1.0f };
		return vec * NegativeOne3;
	}

	Vector3f GetForward() const {
		return MulVec3(MakeVec3(0.0f, 0.0f, -1.0f), Conjugate(vec)); 
	}

	Vector3f GetRight() const {
		return MulVec3(MakeVec3( 1.0f, 0.0f, 0.0f), Conjugate(vec))); 
	}

	Vector3f GetLeft() const {
		return MulVec3(MakeVec3(-1.0f, 0.0f, 0.0f), Conjugate(vec)); 
	}

	Vector3f GetUp() const {
		return MulVec3(MakeVec3( 0.0f, 1.0f, 0.0f), Conjugate(vec))); 
	}
	
    __forceinline Quaternion operator *  (float f) { return {q.x * f, q.y * f, q.z * f, q.w * f } }
	__forceinline Quaternion operator *= (float f) { x *= f; y *= f; z *= f; w *= f; return *this; }
	__forceinline Quaternion operator *  (const Quaternion& b) { Quaternion q; q.vec = Mul(this->vec, b.vec); return q; }
	__forceinline Quaternion operator *= (const Quaternion& b) { this->vec = Mul(this->vec, b.vec); return *this; }
	__forceinline Quaternion operator +  (const Quaternion& b) { Vector4f v = vec; v += b.vec; return BitCast<Quaternion>(v); }
	__forceinline Quaternion operator += (const Quaternion& b) { this->vec += b.vec; return *this; }
};

__forceinline Quaternion MakeQuat(float scale = 0.0f)                     
{
	Quaternion v;
	v.x = v.y = v.z = v.w = scale;
	return v;
}

__forceinline Quaternion MakeQuat(float _x, float _y, float _z, float _w) 
{
	Quaternion v{_x, _y, _z, _w};
	return v; 
}

#endif // AX_SUPPORT_SSE

AX_END_NAMESPACE 