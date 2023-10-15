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

	inline static Vector3f VECTORCALL MulVec3(Vector3f vec, Quaternion quat)
	{
		Vector3f res;
		SSEStoreVector3(&res.x, MulVec3(_mm_setr_ps(vec.x, vec.y, vec.z, 1.0f), quat.vec));
		return res;
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
		__m128 Sign = SSESelect(g_XOne, g_XNegativeOne, Control);

		CosOmega = _mm_mul_ps(CosOmega, Sign);

		Control = _mm_cmplt_ps(CosOmega, OneMINusEpsilon);

		__m128 SinOmega = _mm_mul_ps(CosOmega, CosOmega);
		SinOmega = _mm_sub_ps(g_XOne, SinOmega);
		SinOmega = _mm_sqrt_ps(SinOmega);

		__m128 Omega = _mm_atan2_ps(SinOmega, CosOmega);

		__m128 V01 = _mm_permute_ps(T, _mm_shuffle(2, 3, 0, 1));
		V01 = _mm_and_ps(V01, g_XMaskXY);
		V01 = _mm_xor_ps(V01, SignMask2);
		V01 = _mm_add_ps(g_XIdentityR0, V01);

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

	__forceinline Quaternion Inversed() const
	{
		return Inverse(*this);
	}

	__forceinline static __m128 VECTORCALL Conjugate(const __m128 vec)
	{
		static const __m128 NegativeOne3 = _mm_setr_ps(-1.0f,-1.0f,-1.0f, 1.0f);
		return _mm_mul_ps(vec, NegativeOne3);
	}

	__forceinline static Quaternion VECTORCALL Conjugate(const Quaternion vec)
	{
		static const Quaternion NegativeOne3 = {-1.0f, -1.0f, -1.0f, 1.0f};
		return vec * NegativeOne3;
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

	__forceinline Quaternion operator *  (float b)             const { Quaternion q; q.vec = _mm_mul_ps(this->vec, _mm_set1_ps(b)); return q; }
	__forceinline Quaternion operator *= (float b)                   { this->vec = Mul(this->vec, _mm_set1_ps(b)); return *this; }
	__forceinline Quaternion operator *  (const Quaternion& b) const { Quaternion q; q.vec = Mul(this->vec, b.vec); return q; }
	__forceinline Quaternion operator *= (const Quaternion& b)       { this->vec = Mul(this->vec, b.vec); return *this; }
	__forceinline Quaternion operator +  (const Quaternion& b) const { Quaternion q; q.vec = _mm_add_ps(vec, b.vec); return q; }
	__forceinline Quaternion operator += (const Quaternion& b)       { vec = _mm_add_ps(vec, b.vec); return *this; }
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
		Vector3f matrix[3] {
			Vector3f::Cross(up, direction), up, direction
		};
		QuaternionFromMatrix<3>(&result.x, &matrix[0].x);
		return result;
	}

	static inline Quaternion FromToRotation(const Quaternion& start, const Quaternion& end) { return start.Inversed() * end; }

	inline static Quaternion VECTORCALL Mul(const Quaternion q1, const Quaternion q2) noexcept
	{
		const float x    = q1.x, y    = q1.y, z    = q1.z, w   = q1.w;
		const float num4 = q2.x, num3 = q2.y, num2 = q2.z, num = q2.w;
		
		const float num12 = (y * num2) - (z * num3);
		const float num11 = (z * num4) - (x * num2);
		const float num10 = (x * num3) - (y * num4);
		const float num9  = ((x * num4) + (y * num3)) + (z * num2);

		return {
			((x * num) + (num4 * w)) + num12,
			((y * num) + (num3 * w)) + num11,
			((z * num) + (num2 * w)) + num10,
			(w * num) - num9 
		};
	}

	inline static Vector3f MulVec3(Vector3f vec, Quaternion quat)
	{
		Vector3f qxyz  = MakeVec3(quat.x, quat.y, quat.z);
		Vector3f temp  = Vector3f::Cross(qxyz, vec);
		Vector3f temp1 = vec * quat.z;
		temp += temp1;
		temp1 = Vector3f::Cross(qxyz, temp) * 2.0f;
		return vec + temp1;
	}

	inline Quaternion Slerp(const Quaternion Q0, const Quaternion Q1, float t) noexcept
	{
		const Vector4f T = MakeVec4(t);
		// Result = Q0 * sin((1.0 - t) * Omega) / sin(Omega) + Q1 * sin(t * Omega) / sin(Omega)
		float CosOmega = Vector4f::Dot(Q0.vec, Q1.vec);
		float sign = CosOmega < 0.0f ? -1.0f : 1.0f;
		CosOmega = CosOmega * sign;
		
		bool Control = CosOmega < (1.0f - 0.0001f);
		float SinOmega = CosOmega * CosOmega;
		SinOmega = Sqrt(1.0f - SinOmega);
		
		float Omega  = ATan2(SinOmega, CosOmega);
		Vector4f V01 = MakeVec4(T.z, T.w, T.x, T.y); 
		uint uf = BitCast<uint>(V01.x); // V01 = _mm_and_ps(V01, g_XMaskXY);
		uf ^= 0x80000000;               // V01 = _mm_xor_ps(V01, SignMask2);
		V01.x = BitCast<float>(uf);     // V01 = _mm_add_ps(g_XIdentityR0, V01);
		V01.x += 1.0f;

		Vector4f S0 = V01 * Omega;
		S0.x = Sin(S0.x) / SinOmega; // S0 = _mm_sin_ps(S0);
		S0.y = Sin(S0.y) / SinOmega; // S0 = _mm_div_ps(S0, SinOmega);
		S0 = Control ? S0 : V01; // SSESelect(V01, S0, Control)

		Vector4f S1 = MakeVec4(S0.y);
		S0 = MakeVec4(S0.x);
		S1 = S1 * sign;
		Vector4f Result = Q0.vec * S0;
		S1 = S1 * Q1.vec; // _mm_fmadd_ps(S1, Q1.vec, Result) ?
		Quaternion res;
		res.vec = Result + S1;
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
		const float lengthSq = q.vec.LengthSquared();
		if (lengthSq == 1.0f)
		{
			q = Conjugate(q);
			return q;
		}
		else if (lengthSq >= 0.001f)
		{
			q = Conjugate(q) * (1.0f / lengthSq);
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

	__forceinline static Quaternion VECTORCALL Conjugate(Quaternion v)
	{
		const Quaternion NegativeOne3 = { -1.0f, -1.0f, -1.0f, 1.0f };
		return v * NegativeOne3;
	}

	Vector3f GetForward() const {
		return MulVec3(MakeVec3(0.0f, 0.0f, -1.0f), Conjugate(*this)); 
	}

	Vector3f GetRight() const {
		return MulVec3(MakeVec3( 1.0f, 0.0f, 0.0f), Conjugate(*this)); 
	}

	Vector3f GetLeft() const {
		return MulVec3(MakeVec3(-1.0f, 0.0f, 0.0f), Conjugate(*this)); 
	}

	Vector3f GetUp() const {
		return MulVec3(MakeVec3( 0.0f, 1.0f, 0.0f), Conjugate(*this)); 
	}
	
    __forceinline Quaternion operator *  (float f)             const  { return { x * f, y * f, z * f, w * f }; }
	__forceinline Quaternion operator *= (float f)                    { x *= f; y *= f; z *= f; w *= f; return *this; }
	__forceinline Quaternion operator *  (const Quaternion& b) const  { Quaternion q; q.vec = this->vec * b.vec; return q; }
	__forceinline Quaternion operator *= (const Quaternion& b)        { this->vec = this->vec * b.vec; return *this; }
	__forceinline Quaternion operator +  (const Quaternion& b) const  { Vector4f v = vec; v += b.vec; return BitCast<Quaternion>(v); }
	__forceinline Quaternion operator += (const Quaternion& b)        { this->vec += b.vec; return *this; }
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