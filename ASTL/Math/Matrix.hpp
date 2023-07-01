#pragma once
#include "Quaternion.hpp"

struct Matrix3
{
	union
	{
		float m[3][3] = {0};
		Vector3f vec[3];
		struct { Vector3f x, y, z; };
	};

	const Vector3f& GetForward() const { return vec[2]; }
	const Vector3f& GetUp()      const { return vec[1]; }
	const Vector3f& GetRight()   const { return vec[0]; }

	Matrix3()
	{
		m[0][0] = 1.0f;
		m[1][1] = 1.0f;
		m[2][2] = 1.0f;
	}

	Matrix3(float s)
	{
		for (int i = 0; i < 9; ++i) m[0][i] = s;
	}

	FINLINE static Matrix3 Identity()
	{
		Matrix3 mat;
		mat.m[0][0] = 1.0f;
		mat.m[1][1] = 1.0f;
		mat.m[2][2] = 1.0f;
		return mat;
	}

	inline static Matrix3 LookAt(Vector3f direction, Vector3f up)
	{
		Matrix3 result;
		result.vec[2] = direction;
		Vector3f const& Right = Vector3f::Cross(up, result.vec[2]);
		result.vec[0] = Right * RSqrt(Max(0.00001f, Vector3f::Dot(Right, Right)));
		result.vec[1] = Vector3f::Cross(result.vec[2], result.vec[0]);
		return result;
	}

	inline static Matrix3 Multiply(const Matrix3& a, const Matrix3& b)
	{
		Matrix3 result;
		float3 vx = a.x.xxx() * b.x;
		float3 vy = a.x.yyy() * b.y;
		float3 vz = a.x.zzz() * b.z;

		vx += vz; vx += vy;
		result.x = vx;

		vx = a.y.xxx() * b.x;
		vy = a.y.yyy() * b.y;
		vz = a.y.zzz() * b.z;
		vx += vz; vx += vy;
		result.y = vx;

		vx = a.z.xxx() * b.x;
		vy = a.z.yyy() * b.y;
		vz = a.z.zzz() * b.z;
		vx += vz; vx += vy;
		result.z = vx;
		return result;
	}

	inline static float3 Multiply(const Matrix3& m, const float3& v) {
		return m.x * v.xxx() + m.y * v.yyy() + m.z * v.zzz();
	}

	inline static Matrix3 FromQuaternion(Quaternion q)
	{
		// ai generated, todo fix
		Matrix3 rot;
		float x2 = q.x * q.x, y2 = q.y * q.y;
		float z2 = q.z * q.z, xy = q.x * q.y;
		float xz = q.x * q.z, yz = q.y * q.z;
		float wx = q.w * q.x, wy = q.w * q.y;
		float wz = q.w * q.z;
		rot.m[0][0] = 1.0f - 2.0f * (y2 + z2);
		rot.m[0][1] = 2.0f * (xy - wz);
		rot.m[0][2] = 2.0f * (xz + wy);
		rot.m[1][0] = 2.0f * (xy + wz);
		rot.m[1][1] = 1.0f - 2.0f * (x2 + z2);
		rot.m[1][2] = 2.0f * (yz - wx);
		rot.m[2][0] = 2.0f * (xz - wy);
		rot.m[2][1] = 2.0f * (yz + wx);
		rot.m[2][2] = 1.0f - 2.0f * (x2 + y2);
		return rot;
	}

	Quaternion ToQuaternion() const
	{
		float fourXSquaredMinus1 = m[0][0] - m[1][1] - m[2][2];
		float fourYSquaredMinus1 = m[1][1] - m[0][0] - m[2][2];
		float fourZSquaredMinus1 = m[2][2] - m[0][0] - m[1][1];
		float fourWSquaredMinus1 = m[0][0] + m[1][1] + m[2][2];

		int biggestIndex = 0;
		float fourBiggestSquaredMinus1 = fourWSquaredMinus1;
		if(fourXSquaredMinus1 > fourBiggestSquaredMinus1)
		{
			fourBiggestSquaredMinus1 = fourXSquaredMinus1;
			biggestIndex = 1;
		}
		if(fourYSquaredMinus1 > fourBiggestSquaredMinus1)
		{
			fourBiggestSquaredMinus1 = fourYSquaredMinus1;
			biggestIndex = 2;
		}
		if(fourZSquaredMinus1 > fourBiggestSquaredMinus1)
		{
			fourBiggestSquaredMinus1 = fourZSquaredMinus1;
			biggestIndex = 3;
		}

		float biggestVal = Sqrt(fourBiggestSquaredMinus1 + 1.0f) * 0.5f;
		float mult = 0.25f / biggestVal;
		
		switch(biggestIndex)
		{
			case 0:
				return Quaternion(biggestVal, (m[1][2] - m[2][1]) * mult, (m[2][0] - m[0][2]) * mult, (m[0][1] - m[1][0]) * mult);
			case 1:
				return Quaternion((m[1][2] - m[2][1]) * mult, biggestVal, (m[0][1] + m[1][0]) * mult, (m[2][0] + m[0][2]) * mult);
			case 2:
				return Quaternion((m[2][0] - m[0][2]) * mult, (m[0][1] + m[1][0]) * mult, biggestVal, (m[1][2] + m[2][1]) * mult);
			case 3:
				return Quaternion((m[0][1] - m[1][0]) * mult, (m[2][0] + m[0][2]) * mult, (m[1][2] + m[2][1]) * mult, biggestVal);
		}
	}
};

// todo make non simd version
#ifdef AX_SUPPORT_SSE2

AX_ALIGNED(16) struct Matrix4
{
	union
	{
		struct { __m128 r[4]; };
		struct { float m[4][4]; };
	};

	const __m128& operator [] (int index) const { return r[index]; }
	__m128& operator [] (int index) { return r[index]; }
	
	Vector4f VECTORCALL  operator *  (const Vector3f v)  noexcept { return Vector3Transform(v, *this); };
	Vector4f VECTORCALL  operator *  (const Vector4f& v)   noexcept { return Vector4Transform(v, *this); };

	Matrix4 VECTORCALL  operator *  (const Matrix4& M)  noexcept { return Matrix4::Multiply(*this, M); };
	Matrix4& VECTORCALL operator *= (const Matrix4& M)  noexcept { *this = Matrix4::Multiply(*this, M); return *this; };

	Matrix4() {}
	explicit Matrix4(EForceInit) 
	{
		r[0] = g_XMIdentityR0;
		r[1] = g_XMIdentityR1;
		r[2] = g_XMIdentityR2;
		r[3] = g_XMIdentityR3;
	}

	explicit Matrix4(float s)
	{
		r[0] = r[1] = r[2] = r[3] = _mm_set_ps1(s);
	}

	FINLINE static Matrix4 Identity()
	{
		Matrix4 M;
		M.r[0] = g_XMIdentityR0;
		M.r[1] = g_XMIdentityR1;
		M.r[2] = g_XMIdentityR2;
		M.r[3] = g_XMIdentityR3;
		return M;
	}

	FINLINE static Matrix4 FromPosition(const float x, const float y, const float z)
	{
		Matrix4 M;
		M.r[0] = g_XMIdentityR0;
		M.r[1] = g_XMIdentityR1;
		M.r[2] = g_XMIdentityR2;
		M.r[3] = _mm_set_ps(1.0f, z, y, x);
		return M;
	}

	FINLINE static Matrix4 FromPosition(const Vector3f& vec3)
	{
		return FromPosition(vec3.x, vec3.y, vec3.z);
	}

	FINLINE static Matrix4 CreateScale(const float ScaleX, const float ScaleY, const float ScaleZ)
	{
		Matrix4 M;
		M.r[0] = _mm_set_ps(0.0f, 0.0f, 0.0f, ScaleX);
		M.r[1] = _mm_set_ps(0.0f, 0.0f, ScaleY, 0.0f);
		M.r[2] = _mm_set_ps(0.0f, ScaleZ, 0.0f, 0.0f);
		M.r[3] = g_XMIdentityR3;
		return M;
	}

	FINLINE static Matrix4 CreateScale(const Vector3f& vec3)
	{
		return CreateScale(vec3.x, vec3.y, vec3.z);
	}

	// please assign normalized vectors, returns view matrix
	FINLINE static Matrix4 VECTORCALL LookAtRH(Vector3f eye, Vector3f center, const Vector3f& up)
	{
		__m128 NegEyePosition;
		__m128 D0, D1, D2;
		__m128 R0, R1;

		__m128 EyePosition  = _mm_loadu_ps(&eye.x);
		__m128 EyeDirection = _mm_sub_ps(_mm_setzero_ps(), _mm_loadu_ps(&center.x));
		__m128 UpDirection  = _mm_loadu_ps(&up.x);

		R0 = SSEVectorNormalize(SSEVector3Cross(UpDirection, EyeDirection));
		R1 = SSEVectorNormalize(SSEVector3Cross(EyeDirection, R0));

		NegEyePosition = _mm_sub_ps(_mm_setzero_ps(), EyePosition);

		D0 = SSEVector3Dot(R0, NegEyePosition);
		D1 = SSEVector3Dot(R1, NegEyePosition);
		D2 = SSEVector3Dot(EyeDirection, NegEyePosition);
		Matrix4 M;
		M.r[0] = SSESelect(D0, R0, g_XMSelect1110);
		M.r[1] = SSESelect(D1, R1, g_XMSelect1110);
		M.r[2] = SSESelect(D2, EyeDirection, g_XMSelect1110);
		M.r[3] = g_XMIdentityR3;
		return Matrix4::Transpose(M);
	}

	FINLINE static Matrix4 PerspectiveFovRH(float fov, float width, float height, float zNear, float zFar)
	{
		const float rad = fov;
		const float h = Cos(0.5f * rad) / Sin(0.5f * rad);
		const float w = h * height / width; /// max(width , Height) / min(width , Height)?
		Matrix4 M(ForceInit);
		M.m[0][0] = w;
		M.m[1][1] = h;
		M.m[2][2] = -(zFar + zNear) / (zFar - zNear);
		M.m[2][3] = -1.0f;
		M.m[3][2] = -(2.0f * zFar * zNear) / (zFar - zNear);
		M.m[3][3] = 0.0f;
		return M;
	}

	// https://lxjk.github.io/2017/09/03/Fast-4x4-Matrix-Inverse-with-SSE-SIMD-Explained.html
	#define MakeShuffleMask(x,y,z,w)           (x | (y<<2) | (z<<4) | (w<<6))
	// vec(0, 1, 2, 3) -> (vec[x], vec[y], vec[z], vec[w])
	#define VecSwizzleMask(vec, mask)        _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(vec), mask))
	#define VecSwizzle(vec, x, y, z, w)      VecSwizzleMask(vec, MakeShuffleMask(x,y,z,w))
	#define VecSwizzle1(vec, x)              VecSwizzleMask(vec, MakeShuffleMask(x,x,x,x))
	// special swizzle                   
	#define VecSwizzle_0022(vec)             _mm_moveldup_ps(vec)
	#define VecSwizzle_1133(vec)             _mm_movehdup_ps(vec)
	
	// return (vec1[x], vec1[y], vec2[z], vec2[w])
	#define VecShuffle(vec1, vec2, x,y,z,w)  _mm_shuffle_ps(vec1, vec2, MakeShuffleMask(x,y,z,w))
	// special shuffle                  
	#define VecShuffle_0101(vec1, vec2)      _mm_movelh_ps(vec1, vec2)
	#define VecShuffle_2323(vec1, vec2)      _mm_movehl_ps(vec2, vec1)

	// for row major matrix
    // we use __m128 to represent 2x2 matrix as A = | A0  A1 |
    //                                              | A2  A3 |
    // 2x2 row major Matrix multiply A*B
	static FINLINE __m128 VECTORCALL Mat2Mul(__m128 vec1, __m128 vec2)
	{
		return _mm_add_ps(_mm_mul_ps(vec1, VecSwizzle(vec2, 0, 3, 0, 3)), 
		                  _mm_mul_ps(VecSwizzle(vec1, 1, 0, 3, 2), VecSwizzle(vec2, 2, 1, 2, 1)));
	}
	// 2x2 row major Matrix adjugate multiply (A#)*B
	static FINLINE __m128 VECTORCALL Mat2AdjMul(__m128 vec1, __m128 vec2)
	{
		return _mm_sub_ps(_mm_mul_ps(VecSwizzle(vec1, 3, 3, 0, 0), vec2),
		                  _mm_mul_ps(VecSwizzle(vec1, 1, 1, 2, 2), VecSwizzle(vec2, 2, 3, 0, 1)));

	}
	// 2x2 row major Matrix multiply adjugate A*(B#)
	static FINLINE __m128 VECTORCALL Mat2MulAdj(__m128 vec1, __m128 vec2)
	{
		return _mm_sub_ps(_mm_mul_ps(vec1, VecSwizzle(vec2, 3, 0, 3, 0)),
		                  _mm_mul_ps(VecSwizzle(vec1, 1, 0, 3, 2), VecSwizzle(vec2, 2, 1, 2, 1)));
	}

	// this will not work on camera matrix this is for only transformation matricies
	inline Matrix4 static VECTORCALL InverseTransform(const Matrix4 inM) noexcept
	{
		Matrix4 out;
		// transpose 3x3, we know m03 = m13 = m23 = 0
		__m128 t0 = VecShuffle_0101(inM.r[0], inM.r[1]); // 00, 01, 10, 11
		__m128 t1 = VecShuffle_2323(inM.r[0], inM.r[1]); // 02, 03, 12, 13
		out.r[0]  = VecShuffle(t0, inM.r[2], 0,2,0,3); // 00, 10, 20, 23(=0)
		out.r[1]  = VecShuffle(t0, inM.r[2], 1,3,1,3); // 01, 11, 21, 23(=0)
		out.r[2]  = VecShuffle(t1, inM.r[2], 0,2,2,3); // 02, 12, 22, 23(=0)

		// (SizeSqr(mVec[0]), SizeSqr(mVec[1]), SizeSqr(mVec[2]), 0)
		__m128 sizeSqr;
		sizeSqr =                     _mm_mul_ps(out.r[0], out.r[0]);
		sizeSqr = _mm_add_ps(sizeSqr, _mm_mul_ps(out.r[1], out.r[1]));
		sizeSqr = _mm_add_ps(sizeSqr, _mm_mul_ps(out.r[2], out.r[2]));

		// optional test to avoid divide by 0
		static __m128 const one = _mm_set1_ps(1.f);
		// for each component, if(sizeSqr < SMALL_NUMBER) sizeSqr = 1;
		__m128 rSizeSqr = _mm_blendv_ps(
			_mm_div_ps(one, sizeSqr), one,
			_mm_cmplt_ps(sizeSqr, _mm_set1_ps(1.e-8f))
		);

		out.r[0] = _mm_mul_ps(out.r[0], rSizeSqr);
		out.r[1] = _mm_mul_ps(out.r[1], rSizeSqr);
		out.r[2] = _mm_mul_ps(out.r[2], rSizeSqr);
		// last line
		out.r[3] =                      _mm_mul_ps(out.r[0], VecSwizzle1(inM.r[3], 0));
		out.r[3] = _mm_add_ps(out.r[3], _mm_mul_ps(out.r[1], VecSwizzle1(inM.r[3], 1)));
		out.r[3] = _mm_add_ps(out.r[3], _mm_mul_ps(out.r[2], VecSwizzle1(inM.r[3], 2)));
		out.r[3] = _mm_sub_ps(_mm_setr_ps(0.f, 0.f, 0.f, 1.f), out.r[3]);
		return out;
	}

	inline Matrix4 static VECTORCALL Inverse(const Matrix4 inM) noexcept
	{
		__m128 A = VecShuffle_0101(inM.r[0], inM.r[1]);
		__m128 B = VecShuffle_2323(inM.r[0], inM.r[1]);
		__m128 C = VecShuffle_0101(inM.r[2], inM.r[3]);
		__m128 D = VecShuffle_2323(inM.r[2], inM.r[3]);

		__m128 detSub = _mm_sub_ps(
			_mm_mul_ps(VecShuffle(inM.r[0], inM.r[2], 0, 2, 0, 2), VecShuffle(inM.r[1], inM.r[3], 1, 3, 1, 3)),
			_mm_mul_ps(VecShuffle(inM.r[0], inM.r[2], 1, 3, 1, 3), VecShuffle(inM.r[1], inM.r[3], 0, 2, 0, 2))
		);
		__m128 detA = VecSwizzle1(detSub, 0);
		__m128 detB = VecSwizzle1(detSub, 1);
		__m128 detC = VecSwizzle1(detSub, 2);
		__m128 detD = VecSwizzle1(detSub, 3);

		__m128 D_C  = Mat2AdjMul(D, C);
		__m128 A_B  = Mat2AdjMul(A, B);
		__m128 X_   = _mm_sub_ps(_mm_mul_ps(detD, A), Mat2Mul(B, D_C));
		__m128 W_   = _mm_sub_ps(_mm_mul_ps(detA, D), Mat2Mul(C, A_B));

		__m128 detM = _mm_mul_ps(detA, detD);

		__m128 Y_   = _mm_sub_ps(_mm_mul_ps(detB, C), Mat2MulAdj(D, A_B));
		__m128 Z_   = _mm_sub_ps(_mm_mul_ps(detC, B), Mat2MulAdj(A, D_C));

		detM = _mm_add_ps(detM, _mm_mul_ps(detB, detC));

		__m128 tr = _mm_mul_ps(A_B, VecSwizzle(D_C, 0, 2, 1, 3));
		tr   = _mm_hadd_ps(tr, tr);
		tr   = _mm_hadd_ps(tr, tr);
		detM = _mm_sub_ps(detM, tr);

		const __m128 adjSignMask = _mm_setr_ps(1.f, -1.f, -1.f, 1.f);
		__m128 rDetM = _mm_div_ps(adjSignMask, detM);

		X_ = _mm_mul_ps(X_, rDetM);
		Y_ = _mm_mul_ps(Y_, rDetM);
		Z_ = _mm_mul_ps(Z_, rDetM);
		W_ = _mm_mul_ps(W_, rDetM);

		Matrix4 out;
		out.r[0] = VecShuffle(X_, Y_, 3, 1, 3, 1);
		out.r[1] = VecShuffle(X_, Y_, 2, 0, 2, 0);
		out.r[2] = VecShuffle(Z_, W_, 3, 1, 3, 1);
		out.r[3] = VecShuffle(Z_, W_, 2, 0, 2, 0);
		return out;
	}

	inline Matrix4 static VECTORCALL Multiply(const Matrix4 in1, const Matrix4& in2)
	{
		Matrix4 out;
		__m128 e0 = _mm_shuffle_ps(in2[0], in2[0], _mm_shuffle(0, 0, 0, 0));
		__m128 e1 = _mm_shuffle_ps(in2[0], in2[0], _mm_shuffle(1, 1, 1, 1));
		__m128 e2 = _mm_shuffle_ps(in2[0], in2[0], _mm_shuffle(2, 2, 2, 2));
		__m128 e3 = _mm_shuffle_ps(in2[0], in2[0], _mm_shuffle(3, 3, 3, 3));
		__m128 m0 = _mm_mul_ps(in1[0], e0);
		__m128 m1 = _mm_mul_ps(in1[1], e1);
		__m128 m2 = _mm_mul_ps(in1[2], e2);
		__m128 m3 = _mm_mul_ps(in1[3], e3);
		__m128 a0 = _mm_add_ps(m0, m1);
		__m128 a1 = _mm_add_ps(m2, m3);
		__m128 a2 = _mm_add_ps(a0, a1);
		out[0] = a2;
		
		e0 = _mm_shuffle_ps(in2[1], in2[1], _mm_shuffle(0, 0, 0, 0));
		e1 = _mm_shuffle_ps(in2[1], in2[1], _mm_shuffle(1, 1, 1, 1));
		e2 = _mm_shuffle_ps(in2[1], in2[1], _mm_shuffle(2, 2, 2, 2));
		e3 = _mm_shuffle_ps(in2[1], in2[1], _mm_shuffle(3, 3, 3, 3));
		m0 = _mm_mul_ps(in1[0], e0);
		m1 = _mm_mul_ps(in1[1], e1);
		m2 = _mm_mul_ps(in1[2], e2);
		m3 = _mm_mul_ps(in1[3], e3);
		a0 = _mm_add_ps(m0, m1);
		a1 = _mm_add_ps(m2, m3);
		a2 = _mm_add_ps(a0, a1);
		out[1] = a2;
		
		e0 = _mm_shuffle_ps(in2[2], in2[2], _mm_shuffle(0, 0, 0, 0));
		e1 = _mm_shuffle_ps(in2[2], in2[2], _mm_shuffle(1, 1, 1, 1));
		e2 = _mm_shuffle_ps(in2[2], in2[2], _mm_shuffle(2, 2, 2, 2));
		e3 = _mm_shuffle_ps(in2[2], in2[2], _mm_shuffle(3, 3, 3, 3));
		m0 = _mm_mul_ps(in1[0], e0);
		m1 = _mm_mul_ps(in1[1], e1);
		m2 = _mm_mul_ps(in1[2], e2);
		m3 = _mm_mul_ps(in1[3], e3);
		a0 = _mm_add_ps(m0, m1);
		a1 = _mm_add_ps(m2, m3);
		a2 = _mm_add_ps(a0, a1);
		out[2] = a2;
		
		e0 = _mm_shuffle_ps(in2[3], in2[3], _mm_shuffle(0, 0, 0, 0));
		e1 = _mm_shuffle_ps(in2[3], in2[3], _mm_shuffle(1, 1, 1, 1));
		e2 = _mm_shuffle_ps(in2[3], in2[3], _mm_shuffle(2, 2, 2, 2));
		e3 = _mm_shuffle_ps(in2[3], in2[3], _mm_shuffle(3, 3, 3, 3));
		m0 = _mm_mul_ps(in1[0], e0);
		m1 = _mm_mul_ps(in1[1], e1);
		m2 = _mm_mul_ps(in1[2], e2);
		m3 = _mm_mul_ps(in1[3], e3);
		a0 = _mm_add_ps(m0, m1);
		a1 = _mm_add_ps(m2, m3);
		a2 = _mm_add_ps(a0, a1);
		out[3] = a2;
		return out;
	}

	FINLINE static Matrix4 PositionRotationScale(const Vector3f& position, const Quaternion& rotation, const Vector3f& scale)
	{
		Matrix4 result(ForceInit);
		result *= FromPosition(position);
		result *= FromQuaternion(rotation);
		result *= CreateScale(position);
		return result;
	}

	FINLINE static Vector3f VECTORCALL ExtractPosition(const Matrix4 matrix) noexcept
	{
		Vector3f res;
		_mm_storeu_ps(&res.x, matrix.r[3]);
		return res;
	}

	FINLINE static Vector3f VECTORCALL ExtractScale(const Matrix4 matrix) noexcept
	{
		return Vector3f(SSEVectorLengthf(matrix.r[0]), SSEVectorLengthf(matrix.r[2]), SSEVectorLengthf(matrix.r[1]));
	}

	FINLINE static Matrix4 RotationX(float angleRadians) {
		Matrix4 out_matrix(ForceInit);
		float s, c;
		SinCos(angleRadians, &s, &c);
		out_matrix.m[1][1] = c;
		out_matrix.m[1][2] = s;
		out_matrix.m[2][1] = -s;
		out_matrix.m[2][2] = c;
		return out_matrix;
	}

	FINLINE static Matrix4 RotationY(float angleRadians) {
		Matrix4 out_matrix(ForceInit);
		float s, c;
		SinCos(angleRadians, &s, &c);
		out_matrix.m[0][0] = c;
		out_matrix.m[0][2] = -s;
		out_matrix.m[2][0] = s;
		out_matrix.m[2][2] = c;
		return out_matrix;
	}
	
	FINLINE static Matrix4 RotationZ(float angleRadians) {
		Matrix4 out_matrix(ForceInit);
		float s, c;
		SinCos(angleRadians, &s, &c);
		
		out_matrix.m[0][0] = c;
		out_matrix.m[0][1] = s;
		out_matrix.m[1][0] = -s;
		out_matrix.m[1][1] = c;
		return out_matrix;
	}

	FINLINE static Matrix4 RotationFromEuler(Vector3f eulerRadians) {
		return RotationX(eulerRadians.x) * RotationY(eulerRadians.y) * RotationZ(eulerRadians.z);
	}

	FINLINE Matrix4 static VECTORCALL LookAt(Vector3f eyePosition, Vector3f focusPosition, Vector3f upDirection) noexcept
	{
		__m128 EyePosition = _mm_loadu_ps(&eyePosition.x);
		__m128 FocusPosition = _mm_loadu_ps(&focusPosition.x);
		__m128 UpDirection = _mm_loadu_ps(&upDirection.x);

		// negate because right handed, for left handed no need to negate
		__m128 EyeDirection = _mm_sub_ps(EyePosition, FocusPosition);

		__m128  R2 = SSEVectorNormalize(EyeDirection);
		__m128  R0 = SSEVector3Cross(UpDirection, R2);
		R0 = SSEVectorNormalize(R0);
		__m128  R1 = SSEVector3Cross(R2, R0);
		__m128  NegEyePosition = _mm_mul_ps(EyePosition, g_XMNegativeOne);
		__m128 D0 = SSEVector3Dot(R0, NegEyePosition);
		__m128 D1 = SSEVector3Dot(R1, NegEyePosition);
		__m128 D2 = SSEVector3Dot(R2, NegEyePosition);
		R0 = _mm_and_ps(R0, g_XMMask3);
		R1 = _mm_and_ps(R1, g_XMMask3);
		R2 = _mm_and_ps(R2, g_XMMask3);
		D0 = _mm_and_ps(D0, g_XMMaskW);
		D1 = _mm_and_ps(D1, g_XMMaskW);
		D2 = _mm_and_ps(D2, g_XMMaskW);
		D0 = _mm_or_ps(D0, R0);
		D1 = _mm_or_ps(D1, R1);
		D2 = _mm_or_ps(D2, R2);
		Matrix4 M;
		M.r[0] = D0;
		M.r[1] = D1;
		M.r[2] = D2;
		M.r[3] = g_XMIdentityR3;
		M = Matrix4::Transpose(M);
		return M;
	}

	static Matrix3 VECTORCALL ConvertToMatrix3(const Matrix4 M)
	{
		Matrix3 result;
		SSEStoreVector3(&result.x.x, M.r[0]);
		SSEStoreVector3(&result.y.x, M.r[1]);
		SSEStoreVector3(&result.z.x, M.r[2]);
		return result;
	}

	static Quaternion VECTORCALL ExtractRotation(const Matrix4 M, bool rowNormalize = true) 
	{
		int i, j, k = 0;
		float root, trace = M.m[0][0] + M.m[1][1] + M.m[2][2];
		Quaternion Orientation;

		if (trace > 0.0f)
		{
			root = Sqrt(trace + 1.0f);
			Orientation.w = 0.5f * root;
			root = 0.5f / root;
			Orientation.x = root * (M.m[1][2] - M.m[2][1]);
			Orientation.y = root * (M.m[2][0] - M.m[0][2]);
			Orientation.z = root * (M.m[0][1] - M.m[1][0]);
		}
		else
		{
			static int Next[3] = { 1, 2, 0 };
			i = 0;
			if (M.m[1][1] > M.m[0][0]) i = 1;
			if (M.m[2][2] > M.m[i][i]) i = 2;
			j = Next[i];
			k = Next[j];

			root = Sqrt(M.m[i][i] - M.m[j][j] - M.m[k][k] + 1.0f);

			Orientation[i] = 0.5f * root;
			root = 0.5f / root;
			Orientation[j] = root * (M.m[i][j] + M.m[j][i]);
			Orientation[k] = root * (M.m[i][k] + M.m[k][i]);
			Orientation.w  = root * (M.m[j][k] - M.m[k][j]);
		} 
		return Orientation;
	}
	static Matrix4 VECTORCALL FromQuaternion(const Quaternion quaternion)
	{
		Matrix4 M;
		static const Vector432F  Constant1110 = {1.0f, 1.0f, 1.0f, 0.0f};
		__m128 q = quaternion.vec;

		__m128 Q0 = _mm_add_ps(q, q);
		__m128 Q1 = _mm_mul_ps(q,Q0);

		__m128 V0 = _mm_shuffle_ps(Q1,Q1,_mm_shuffle(3,0,0,1));
		V0 = _mm_and_ps(V0,g_XMMask3);
		__m128 V1 = _mm_shuffle_ps(Q1,Q1,_mm_shuffle(3,1,2,2));
		V1 = _mm_and_ps(V1,g_XMMask3);
		__m128  R0 = _mm_sub_ps(Constant1110,V0);
		R0 = _mm_sub_ps(R0, V1);

		V0 = _mm_shuffle_ps(q, q,_mm_shuffle(3,1,0,0));
		V1 = _mm_shuffle_ps(Q0,Q0,_mm_shuffle(3,2,1,2));
		V0 = _mm_mul_ps(V0, V1);

		V1 = _mm_shuffle_ps(q, q,_mm_shuffle(3,3,3,3));
		__m128 V2 = _mm_shuffle_ps(Q0,Q0,_mm_shuffle(3,0,2,1));
		V1 = _mm_mul_ps(V1, V2);

		__m128 R1 = _mm_add_ps(V0, V1);
		__m128 R2 = _mm_sub_ps(V0, V1);

		V0 = _mm_shuffle_ps(R1,R2,_mm_shuffle(1,0,2,1));
		V0 = _mm_shuffle_ps(V0,V0,_mm_shuffle(1,3,2,0));
		V1 = _mm_shuffle_ps(R1,R2,_mm_shuffle(2,2,0,0));
		V1 = _mm_shuffle_ps(V1,V1,_mm_shuffle(2,0,2,0));

		Q1 = _mm_shuffle_ps(R0,V0,_mm_shuffle(1,0,3,0));
		Q1 = _mm_shuffle_ps(Q1,Q1,_mm_shuffle(1,3,2,0));
		M.r[0] = Q1;
		Q1 = _mm_shuffle_ps(R0,V0,_mm_shuffle(3,2,3,1));
		Q1 = _mm_shuffle_ps(Q1,Q1,_mm_shuffle(1,3,0,2));
		M.r[1] = Q1;
		Q1 = _mm_shuffle_ps(V1,R0,_mm_shuffle(3,2,1,0));
		M.r[2] = Q1;
		M.r[3] = g_XMIdentityR3;
		return M;
	}

	FINLINE static Matrix4 VECTORCALL Transpose(const Matrix4 M)
	{
		const __m128 vTemp1 = _mm_shuffle_ps(M.r[0], M.r[1], _mm_shuffle(1, 0, 1, 0));
		const __m128 vTemp3 = _mm_shuffle_ps(M.r[0], M.r[1], _mm_shuffle(3, 2, 3, 2));
		const __m128 vTemp2 = _mm_shuffle_ps(M.r[2], M.r[3], _mm_shuffle(1, 0, 1, 0));
		const __m128 vTemp4 = _mm_shuffle_ps(M.r[2], M.r[3], _mm_shuffle(3, 2, 3, 2));
		Matrix4 mResult;
		mResult.r[0] = _mm_shuffle_ps(vTemp1, vTemp2, _mm_shuffle(2, 0, 2, 0));
		mResult.r[1] = _mm_shuffle_ps(vTemp1, vTemp2, _mm_shuffle(3, 1, 3, 1));
		mResult.r[2] = _mm_shuffle_ps(vTemp3, vTemp4, _mm_shuffle(2, 0, 2, 0));
		mResult.r[3] = _mm_shuffle_ps(vTemp3, vTemp4, _mm_shuffle(3, 1, 3, 1));
		return mResult;
	}

	FINLINE static __m128 VECTORCALL Vector3Transform(const Vector3f V, const Matrix4& M)
	{
		__m128 vec = _mm_loadu_ps(&V.x);
		__m128 vResult = _mm_shuffle_ps(vec, vec, _mm_shuffle(0, 0, 0, 0));
		vResult = _mm_mul_ps(vResult, M.r[0]);
		__m128 vTemp = _mm_shuffle_ps(vec, vec, _mm_shuffle(1, 1, 1, 1));
		vTemp = _mm_mul_ps(vTemp, M.r[1]);
		vResult = _mm_add_ps(vResult, vTemp);
		vTemp = _mm_shuffle_ps(vec, vec, _mm_shuffle(2, 2, 2, 2));
		vTemp = _mm_mul_ps(vTemp, M.r[2]);
		vResult = _mm_add_ps(vResult, vTemp);
		vResult = _mm_add_ps(vResult, M.r[3]);
		return vResult;
	}

	FINLINE static Vector4f VECTORCALL Vector4Transform(Vector4f v, const Matrix4& m)
	{
		__m128 v0 = _mm_mul_ps(m.r[0], _mm_permute_ps(v.vec, _mm_shuffle(0, 0, 0, 0)));
		__m128 v1 = _mm_mul_ps(m.r[1], _mm_permute_ps(v.vec, _mm_shuffle(1, 1, 1, 1)));
		__m128 v2 = _mm_mul_ps(m.r[2], _mm_permute_ps(v.vec, _mm_shuffle(2, 2, 2, 2)));
		__m128 v3 = _mm_mul_ps(m.r[3], _mm_permute_ps(v.vec, _mm_shuffle(3, 3, 3, 3)));
		__m128 a0 = _mm_add_ps(v0, v1);
		__m128 a1 = _mm_add_ps(v2, v3);
		__m128 a2 = _mm_add_ps(a0, a1);
		return a2;
	}
};
 
FINLINE static __m128 VECTORCALL Vector4Transform(__m128 v, const Matrix4& m)
{
	__m128 v0 = _mm_mul_ps(m.r[0], _mm_permute_ps(v, _mm_shuffle(0, 0, 0, 0)));
	__m128 v1 = _mm_mul_ps(m.r[1], _mm_permute_ps(v, _mm_shuffle(1, 1, 1, 1)));
	__m128 v2 = _mm_mul_ps(m.r[2], _mm_permute_ps(v, _mm_shuffle(2, 2, 2, 2)));
	__m128 v3 = _mm_mul_ps(m.r[3], _mm_permute_ps(v, _mm_shuffle(3, 3, 3, 3)));
	__m128 a0 = _mm_add_ps(v0, v1);
	__m128 a1 = _mm_add_ps(v2, v3);
	__m128 a2 = _mm_add_ps(a0, a1);
	return a2;
}

FINLINE void InitializeMatrix4(Matrix4& mat, float s) 
{
	mat.r[0] = mat.r[1] = mat.r[2] = mat.r[3] = _mm_set_ps1(s);
}

FINLINE void VECTORCALL InitializeMatrix4(Matrix4& r, __m128 x, __m128 y, const __m128& z, const __m128& w)
{
	r[0] = x; r[1] = y; r[2] = z; r[3] = w;
}

#else // sse is not supported

// todo create noın sse matrix
AX_ALIGNED(16) struct Matrix4
{
	struct
	{
		float m[4][4];
	};

	const __m128& operator [] (int index) const { return r[index]; }
	__m128& operator [] (int index) { return r[index]; }
	
	Vector4f VECTORCALL  operator *  (const Vector3f v)  noexcept { return Vector3Transform(v, *this); };
	Vector4f VECTORCALL  operator *  (const Vector4f& v)   noexcept { return Vector4Transform(v, *this); };

	Matrix4 VECTORCALL  operator *  (const Matrix4& M)  noexcept { return Matrix4::Multiply(*this, M); };
	Matrix4& VECTORCALL operator *= (const Matrix4& M)  noexcept { *this = Matrix4::Multiply(*this, M); return *this; };

	Matrix4() {}
	explicit Matrix4(EForceInit) 
	{
		r[0] = g_XMIdentityR0;
		r[1] = g_XMIdentityR1;
		r[2] = g_XMIdentityR2;
		r[3] = g_XMIdentityR3;
	}

	explicit Matrix4(float s)
	{
		r[0] = r[1] = r[2] = r[3] = _mm_set_ps1(s);
	}

	FINLINE static Matrix4 Identity()
	{
		Matrix4 M;
		M.r[0] = g_XMIdentityR0;
		M.r[1] = g_XMIdentityR1;
		M.r[2] = g_XMIdentityR2;
		M.r[3] = g_XMIdentityR3;
		return M;
	}

	FINLINE static Matrix4 FromPosition(const float x, const float y, const float z)
	{
		Matrix4 M;
		M.r[0] = g_XMIdentityR0;
		M.r[1] = g_XMIdentityR1;
		M.r[2] = g_XMIdentityR2;
		M.r[3] = _mm_set_ps(1.0f, z, y, x);
		return M;
	}

	FINLINE static Matrix4 FromPosition(const Vector3f& vec3)
	{
		return FromPosition(vec3.x, vec3.y, vec3.z);
	}

	FINLINE static Matrix4 CreateScale(const float ScaleX, const float ScaleY, const float ScaleZ)
	{
		Matrix4 M;
		M.r[0] = _mm_set_ps(0.0f, 0.0f, 0.0f, ScaleX);
		M.r[1] = _mm_set_ps(0.0f, 0.0f, ScaleY, 0.0f);
		M.r[2] = _mm_set_ps(0.0f, ScaleZ, 0.0f, 0.0f);
		M.r[3] = g_XMIdentityR3;
		return M;
	}

	FINLINE static Matrix4 CreateScale(const Vector3f& vec3)
	{
		return CreateScale(vec3.x, vec3.y, vec3.z);
	}

	// please assign normalized vectors, returns view matrix
	FINLINE static Matrix4 VECTORCALL LookAtRH(Vector3f eye, Vector3f center, const Vector3f& up)
	{
		return ;
	}

	FINLINE static Matrix4 PerspectiveFovRH(float fov, float width, float height, float zNear, float zFar)
	{
		return M;
	}

	inline Matrix4 static VECTORCALL InverseTransform(const Matrix4 inM) noexcept
	{
		return out;
	}

	inline Matrix4 static VECTORCALL Inverse(const Matrix4 inM) noexcept
	{
		return out;
	}

	inline Matrix4 static VECTORCALL Multiply(const Matrix4 in1, const Matrix4& in2)
	{
		return out;
	}

	FINLINE static Matrix4 PositionRotationScale(const Vector3f& position, const Quaternion& rotation, const Vector3f& scale)
	{
		Matrix4 result(ForceInit);
		result *= FromPosition(position);
		result *= FromQuaternion(rotation);
		result *= CreateScale(position);
		return result;
	}

	FINLINE static Vector3f VECTORCALL ExtractPosition(const Matrix4 matrix) noexcept
	{
		return res;
	}

	FINLINE static Vector3f VECTORCALL ExtractScale(const Matrix4 matrix) noexcept
	{
		return Vector3f(SSEVectorLengthf(matrix.r[0]), SSEVectorLengthf(matrix.r[2]), SSEVectorLengthf(matrix.r[1]));
	}

	FINLINE static Matrix4 RotationX(float angleRadians) {
		Matrix4 out_matrix(ForceInit);
		float s, c;
		SinCos(angleRadians, &s, &c);
		
		out_matrix.m[1][1] = c;
		out_matrix.m[1][2] = s;
		out_matrix.m[2][1] = -s;
		out_matrix.m[2][2] = c;
		return out_matrix;
	}

	FINLINE static Matrix4 RotationY(float angleRadians) {
		Matrix4 out_matrix(ForceInit);
		float s, c;
		SinCos(angleRadians, &s, &c);
		out_matrix.m[0][0] = c;
		out_matrix.m[0][2] = -s;
		out_matrix.m[2][0] = s;
		out_matrix.m[2][2] = c;
		return out_matrix;
	}
	
	FINLINE static Matrix4 RotationZ(float angleRadians) {
		Matrix4 out_matrix(ForceInit);
		float s, c;
		SinCos(angleRadians, &s, &c);
		out_matrix.m[0][0] = c;
		out_matrix.m[0][1] = s;
		out_matrix.m[1][0] = -s;
		out_matrix.m[1][1] = c;
		return out_matrix;
	}

	FINLINE static Matrix4 RotationFromEuler(Vector3f eulerRadians) {
		return RotationX(eulerRadians.x) * RotationY(eulerRadians.y) * RotationZ(eulerRadians.z);
	}

	FINLINE Matrix4 static VECTORCALL LookAt(Vector3f eyePosition, Vector3f focusPosition, Vector3f upDirection) noexcept
	{
		return M;
	}

	static Matrix3 VECTORCALL ConvertToMatrix3(const Matrix4 M)
	{
		Matrix3 result;
		SSEStoreVector3(&result.x.x, M.r[0]);
		SSEStoreVector3(&result.y.x, M.r[1]);
		SSEStoreVector3(&result.z.x, M.r[2]);
		return result;
	}

	static Quaternion VECTORCALL ExtractRotation(const Matrix4 M, bool rowNormalize = true) noexcept
	{
		int i, j, k = 0;
		float root, trace = M.m[0][0] + M.m[1][1] + M.m[2][2];
		Quaternion Orientation;

		if (trace > 0.0f)
		{
			root = Sqrt(trace + 1.0f);
			Orientation.w = 0.5f * root;
			root = 0.5f / root;
			Orientation.x = root * (M.m[1][2] - M.m[2][1]);
			Orientation.y = root * (M.m[2][0] - M.m[0][2]);
			Orientation.z = root * (M.m[0][1] - M.m[1][0]);
		}
		else
		{
			static int Next[3] = { 1, 2, 0 };
			i = 0;
			if (M.m[1][1] > M.m[0][0]) i = 1;
			if (M.m[2][2] > M.m[i][i]) i = 2;
			j = Next[i];
			k = Next[j];

			root = Sqrt(M.m[i][i] - M.m[j][j] - M.m[k][k] + 1.0f);

			Orientation[i] = 0.5f * root;
			root = 0.5f / root;
			Orientation[j] = root * (M.m[i][j] + M.m[j][i]);
			Orientation[k] = root * (M.m[i][k] + M.m[k][i]);
			Orientation.w  = root * (M.m[j][k] - M.m[k][j]);
		} 
		return Orientation;
	}

	static Matrix4 VECTORCALL FromQuaternion(const Quaternion quaternion)
	{
		return M;
	}

	FINLINE static Matrix4 VECTORCALL Transpose(const Matrix4 M)
	{
		return mResult;
	}

	FINLINE static __m128 VECTORCALL Vector3Transform(const Vector3f V, const Matrix4& M) noexcept
	{
		return vResult;
	}

	FINLINE static Vector4f VECTORCALL Vector4Transform(Vector4f v, const Matrix4& m)
	{
		return a2;
	}
};
 
FINLINE static Vector4f VECTORCALL Vector4Transform(Vector4f v, const Matrix4& m)
{
	return a2;
}

FINLINE void InitializeMatrix4(Matrix4& mat, float s) 
{

}

FINLINE void InitializeMatrix4(Matrix4& mat, Vector4f a, Vector4f b, Vector4f c, Vector4f d) 
{

}

#endif // AX_SUPPORT_SSE2