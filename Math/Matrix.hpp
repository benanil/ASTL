#pragma once
#include "Quaternion.hpp"

AX_NAMESPACE 

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

	__forceinline static Matrix3 Make(float x, float y, float z,
	                                  float a, float b, float c,
	                                  float u, float v, float s)
	{
		Matrix3 M;
		M.m[0][0] = x; M.m[0][1] = y; M.m[0][2] = z; 
		M.m[1][0] = a; M.m[1][1] = b; M.m[1][2] = c; 
		M.m[2][0] = u; M.m[2][1] = v; M.m[2][2] = s; 
		return M;
	}

	__forceinline static Matrix3 Identity()
	{
		return Make(1.0f, 0.0f, 0.0f,
		            0.0f, 1.0f, 0.0f,
		            0.0f, 0.0f, 1.0f);
	}

	inline static Matrix3 LookAt(Vector3f direction, Vector3f up)
	{
		Matrix3 result;
		result.vec[2] = direction;
		Vector3f const& Right = Vector3f::Cross(up, result.vec[2]);
		result.vec[0] = Right * RSqrt(MAX(0.00001f, Vector3f::Dot(Right, Right)));
		result.vec[1] = Vector3f::Cross(result.vec[2], result.vec[0]);
		return result;
	}

	inline static Matrix3 Multiply(const Matrix3& a, const Matrix3& b)
	{
		Matrix3 result;
		float3 vx = b.x * a.x.x;
		float3 vy = b.y * a.x.y;
		float3 vz = b.z * a.x.z;

		vx += vz; vx += vy;
		result.x = vx;

		vx = b.x * a.y.x;
		vy = b.y * a.y.y;
		vz = b.z * a.y.z;
		vx += vz; vx += vy;
		result.y = vx;

		vx = b.x * a.z.x;
		vy = b.y * a.z.y;
		vz = b.z * a.z.z;
		vx += vz; vx += vy;
		result.z = vx;
		return result;
	}

	inline static float3 Multiply(const Matrix3& m, const float3& v) {
		return m.x * v.x + m.y * v.y + m.z * v.z;
	}

	inline static Matrix3 FromQuaternion(const Quaternion& quat)
	{
		const float num9 = quat.x * quat.x, num8 = quat.y * quat.y, num7 = quat.z * quat.z, 
		            num6 = quat.x * quat.y, num5 = quat.z * quat.w, num4 = quat.z * quat.x,
		            num3 = quat.y * quat.w, num2 = quat.y * quat.z, num  = quat.x * quat.w;
		return Make(
			1.0f - (2.0f * (num8 + num7)), 2.0f * (num6 + num5), 2.0f * (num4 - num3),
			2.0f * (num6 - num5), 1.0f - (2.0f * (num7 + num9)), 2.0f * (num2 + num) ,
			2.0f * (num4 + num3), 2.0f * (num2 - num), 1.0f - (2.0f * (num8 + num9)) 
		);
	}

	Quaternion ToQuaternion() const
	{
		Quaternion Orientation;
		QuaternionFromMatrix(&Orientation.x, vec);
		return Orientation;
	}
};

// todo make non simd version
#ifdef AX_SUPPORT_SSE

__forceinline static __m128 VECTORCALL Vector4Transform(__m128 v, const __m128 r[4])
{
	return _mm_add_ps(
		_mm_add_ps(_mm_mul_ps(r[0], VecSwizzle1(v, 0)), _mm_mul_ps(r[1], VecSwizzle1(v, 1))),
		_mm_add_ps(_mm_mul_ps(r[2], VecSwizzle1(v, 2)), _mm_mul_ps(r[3], VecSwizzle1(v, 3)))
	);
}

struct alignas(16) Matrix4
{
	union
	{
		struct { __m128   r[4]; };
		struct { Vector4f v[4]; };
		struct { float m[4][4]; };
	};

	const Vector4f& operator [] (int index) const { return v[index]; }
	Vector4f& operator [] (int index) { return v[index]; }

	Vector3f VECTORCALL  operator * (const Vector3f v)  noexcept { Vector4f x; x.vec = Vector3Transform(v, *this); return x.xyz(); };
	Vector4f VECTORCALL  operator * (const Vector4f& v) noexcept { Vector4f x; x.vec = Vector4Transform(v.vec, r); return x; };

	Matrix4 VECTORCALL  operator *  (const Matrix4& M)  noexcept { return Matrix4::Multiply(*this, M); };
	Matrix4& VECTORCALL operator *= (const Matrix4& M)  noexcept { *this = Matrix4::Multiply(*this, M); return *this; };

	__forceinline static Matrix4 Make(float x, float y, float z, float w,
	                                  float a, float b, float c, float d,
	                                  float u, float v, float s, float t,
	                                  float j, float k, float l, float m)
	{
		Matrix4 M;
		M.r[0] = _mm_set_ps(x, y, z, w);
		M.r[1] = _mm_set_ps(a, b, c, d);
		M.r[2] = _mm_set_ps(u, v, s, t);
		M.r[3] = _mm_set_ps(j, k, l, m);
		return M;
	}

	__forceinline static Matrix4 Identity()
	{
		Matrix4 M;
		M.r[0] = g_XMIdentityR0;
		M.r[1] = g_XMIdentityR1;
		M.r[2] = g_XMIdentityR2;
		M.r[3] = g_XMIdentityR3;
		return M;
	}

	__forceinline static Matrix4 FromPosition(const float x, const float y, const float z)
	{
		Matrix4 M;
		M.r[0] = g_XMIdentityR0;
		M.r[1] = g_XMIdentityR1;
		M.r[2] = g_XMIdentityR2;
		M.r[3] = _mm_set_ps(1.0f, z, y, x);
		return M;
	}

	__forceinline static Matrix4 FromPosition(const Vector3f& vec3)
	{
		return FromPosition(vec3.x, vec3.y, vec3.z);
	}

	__forceinline static Matrix4 CreateScale(const float ScaleX, const float ScaleY, const float ScaleZ)
	{
		Matrix4 M;
		M.r[0] = _mm_set_ps(0.0f, 0.0f, 0.0f, ScaleX);
		M.r[1] = _mm_set_ps(0.0f, 0.0f, ScaleY, 0.0f);
		M.r[2] = _mm_set_ps(0.0f, ScaleZ, 0.0f, 0.0f);
		M.r[3] = g_XMIdentityR3;
		return M;
	}

	__forceinline static Matrix4 CreateScale(Vector3f vec3)
	{
		return CreateScale(vec3.x, vec3.y, vec3.z);
	}

	__forceinline static Matrix4 CreateRotation(Vector3f right, Vector3f up, Vector3f forward)
	{
		Matrix4 m;
		m.r[0] = _mm_and_ps(g_XMMask3, _mm_loadu_ps(&right.x));
		m.r[1] = _mm_and_ps(g_XMMask3, _mm_loadu_ps(&up.x));
		m.r[2] = _mm_and_ps(g_XMMask3, _mm_loadu_ps(&forward.x));
		m.r[3] = _mm_setr_ps(0.0f, 0.0f, 0.0f, 1.0f);
		return m;
	}

	const Vector3f& GetForward() const { return v[2].xyz(); }
	const Vector3f& GetUp()      const { return v[1].xyz(); }
	const Vector3f& GetRight()   const { return v[0].xyz(); }

	// please assign normalized vectors, returns view matrix
	__forceinline static Matrix4 VECTORCALL LookAtRH(Vector3f eye, Vector3f center, const Vector3f& up)
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
		M.r[0] = SSESelect(D0, R0, g_XMSelect1110); // no need select ?
		M.r[1] = SSESelect(D1, R1, g_XMSelect1110);
		M.r[2] = SSESelect(D2, EyeDirection, g_XMSelect1110);
		M.r[3] = g_XMIdentityR3;
		return Matrix4::Transpose(M);
	}

	__forceinline static Matrix4 PerspectiveFovRH(float fov, float width, float height, float zNear, float zFar)
	{
		const float rad = fov;
		const float h = Cos(0.5f * rad) / Sin(0.5f * rad);
		const float w = h * height / width; /// max(width , Height) / min(width , Height)?
		Matrix4 M = Identity();
		M.m[0][0] = w;
		M.m[1][1] = h;
		M.m[2][2] = -(zFar + zNear) / (zFar - zNear);
		M.m[2][3] = -1.0f;
		M.m[3][2] = -(2.0f * zFar * zNear) / (zFar - zNear);
		M.m[3][3] = 0.0f;
		return M;
	}

	// https://lxjk.github.io/2017/09/03/Fast-4x4-Matrix-Inverse-with-SSE-SIMD-Explained.html
	// for row major matrix
    // we use __m128 to represent 2x2 matrix as A = | A0  A1 |
    //                                              | A2  A3 |
    // 2x2 row major Matrix multiply A*B
	static __forceinline __m128 VECTORCALL Mat2Mul(__m128 vec1, __m128 vec2)
	{
		return _mm_add_ps(_mm_mul_ps(vec1, VecSwizzle(vec2, 0, 3, 0, 3)), 
		                  _mm_mul_ps(VecSwizzle(vec1, 1, 0, 3, 2), VecSwizzle(vec2, 2, 1, 2, 1)));
	}
	// 2x2 row major Matrix adjugate multiply (A#)*B
	static __forceinline __m128 VECTORCALL Mat2AdjMul(__m128 vec1, __m128 vec2)
	{
		return _mm_sub_ps(_mm_mul_ps(VecSwizzle(vec1, 3, 3, 0, 0), vec2),
		                  _mm_mul_ps(VecSwizzle(vec1, 1, 1, 2, 2), VecSwizzle(vec2, 2, 3, 0, 1)));

	}
	// 2x2 row major Matrix multiply adjugate A*(B#)
	static __forceinline __m128 VECTORCALL Mat2MulAdj(__m128 vec1, __m128 vec2)
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
		__m128 m0 = _mm_mul_ps(in1.r[0], VecSwizzle1(in2.r[0], 0));
		__m128 m1 = _mm_mul_ps(in1.r[1], VecSwizzle1(in2.r[0], 1));
		__m128 m2 = _mm_mul_ps(in1.r[2], VecSwizzle1(in2.r[0], 2));
		__m128 m3 = _mm_mul_ps(in1.r[3], VecSwizzle1(in2.r[0], 3));
		out.r[0] = _mm_add_ps(_mm_add_ps(m0, m1), _mm_add_ps(m2, m3));

		m0 = _mm_mul_ps(in1.r[0], VecSwizzle1(in2.r[1], 0));
		m1 = _mm_mul_ps(in1.r[1], VecSwizzle1(in2.r[1], 1));
		m2 = _mm_mul_ps(in1.r[2], VecSwizzle1(in2.r[1], 2));
		m3 = _mm_mul_ps(in1.r[3], VecSwizzle1(in2.r[1], 3));
		out.r[1] = _mm_add_ps(_mm_add_ps(m0, m1), _mm_add_ps(m2, m3));

		m0 = _mm_mul_ps(in1.r[0], VecSwizzle1(in2.r[2], 0));
		m1 = _mm_mul_ps(in1.r[1], VecSwizzle1(in2.r[2], 1));
		m2 = _mm_mul_ps(in1.r[2], VecSwizzle1(in2.r[2], 2));
		m3 = _mm_mul_ps(in1.r[3], VecSwizzle1(in2.r[2], 3));
		out.r[2] = _mm_add_ps(_mm_add_ps(m0, m1), _mm_add_ps(m2, m3));

		m0 = _mm_mul_ps(in1.r[0], VecSwizzle1(in2.r[3], 0));
		m1 = _mm_mul_ps(in1.r[1], VecSwizzle1(in2.r[3], 1));
		m2 = _mm_mul_ps(in1.r[2], VecSwizzle1(in2.r[3], 2));
		m3 = _mm_mul_ps(in1.r[3], VecSwizzle1(in2.r[3], 3));
		out.r[3] = _mm_add_ps(_mm_add_ps(m0, m1), _mm_add_ps(m2, m3));
		return out;
	}

	__forceinline static Matrix4 PositionRotationScale(const Vector3f& position, const Quaternion& rotation, const Vector3f& scale)
	{
		Matrix4 result = Identity();
		result *= FromPosition(position);
		result *= FromQuaternion(rotation);
		result *= CreateScale(position);
		return result;
	}

	__forceinline static Vector3f VECTORCALL ExtractPosition(const Matrix4 matrix) noexcept
	{
		Vector3f res;
		_mm_storeu_ps(&res.x, matrix.r[3]);
		return res;
	}

	__forceinline static Vector3f VECTORCALL ExtractScale(const Matrix4 matrix) noexcept
	{
		return MakeVec3(SSEVectorLengthf(matrix.r[0]), SSEVectorLengthf(matrix.r[2]), SSEVectorLengthf(matrix.r[1]));
	}

	__forceinline static Matrix4 RotationX(float angleRadians) {
		Matrix4 out_matrix = Identity();
		float s, c;
		SinCos(angleRadians, &s, &c);
		out_matrix.m[1][1] = c;
		out_matrix.m[1][2] = s;
		out_matrix.m[2][1] = -s;
		out_matrix.m[2][2] = c;
		return out_matrix;
	}

	__forceinline static Matrix4 RotationY(float angleRadians) {
		Matrix4 out_matrix = Identity();
		float s, c;
		SinCos(angleRadians, &s, &c);
		out_matrix.m[0][0] = c;
		out_matrix.m[0][2] = -s;
		out_matrix.m[2][0] = s;
		out_matrix.m[2][2] = c;
		return out_matrix;
	}
	
	__forceinline static Matrix4 RotationZ(float angleRadians) {
		Matrix4 out_matrix = Identity();
		float s, c;
		SinCos(angleRadians, &s, &c);
		out_matrix.m[0][0] = c;
		out_matrix.m[0][1] = s;
		out_matrix.m[1][0] = -s;
		out_matrix.m[1][1] = c;
		return out_matrix;
	}

	__forceinline static Matrix4 RotationFromEuler(Vector3f eulerRadians) 
	{
		return FromQuaternion(Quaternion::FromEuler(eulerRadians));
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
		Quaternion res;
		QuaternionFromMatrix(&res.x, M.c);
		return res;
	}

	static Matrix4 VECTORCALL FromQuaternion(const Quaternion quaternion)
	{
		Matrix4 M;
		const __m128  Constant1110 = _mm_setr_ps(1.0f, 1.0f, 1.0f, 0.0f);
		__m128 q = quaternion.vec;

		__m128 Q0 = _mm_add_ps(q, q);
		__m128 Q1 = _mm_mul_ps(q,Q0);

		__m128 V0 = _mm_shuffle_ps(Q1,Q1,_mm_shuffle(3,0,0,1));
		V0 = _mm_and_ps(V0,g_XMMask3);
		__m128 V1 = _mm_shuffle_ps(Q1,Q1,_mm_shuffle(3,1,2,2));
		V1 = _mm_and_ps(V1,g_XMMask3);
		__m128  R0 = _mm_sub_ps(Constant1110, V0);
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

	__forceinline static Matrix4 VECTORCALL Transpose(const Matrix4 M)
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

	__forceinline static Vector3f VECTORCALL Vector3Transform(const Vector3f V, const Matrix4& M)
	{
		Vector3f result;
		SSEVector3Store(&result.x, Vector3Transform(_mm_loadu_ps(&V.x), M));
		return result;
	}

	__forceinline static __m128 VECTORCALL Vector3Transform(__m128 vec, const Matrix4& M)
	{
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
};
 
__forceinline void InitializeMatrix4(Matrix4& mat, float s) 
{
	mat.r[0] = mat.r[1] = mat.r[2] = mat.r[3] = _mm_set_ps1(s);
}

__forceinline void VECTORCALL InitializeMatrix4(Matrix4& r, __m128 x, __m128 y, const __m128& z, const __m128& w)
{
	r.r[0] = x; r.r[1] = y; r.r[2] = z; r.r[3] = w;
}

#else // sse is not supported

struct Matrix4
{
	struct
	{
		float    m[4][4];
		Vector4f r[4];
	};

	const Vector4f& operator [] (int index) const { return r[index]; }
	Vector4f& operator [] (int index) { return r[index]; }
	
	Vector3f operator *  (const Vector3f v)  noexcept { return Vector3Transform(v, *this); };
	Vector4f operator *  (const Vector4f& v) noexcept { return Vector4Transform(v, *this); };
    Matrix4  operator *  (const Matrix4& M)  noexcept { return Matrix4::Multiply(*this, M); };
	Matrix4& operator *= (const Matrix4& M)  noexcept { *this = Matrix4::Multiply(*this, M); return *this; };

	__forceinline static Matrix4 Make(float x, float y, float z, float w,
	                                  float a, float b, float c, float d,
	                                  float u, float v, float s, float t,
	                                  float j, float k, float l, float m)
	{
		Matrix4 M;
		M.m[0][0] = x; M.m[0][1] = y; M.m[0][2] = z; M.m[0][3] = w; 
		M.m[1][0] = a; M.m[1][1] = b; M.m[1][2] = c; M.m[1][3] = d; 
		M.m[2][0] = u; M.m[2][1] = v; M.m[2][2] = s; M.m[2][3] = t; 
		M.m[3][0] = j; M.m[3][1] = k; M.m[3][2] = l; M.m[3][3] = m; 
		return M;
	}

	__forceinline static Matrix4 Identity()
	{
		Matrix4 M;
		SmallMemSet(&M, 0, sizeof(Matrix4));
		M.m[0][0] = M.m[1][1] = M.m[2][2] = M.m[3][3] = 1.0f; 
		return M;
	}

	__forceinline static Matrix4 FromPosition(float x, float y, float z)
	{
		Matrix4 M;
		SmallMemSet(&M, 0, sizeof(Matrix4));
		M.m[3][0] = x; M.m[3][1] = y; M.m[3][2] = z; M.m[3][3] = 1.0f; 
		return M;
	}

	__forceinline static Matrix4 FromPosition(const Vector3f& vec3)
	{
		return FromPosition(vec3.x, vec3.y, vec3.z);
	}

	__forceinline static Matrix4 CreateScale(const float x, const float y, const float z)
	{
		Matrix4 M;
		M.m[0][0] =    x; M.m[0][1] = 0.0f; M.m[0][2] = 0.0f; M.m[0][3] = 0.0f;
		M.m[1][0] = 0.0f; M.m[1][1] =    y; M.m[1][2] = 0.0f; M.m[1][3] = 0.0f;
		M.m[2][0] = 0.0f; M.m[2][1] = 0.0f; M.m[2][2] =    z; M.m[2][3] = 0.0f;
		M.m[3][0] = 0.0f; M.m[3][1] = 0.0f; M.m[3][2] = 0.0f; M.m[3][3] = 1.0f;
		return M;
	}

	__forceinline static Matrix4 CreateScale(const Vector3f& vec3)
	{
		return CreateScale(vec3.x, vec3.y, vec3.z);
	}

	__forceinline static Matrix4 CreateRotation(Vector3f right, Vector3f up, Vector3f forward)
	{
		Matrix4 m;
		m.r[0][0] = right.x  ; m.r[0][1] = right.y  ; m.r[0][2] = right.z  ; m.r[0][3] = 0.0f;
		m.r[1][0] = up.x     ; m.r[1][1] = up.y     ; m.r[1][2] = up.z     ; m.r[1][3] = 0.0f;
		m.r[2][0] = forward.x; m.r[2][1] = forward.y; m.r[2][2] = forward.z; m.r[2][3] = 0.0f;
		m.r[3][0] = m.r[3][1] = m.r[3][2] = 0.0f; m.r[3][3] = 1.0f;
		return m;
	}

	// please assign normalized vectors, returns view matrix
	__forceinline static Matrix4 LookAtRH(Vector3f eye, Vector3f center, const Vector3f& up)
	{
		Vector3f eyeDir = -eye;
		Vector3f r0 = Vector3f::Normalize(Vector3f::Cross(up, eyeDir));
		Vector3f r1 = Vector3f::Normalize(Vector3f::Cross(eye, r0));
		
		Vector3f negEyePosition = -eye;
		
		float d0 = Vector3f::Dot(r0, negEyePosition);
		float d1 = Vector3f::Dot(r1, negEyePosition);
		float d2 = Vector3f::Dot(eyeDir, negEyePosition);
		
		Matrix4 M;
		M.m[0][0] = r0.x; M.m[0][1] = r0.y; M.m[0][2] = r0.z; M.m[0][3] = d0;
		M.m[1][0] = r1.x; M.m[2][1] = r1.y; M.m[2][2] = r1.z; M.m[2][3] = d1;
		M.m[2][0] = eyeDir.x; M.m[2][1] = eyeDir.y; M.m[2][2] = eyeDir.z; M.m[2][3] = d2;
		M.m[3][0] = 0.0f; M.m[3][1] = 0.0f; M.m[3][2] = 0.0f; M.m[3][3] = 1.0f;
		return Matrix4::Transpose(M);
	}

	__forceinline static Matrix4 PerspectiveFovRH(float fov, float width, float height, float zNear, float zFar)
	{
		const float rad = fov;
		const float h = Cos(0.5f * rad) / Sin(0.5f * rad);
		const float w = h * height / width; /// max(width , Height) / min(width , Height)?
		Matrix4 M = Identity();
		M.m[0][0] = w;
		M.m[1][1] = h;
		M.m[2][2] = -(zFar + zNear) / (zFar - zNear);
		M.m[2][3] = -1.0f;
		M.m[3][2] = -(2.0f * zFar * zNear) / (zFar - zNear);
		M.m[3][3] = 0.0f;
		return M;
	}

	// this will not work on camera matrix this is for only transformation matricies
	inline Matrix4 static InverseTransform(Matrix4 inM)
	{
		Matrix4 out;
		// transpose 3x3, we know m03 = m13 = m23 = 0
		Vector4f t0 = VecShuffle_0101(inM.r[0], inM.r[1]); // 00, 01, 10, 11
		Vector4f t1 = VecShuffle_2323(inM.r[0], inM.r[1]); // 02, 03, 12, 13
		out.r[0]  = VecShuffle(t0, inM.r[2], 0,2,0,3); // 00, 10, 20, 23(=0)
		out.r[1]  = VecShuffle(t0, inM.r[2], 1,3,1,3); // 01, 11, 21, 23(=0)
		out.r[2]  = VecShuffle(t1, inM.r[2], 0,2,2,3); // 02, 12, 22, 23(=0)

		// (SizeSqr(mVec[0]), SizeSqr(mVec[1]), SizeSqr(mVec[2]), 0)
		Vector4f sizeSqr;
		sizeSqr =            out.r[0] * out.r[0];
		sizeSqr = sizeSqr + (out.r[1] * out.r[1]);
		sizeSqr = sizeSqr + (out.r[2] * out.r[2]);

		// optional test to avoid divide by 0
		static Vector4f const one = MakeVec4(1.0f);
		// for each component, if(sizeSqr < SMALL_NUMBER) sizeSqr = 1;
		Vector4f rSizeSqr = one / sizeSqr;
		rSizeSqr.x = rSizeSqr.x < 1.e-8f ? 1.0f : rSizeSqr.x;
		rSizeSqr.y = rSizeSqr.y < 1.e-8f ? 1.0f : rSizeSqr.y;
		rSizeSqr.z = rSizeSqr.z < 1.e-8f ? 1.0f : rSizeSqr.z;
		rSizeSqr.w = rSizeSqr.w < 1.e-8f ? 1.0f : rSizeSqr.w;

		out.r[0] = out.r[0] * rSizeSqr;
		out.r[1] = out.r[1] * rSizeSqr;
		out.r[2] = out.r[2] * rSizeSqr;
		// last line
		out.r[3] =            (out.r[0] * inM.r[3][0]);
		out.r[3] = out.r[3] + (out.r[1] * inM.r[3][1]);
		out.r[3] = out.r[3] + (out.r[2] * inM.r[3][2]);
		out.r[3] = MakeVec4(0.f, 0.f, 0.f, 1.f) - out.r[3];
		return out;
	}

	inline Matrix4 static Inverse(const Matrix4 matrix)
	{
		float v0 = matrix.m[2][0] * matrix.m[3][1] - matrix.m[2][1] * matrix.m[3][0];
		float v1 = matrix.m[2][0] * matrix.m[3][2] - matrix.m[2][2] * matrix.m[3][0];
		float v2 = matrix.m[2][0] * matrix.m[3][3] - matrix.m[2][3] * matrix.m[3][0];
		float v3 = matrix.m[2][1] * matrix.m[3][2] - matrix.m[2][2] * matrix.m[3][1];
		float v4 = matrix.m[2][1] * matrix.m[3][3] - matrix.m[2][3] * matrix.m[3][1];
		float v5 = matrix.m[2][2] * matrix.m[3][3] - matrix.m[2][3] * matrix.m[3][2];

		float i00 =  (v5 * matrix.m[1][1] - v4 * matrix.m[1][2] + v3 * matrix.m[1][3]);
		float i10 = -(v5 * matrix.m[1][0] - v2 * matrix.m[1][2] + v1 * matrix.m[1][3]);
		float i20 =  (v4 * matrix.m[1][0] - v2 * matrix.m[1][1] + v0 * matrix.m[1][3]);
		float i30 = -(v3 * matrix.m[1][0] - v1 * matrix.m[1][1] + v0 * matrix.m[1][2]);

		const float invDet = 1.0f / (i00 * matrix.m[0][0] + i10 * matrix.m[0][1] + i20 * matrix.m[0][2] + i30 * matrix.m[0][3]);

		i00 *= invDet; i10 *= invDet; i20 *= invDet; i30 *= invDet;

		const float i01 = -(v5 * matrix.m[0][1] - v4 * matrix.m[0][2] + v3 * matrix.m[0][3]) * invDet;
		const float i11 =  (v5 * matrix.m[0][0] - v2 * matrix.m[0][2] + v1 * matrix.m[0][3]) * invDet;
		const float i21 = -(v4 * matrix.m[0][0] - v2 * matrix.m[0][1] + v0 * matrix.m[0][3]) * invDet;
		const float i31 =  (v3 * matrix.m[0][0] - v1 * matrix.m[0][1] + v0 * matrix.m[0][2]) * invDet;

		v0 = matrix.m[1][0] * matrix.m[3][1] - matrix.m[1][1] * matrix.m[3][0];
		v1 = matrix.m[1][0] * matrix.m[3][2] - matrix.m[1][2] * matrix.m[3][0];
		v2 = matrix.m[1][0] * matrix.m[3][3] - matrix.m[1][3] * matrix.m[3][0];
		v3 = matrix.m[1][1] * matrix.m[3][2] - matrix.m[1][2] * matrix.m[3][1];
		v4 = matrix.m[1][1] * matrix.m[3][3] - matrix.m[1][3] * matrix.m[3][1];
		v5 = matrix.m[1][2] * matrix.m[3][3] - matrix.m[1][3] * matrix.m[3][2];

		const float i02 =  (v5 * matrix.m[0][1] - v4 * matrix.m[0][2] + v3 * matrix.m[0][3]) * invDet;
		const float i12 = -(v5 * matrix.m[0][0] - v2 * matrix.m[0][2] + v1 * matrix.m[0][3]) * invDet;
		const float i22 =  (v4 * matrix.m[0][0] - v2 * matrix.m[0][1] + v0 * matrix.m[0][3]) * invDet;
		const float i32 = -(v3 * matrix.m[0][0] - v1 * matrix.m[0][1] + v0 * matrix.m[0][2]) * invDet;

		v0 = matrix.m[2][1] * matrix.m[1][0] - matrix.m[2][0] * matrix.m[1][1];
		v1 = matrix.m[2][2] * matrix.m[1][0] - matrix.m[2][0] * matrix.m[1][2];
		v2 = matrix.m[2][3] * matrix.m[1][0] - matrix.m[2][0] * matrix.m[1][3];
		v3 = matrix.m[2][2] * matrix.m[1][1] - matrix.m[2][1] * matrix.m[1][2];
		v4 = matrix.m[2][3] * matrix.m[1][1] - matrix.m[2][1] * matrix.m[1][3];
		v5 = matrix.m[2][3] * matrix.m[1][2] - matrix.m[2][2] * matrix.m[1][3];

		const float i03 = -(v5 * matrix.m[0][1] - v4 * matrix.m[0][2] + v3 * matrix.m[0][3]) * invDet;
		const float i13 =  (v5 * matrix.m[0][0] - v2 * matrix.m[0][2] + v1 * matrix.m[0][3]) * invDet;
		const float i23 = -(v4 * matrix.m[0][0] - v2 * matrix.m[0][1] + v0 * matrix.m[0][3]) * invDet;
		const float i33 =  (v3 * matrix.m[0][0] - v1 * matrix.m[0][1] + v0 * matrix.m[0][2]) * invDet;

		return Make(
			i00, i01, i02, i03,
			i10, i11, i12, i13,
			i20, i21, i22, i23,
			i30, i31, i32, i33
		);
	}

	inline Matrix4 static Multiply(const Matrix4 in1, const Matrix4& in2)
	{
		Matrix4 out;
		out.r[0][0] = in1[0].x * in2[0].x + in1[0].x * in2[0].y + in1[0].x * in2[0].z + in1[0].x * in2[0].w;
		out.r[0][1] = in1[1].y * in2[0].x + in1[1].y * in2[0].y + in1[1].y * in2[0].z + in1[1].y * in2[0].w;
		out.r[0][2] = in1[2].z * in2[0].x + in1[2].z * in2[0].y + in1[2].z * in2[0].z + in1[2].z * in2[0].w;
		out.r[0][3] = in1[3].w * in2[0].x + in1[3].w * in2[0].y + in1[3].w * in2[0].z + in1[3].w * in2[0].w;
		
		out.r[1][0] = in1[0].x * in2[1].x + in1[0].x * in2[1].y + in1[0].x * in2[1].z + in1[0].x * in2[1].w;
		out.r[1][1] = in1[1].y * in2[1].x + in1[1].y * in2[1].y + in1[1].y * in2[1].z + in1[1].y * in2[1].w;
		out.r[1][2] = in1[2].z * in2[1].x + in1[2].z * in2[1].y + in1[2].z * in2[1].z + in1[2].z * in2[1].w;
		out.r[1][3] = in1[3].w * in2[1].x + in1[3].w * in2[1].y + in1[3].w * in2[1].z + in1[3].w * in2[1].w;
		
		out.r[2][0] = in1[0].x * in2[2].x + in1[0].x * in2[2].y + in1[0].x * in2[2].z + in1[0].x * in2[2].w;
		out.r[2][1] = in1[1].y * in2[2].x + in1[1].y * in2[2].y + in1[1].y * in2[2].z + in1[1].y * in2[2].w;
		out.r[2][2] = in1[2].z * in2[2].x + in1[2].z * in2[2].y + in1[2].z * in2[2].z + in1[2].z * in2[2].w;
		out.r[2][3] = in1[3].w * in2[2].x + in1[3].w * in2[2].y + in1[3].w * in2[2].z + in1[3].w * in2[2].w;
		
		out.r[3][0] = in1[0].x * in2[3].x + in1[0].x * in2[3].y + in1[0].x * in2[3].z + in1[0].x * in2[3].w;
		out.r[3][1] = in1[1].y * in2[3].x + in1[1].y * in2[3].y + in1[1].y * in2[3].z + in1[1].y * in2[3].w;
		out.r[3][2] = in1[2].z * in2[3].x + in1[2].z * in2[3].y + in1[2].z * in2[3].z + in1[2].z * in2[3].w;
		out.r[3][3] = in1[3].w * in2[3].x + in1[3].w * in2[3].y + in1[3].w * in2[3].z + in1[3].w * in2[3].w;
		return out;
	}

	__forceinline static Matrix4 PositionRotationScale(const Vector3f& position, const Quaternion& rotation, const Vector3f& scale)
	{
		Matrix4 result = Matrix4::Identity();
		result *= FromPosition(position);
		result *= FromQuaternion(rotation);
		result *= CreateScale(position);
		return result;
	}

	__forceinline static Vector3f ExtractPosition(const Matrix4 matrix) noexcept
	{
		return { matrix.m[3][0], matrix.m[3][1], matrix.m[3][2]};
	}

	__forceinline static Vector3f ExtractScale(const Matrix4 matrix) noexcept
	{
		return MakeVec3(Vector3f::Length(matrix.r[0].xyz()),
		                Vector3f::Length(matrix.r[2].xyz()),
		                Vector3f::Length(matrix.r[1].xyz()));
	}

	__forceinline static Matrix4 RotationX(float angleRadians) {
		Matrix4 out_matrix = Identity();
		float s, c;
		SinCos(angleRadians, &s, &c);
		out_matrix.m[1][1] = c; out_matrix.m[1][2] = s; out_matrix.m[2][1] = -s; out_matrix.m[2][2] = c;
		return out_matrix;
	}

	__forceinline static Matrix4 RotationY(float angleRadians) {
		Matrix4 out_matrix = Identity();
		float s, c;
		SinCos(angleRadians, &s, &c);
		out_matrix.m[0][0] = c; out_matrix.m[0][2] = -s; out_matrix.m[2][0] = s; out_matrix.m[2][2] = c;
		return out_matrix;
	}
	
	__forceinline static Matrix4 RotationZ(float angleRadians) {
		Matrix4 out_matrix = Identity();
		float s, c;
		SinCos(angleRadians, &s, &c);
		out_matrix.m[0][0] = c; out_matrix.m[0][1] = s; out_matrix.m[1][0] = -s; out_matrix.m[1][1] = c;
		return out_matrix;
	}

	__forceinline static Matrix4 RotationFromEuler(Vector3f eulerRadians) 
	{
		return FromQuaternion(Quaternion::FromEuler(eulerRadians));
	}

	__forceinline static Matrix3 ConvertToMatrix3(const Matrix4 M)
	{
		Matrix3 result;
		result.m[0][0] = M.m[0][0]; result.m[0][1] = M.m[0][1]; result.m[0][2] = M.m[0][2];
		result.m[1][0] = M.m[1][0]; result.m[1][1] = M.m[1][1]; result.m[1][2] = M.m[1][2];
		result.m[2][0] = M.m[2][0]; result.m[2][1] = M.m[2][1]; result.m[2][2] = M.m[2][2];
		return result;
	}

	static Quaternion ExtractRotation(const Matrix4 M) noexcept
	{
		Quaternion res;
		QuaternionFromMatrix(&res.x, M.r);
		return res;
	}

	static Matrix4 FromQuaternion(Quaternion q)
	{
		const float num9 = q.x * q.x, num8 = q.y * q.y, num7 = q.z * q.z,
		            num6 = q.x * q.y, num5 = q.z * q.w, num4 = q.z * q.x,
		            num3 = q.y * q.w, num2 = q.y * q.z, num  = q.x * q.w;
		return Make(
			1.0f - (2.0f * (num8 + num7)), 2.0f * (num6 + num5), 2.0f * (num4 - num3), 0.0f,
			2.0f * (num6 - num5), 1.0f - (2.0f * (num7 + num9)), 2.0f * (num2 + num) , 0.0f,
			2.0f * (num4 + num3), 2.0f * (num2 - num), 1.0f - (2.0f * (num8 + num9)) , 0.0f,
			0.0f                , 0.0f               , 0.0f                          , 1.0f
		);
	}

	__forceinline static Matrix4 Transpose(Matrix4 M)
	{
		Matrix4 P;
		P.m[0][0] = M.m[0][0]; P.m[0][1] = M.m[1][0]; P.m[0][2] = M.m[2][0]; P.m[0][3] = M.m[3][0];
		P.m[1][0] = M.m[0][1]; P.m[1][1] = M.m[1][1]; P.m[1][2] = M.m[2][1]; P.m[1][3] = M.m[3][1];
		P.m[2][0] = M.m[0][2]; P.m[2][1] = M.m[1][2]; P.m[2][2] = M.m[2][2]; P.m[2][3] = M.m[3][2];
		P.m[3][0] = M.m[0][3]; P.m[3][1] = M.m[1][3]; P.m[3][2] = M.m[2][3]; P.m[3][3] = M.m[3][3];
		return M;
	}

	__forceinline static Vector3f Vector3Transform(const Vector3f v, const Matrix4& m) noexcept
	{
		return MakeVec3(
			v.x * m.m[0][0] + v.y * m.m[1][0] + v.z * m.m[2][0],
			v.x * m.m[0][1] + v.y * m.m[1][1] + v.z * m.m[2][1],
			v.x * m.m[0][2] + v.y * m.m[1][2] + v.z * m.m[2][2]
		);
	}

	__forceinline static Vector4f Vector4Transform(Vector4f v, const Matrix4& m)
	{
		return MakeVec4(
			v.x * m.m[0][0] + v.y * m.m[1][0] + v.z * m.m[2][0] + v.w * m.m[3][0],
			v.x * m.m[0][1] + v.y * m.m[1][1] + v.z * m.m[2][1] + v.w * m.m[3][1],
			v.x * m.m[0][2] + v.y * m.m[1][2] + v.z * m.m[2][2] + v.w * m.m[3][2],
			v.x * m.m[0][3] + v.y * m.m[1][3] + v.z * m.m[2][3] + v.w * m.m[3][3]
		);
	}
};
 
__forceinline void InitializeMatrix4(Matrix4& mat, float s) 
{
	for (int i = 0; i < 16; i++)
		mat.m[0][i] = s;
}

__forceinline void InitializeMatrix4(Matrix4& mat, Vector4f a, Vector4f b, Vector4f c, Vector4f d) 
{
	mat.r[0] = a;
	mat.r[1] = b;
	mat.r[2] = c;
	mat.r[3] = d;
}

#endif // AX_SUPPORT_SSE

AX_END_NAMESPACE 