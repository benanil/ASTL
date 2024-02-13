#pragma once
#include "Quaternion.hpp"
#include <math.h> // < sinf, cosf for more preciseness

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
	
	float* GetPtr()              { return &m[0][0]; }
	const float* GetPtr() const  { return &m[0][0]; }

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

	__forceinline static Matrix3 TBN(Vector3f normal, Vector3f tangent, Vector3f bitangent)
	{
		Matrix3 M;
		M.vec[0] = normal;
		M.vec[1] = tangent;
		M.vec[2] = bitangent;
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
		result.x = vx + vy + vz;

		vx = b.x * a.y.x;
		vy = b.y * a.y.y;
		vz = b.z * a.y.z;
		result.y = vx + vy + vz;

		vx = b.x * a.z.x;
		vy = b.y * a.z.y;
		vz = b.z * a.z.z;
		result.z = vx + vy + vz;
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
		QuaternionFromMatrix<3>(&Orientation.x, &m[0][0]);
		return Orientation;
	}
};

__forceinline static vec_t VECTORCALL Vector4Transform(vec_t v, const vec_t r[4])
{
	return VecAdd(
		VecAdd(VecMul(r[0], VecSplatX(v)), VecMul(r[1], VecSplatY(v))),
		VecAdd(VecMul(r[2], VecSplatZ(v)), VecMul(r[3], VecSplatW(v)))
	);
}

struct alignas(16) Matrix4
{
	union
	{
		struct { vec_t r[4]; };
		struct { float m[4][4]; };
	};

	const vec_t& operator [] (int index) const { return r[index]; }
	vec_t& operator [] (int index) { return r[index]; }

	Vector3f VECTORCALL  operator * (const Vector3f v)  noexcept { return Vector3Transform(v, *this); };
	vec_t VECTORCALL  operator * (const vec_t& v) noexcept { vec_t x; x = ::Vector4Transform(v, r); return x; };

	Matrix4 VECTORCALL  operator *  (const Matrix4& M)  noexcept { return Matrix4::Multiply(M, *this); };
	Matrix4& VECTORCALL operator *= (const Matrix4& M)  noexcept { *this = Matrix4::Multiply(M, *this); return *this; };

	float* GetPtr()              { return &m[0][0]; }
	const float* GetPtr() const  { return &m[0][0]; }

	__forceinline static Matrix4 Make(float x, float y, float z, float w,
	                                  float a, float b, float c, float d,
	                                  float u, float v, float s, float t,
	                                  float j, float k, float l, float m)
	{
		Matrix4 M;
		M.r[0] = VecSetR(x, y, z, w);
		M.r[1] = VecSetR(a, b, c, d);
		M.r[2] = VecSetR(u, v, s, t);
		M.r[3] = VecSetR(j, k, l, m);
		return M;
	}

	__forceinline static Matrix4 Identity()
	{
		Matrix4 M;
		M.r[0] = VecIdentityR0;
		M.r[1] = VecIdentityR1;
		M.r[2] = VecIdentityR2;
		M.r[3] = VecIdentityR3;
		return M;
	}

	__forceinline static Matrix4 FromPosition(const float x, const float y, const float z)
	{
		Matrix4 M;
		M.r[0] = VecIdentityR0;
		M.r[1] = VecIdentityR1;
		M.r[2] = VecIdentityR2;
		M.r[3] = VecSet(1.0f, z, y, x);
		return M;
	}

	__forceinline static Matrix4 FromPosition(const float* vec3)
	{
		return FromPosition(vec3[0], vec3[1], vec3[2]);
	}

	__forceinline static Matrix4 FromPosition(const Vector3f& vec3)
	{
		return FromPosition(vec3.x, vec3.y, vec3.z);
	}

	__forceinline static Matrix4 CreateScale(const float ScaleX, const float ScaleY, const float ScaleZ)
	{
		Matrix4 M;
		M.r[0] = VecSet(0.0f, 0.0f, 0.0f, ScaleX);
		M.r[1] = VecSet(0.0f, 0.0f, ScaleY, 0.0f);
		M.r[2] = VecSet(0.0f, ScaleZ, 0.0f, 0.0f);
		M.r[3] = VecIdentityR3;
		return M;
	}

	__forceinline static Matrix4 CreateScale(Vector3f vec3)
	{
		return CreateScale(vec3.x, vec3.y, vec3.z);
	}
	
	__forceinline static Matrix4 CreateScale(float* vec3)
	{
		return CreateScale(vec3[0], vec3[1], vec3[2]);
	}

	__forceinline static Matrix4 CreateScale(float scale)
	{
		return CreateScale(scale, scale, scale);
	}

	__forceinline static Matrix4 CreateRotation(Vector3f right, Vector3f up, Vector3f forward)
	{
		Matrix4 m;
		veci_t mask = VecMask3;
		m.r[0] = VecMask(VecLoad(&right.x)  , mask);
		m.r[1] = VecMask(VecLoad(&up.x)     , mask);
		m.r[2] = VecMask(VecLoad(&forward.x), mask);
		m.r[3] = VecSetR(0.0f, 0.0f, 0.0f, 1.0f);
		return m;
	}

	const Vector3f GetForward() const { Vector3f res; Vec3Store(&res.x, r[2]); return res; }
	const Vector3f GetUp()      const { Vector3f res; Vec3Store(&res.x, r[1]); return res; }
	const Vector3f GetRight()   const { Vector3f res; Vec3Store(&res.x, r[0]); return res; }

	// please assign normalized vectors, returns view matrix
	// creates view matrix
	__forceinline static Matrix4 VECTORCALL LookAtRH(Vector3f eye, Vector3f center, const Vector3f& up)
	{
		vec_t EyePosition  = VecLoad(&eye.x);
		vec_t EyeDirection = VecSub(VecZero(), VecLoad(&center.x));
		vec_t UpDirection  = VecLoad(&up.x);

		vec_t R0 = Vec3Norm(Vec3Cross(UpDirection, EyeDirection));
		vec_t R1 = Vec3Norm(Vec3Cross(EyeDirection, R0));

		vec_t NegEyePosition = VecSub(VecZero(), EyePosition);

		vec_t D0 = Vec3Dot(R0, NegEyePosition);
		vec_t D1 = Vec3Dot(R1, NegEyePosition);
		vec_t D2 = Vec3Dot(EyeDirection, NegEyePosition);
		Matrix4 M;
		veci_t test = VecSelect1110;
		M.r[0] = VecSelect(D0, R0, test); // no need select ?
		M.r[1] = VecSelect(D1, R1, test);
		M.r[2] = VecSelect(D2, EyeDirection, test);
		M.r[3] = VecIdentityR3;
		return Matrix4::Transpose(M);
	}

	__forceinline static Matrix4 PerspectiveFovRH(float fov, float width, float height, float zNear, float zFar)
	{
		const float rad = fov;
		const float h = cosf(0.5f * rad) / sinf(0.5f * rad);
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

	__forceinline static Matrix4 OrthoRH(float left, float right, float bottom, float top, float zNear, float zFar)
	{
		Matrix4 Result = Identity();
		Result.m[0][0] =  2.0f / (right - left);
		Result.m[1][1] =  2.0f / (top - bottom);
		Result.m[2][2] = -2.0f / (zFar - zNear);
		Result.m[3][0] = -(right + left) / (right - left);
		Result.m[3][1] = -(top + bottom) / (top - bottom);
		Result.m[3][2] = -(zFar + zNear) / (zFar - zNear);
		return Result;
	}

	// https://lxjk.github.io/2017/09/03/Fast-4x4-Matrix-Inverse-with-SSE-SIMD-Explained.html
	// for row major matrix
    // we use vec_t to represent 2x2 matrix as A = | A0  A1 |
    //                                             | A2  A3 |
    // 2x2 row major Matrix multiply A*B
	static __forceinline vec_t VECTORCALL Mat2Mul(vec_t vec1, vec_t vec2)
	{
		return VecAdd(VecMul(vec1, VecSwizzle(vec2, 0, 3, 0, 3)), 
		              VecMul(VecSwizzle(vec1, 1, 0, 3, 2), VecSwizzle(vec2, 2, 1, 2, 1)));
	}
	// 2x2 row major Matrix adjugate multiply (A#)*B
	static __forceinline vec_t VECTORCALL Mat2AdjMul(vec_t vec1, vec_t vec2)
	{
		return VecSub(VecMul(VecSwizzle(vec1, 3, 3, 0, 0), vec2),
		              VecMul(VecSwizzle(vec1, 1, 1, 2, 2), VecSwizzle(vec2, 2, 3, 0, 1)));
	}
	// 2x2 row major Matrix multiply adjugate A*(B#)
	static __forceinline vec_t VECTORCALL Mat2MulAdj(vec_t vec1, vec_t vec2)
	{
		return VecSub(VecMul(vec1, VecSwizzle(vec2, 3, 0, 3, 0)),
		              VecMul(VecSwizzle(vec1, 1, 0, 3, 2), VecSwizzle(vec2, 2, 1, 2, 1)));
	}

	// this will not work on camera matrix this is for only transformation matricies
	inline Matrix4 static VECTORCALL InverseTransform(const Matrix4 inM) noexcept
	{
		Matrix4 out;
		// transpose 3x3, we know m03 = m13 = m23 = 0
		vec_t t0 = VecShuffle_0101(inM.r[0], inM.r[1]); // 00, 01, 10, 11
		vec_t t1 = VecShuffle_2323(inM.r[0], inM.r[1]); // 02, 03, 12, 13
		out.r[0] = VecShuffle(t0, inM.r[2], 0, 2, 0, 3); // 00, 10, 20, 23(=0)
		out.r[1] = VecShuffle(t0, inM.r[2], 1, 3, 1, 3); // 01, 11, 21, 23(=0)
		out.r[2] = VecShuffle(t1, inM.r[2], 0, 2, 2, 3); // 02, 12, 22, 23(=0)

		// (SizeSqr(mVec[0]), SizeSqr(mVec[1]), SizeSqr(mVec[2]), 0)
		vec_t sizeSqr;
		sizeSqr =                 VecMul(out.r[0], out.r[0]);
		sizeSqr = VecAdd(sizeSqr, VecMul(out.r[1], out.r[1]));
		sizeSqr = VecAdd(sizeSqr, VecMul(out.r[2], out.r[2]));

		// optional test to avoid divide by 0
		const vec_t one = VecOne();
		// for each component, if(sizeSqr < SMALL_NUMBER) sizeSqr = 1;
		vec_t rSizeSqr = VecSelect(
			VecDiv(one, sizeSqr), one,
			VecCmpLt(sizeSqr, VecSet1(1.e-8f))
		);

		out.r[0] = VecMul(out.r[0], rSizeSqr);
		out.r[1] = VecMul(out.r[1], rSizeSqr);
		out.r[2] = VecMul(out.r[2], rSizeSqr);
		// last line
		out.r[3] =                  VecMul(out.r[0], VecSplatX(inM.r[3]));
		out.r[3] = VecAdd(out.r[3], VecMul(out.r[1], VecSplatY(inM.r[3])));
		out.r[3] = VecAdd(out.r[3], VecMul(out.r[2], VecSplatZ(inM.r[3])));
		out.r[3] = VecSub(VecSetR(0.f, 0.f, 0.f, 1.f), out.r[3]);
		return out;
	}

	inline Matrix4 static VECTORCALL Inverse(const Matrix4 matrix) noexcept
	{
		#if defined(AX_ARM)
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
		#else
		vec_t A = VecShuffle_0101(matrix.r[0], matrix.r[1]);
		vec_t B = VecShuffle_2323(matrix.r[0], matrix.r[1]);
		vec_t C = VecShuffle_0101(matrix.r[2], matrix.r[3]);
		vec_t D = VecShuffle_2323(matrix.r[2], matrix.r[3]);

		vec_t detSub = VecSub(
			VecMul(VecShuffle(matrix.r[0], matrix.r[2], 0, 2, 0, 2), VecShuffle(matrix.r[1], matrix.r[3], 1, 3, 1, 3)),
			VecMul(VecShuffle(matrix.r[0], matrix.r[2], 1, 3, 1, 3), VecShuffle(matrix.r[1], matrix.r[3], 0, 2, 0, 2))
		);
		vec_t detA = VecSplatX(detSub);
		vec_t detB = VecSplatY(detSub);
		vec_t detC = VecSplatZ(detSub);
		vec_t detD = VecSplatW(detSub);

		vec_t D_C  = Mat2AdjMul(D, C);
		vec_t A_B  = Mat2AdjMul(A, B);
		vec_t X_   = VecSub(VecMul(detD, A), Mat2Mul(B, D_C));
		vec_t W_   = VecSub(VecMul(detA, D), Mat2Mul(C, A_B));
		
		vec_t detM = VecMul(detA, detD);
		vec_t Y_   = VecSub(VecMul(detB, C), Mat2MulAdj(D, A_B));
		vec_t Z_   = VecSub(VecMul(detC, B), Mat2MulAdj(A, D_C));

		detM = VecAdd(detM, VecMul(detB, detC));

		vec_t tr = VecMul(A_B, VecSwizzle(D_C, 0, 2, 1, 3));
		tr   = VecHadd(tr, tr);
		tr   = VecHadd(tr, tr);
		detM = VecSub(detM, tr);

		const vec_t adjSignMask = VecSetR(1.f, -1.f, -1.f, 1.f);
		vec_t rDetM = VecDiv(adjSignMask, detM);
		X_ = VecMul(X_, rDetM);
		Y_ = VecMul(Y_, rDetM);
		Z_ = VecMul(Z_, rDetM);
		W_ = VecMul(W_, rDetM);

		Matrix4 out;
		out.r[0] = VecShuffle(X_, Y_, 3, 1, 3, 1);
		out.r[1] = VecShuffle(X_, Y_, 2, 0, 2, 0);
		out.r[2] = VecShuffle(Z_, W_, 3, 1, 3, 1);
		out.r[3] = VecShuffle(Z_, W_, 2, 0, 2, 0);
		return out;
		#endif
	}

	inline Matrix4 static VECTORCALL Multiply(const Matrix4 in1, const Matrix4& in2)
	{
		Matrix4 out;
		vec_t m0 = VecMul(in1.r[0], VecSplatX(in2.r[0]));
		vec_t m1 = VecMul(in1.r[1], VecSplatY(in2.r[0]));
		vec_t m2 = VecMul(in1.r[2], VecSplatZ(in2.r[0]));
		vec_t m3 = VecMul(in1.r[3], VecSplatW(in2.r[0]));
		out.r[0] = VecAdd(VecAdd(m0, m1), VecAdd(m2, m3));

		m0 = VecMul(in1.r[0], VecSplatX(in2.r[1]));
		m1 = VecMul(in1.r[1], VecSplatY(in2.r[1]));
		m2 = VecMul(in1.r[2], VecSplatZ(in2.r[1]));
		m3 = VecMul(in1.r[3], VecSplatW(in2.r[1]));
		out.r[1] = VecAdd(VecAdd(m0, m1), VecAdd(m2, m3));

		m0 = VecMul(in1.r[0], VecSplatX(in2.r[2]));
		m1 = VecMul(in1.r[1], VecSplatY(in2.r[2]));
		m2 = VecMul(in1.r[2], VecSplatZ(in2.r[2]));
		m3 = VecMul(in1.r[3], VecSplatW(in2.r[2]));
		out.r[2] = VecAdd(VecAdd(m0, m1), VecAdd(m2, m3));

		m0 = VecMul(in1.r[0], VecSplatX(in2.r[3]));
		m1 = VecMul(in1.r[1], VecSplatY(in2.r[3]));
		m2 = VecMul(in1.r[2], VecSplatZ(in2.r[3]));
		m3 = VecMul(in1.r[3], VecSplatW(in2.r[3]));
		out.r[3] = VecAdd(VecAdd(m0, m1), VecAdd(m2, m3));
		return out;
	}

	__forceinline static Matrix4 PositionRotationScale(Vector3f position, Quaternion rotation, const Vector3f& scale)
	{
		return CreateScale(scale) * FromQuaternion(rotation) * FromPosition(position);
	}
	
	__forceinline static Matrix4 PositionRotationScale(float* position, float* rotation, float* scale)
	{
		return CreateScale(scale) * FromQuaternion(rotation) * FromPosition(position);
	}

	__forceinline static Vector3f VECTORCALL ExtractPosition(const Matrix4 matrix) noexcept
	{
		Vector3f res;
		res.x = matrix.m[3][0];
		res.y = matrix.m[3][1];
		res.z = matrix.m[3][2];
		return res;
	}

	static Quaternion VECTORCALL ExtractRotation(const Matrix4 M, bool rowNormalize = true) 
	{
		Quaternion res;
		QuaternionFromMatrix(&res.x, &M.m[0][0]);
		return res;
	}

	__forceinline static Vector3f VECTORCALL ExtractScale(const Matrix4 matrix) noexcept
	{
		return MakeVec3(Vec3Lenf(matrix.r[0]), Vec3Lenf(matrix.r[2]), Vec3Lenf(matrix.r[1]));
	}

	void SetPosition(Vector3f position) noexcept
	{
		m[3][0] = position.x;
		m[3][1] = position.y;
		m[3][2] = position.z;
	}

	Vector3f GetPosition() { return ExtractPosition(*this); }

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
	
	__forceinline static Matrix4 RotationFromEuler(float x, float y, float z)
	{
	  return FromQuaternion(Quaternion::FromEuler(MakeVec3(x, y, z)));
	}

	static Matrix3 VECTORCALL ConvertToMatrix3(const Matrix4 M)
	{
		Matrix3 result;
		Vec3Store(&result.x.x, M.r[0]);
		Vec3Store(&result.y.x, M.r[1]);
		Vec3Store(&result.z.x, M.r[2]);
		return result;
	}

	static Matrix4 VECTORCALL FromQuaternion(const Quaternion quaternion)
	{
		#if defined(AX_ARM)
		Quaternion q = quaternion;
		const float num9 = q.x * q.x, num8 = q.y * q.y, num7 = q.z * q.z,
		            num6 = q.x * q.y, num5 = q.z * q.w, num4 = q.z * q.x,
		            num3 = q.y * q.w, num2 = q.y * q.z, num  = q.x * q.w;
		return Make(
			1.0f - (2.0f * (num8 + num7)), 2.0f * (num6 + num5), 2.0f * (num4 - num3), 0.0f,
			2.0f * (num6 - num5), 1.0f - (2.0f * (num7 + num9)), 2.0f * (num2 + num) , 0.0f,
			2.0f * (num4 + num3), 2.0f * (num2 - num), 1.0f - (2.0f * (num8 + num9)) , 0.0f,
			0.0f                , 0.0f               , 0.0f                          , 1.0f
		);
		#else
		Matrix4 M;
		const vec_t  Constant1110 = VecSetR(1.0f, 1.0f, 1.0f, 0.0f);
		vec_t q = quaternion.vec;

		vec_t Q0 = VecAdd(q, q);
		vec_t Q1 = VecMul(q,Q0);

		vec_t V0 = VecShuffleR(Q1, Q1, 3,0,0,1);
		V0 = VecMask(V0, VecMask3);
		vec_t V1 = VecShuffleR(Q1, Q1, 3,1,2,2);
		V1 = VecMask(V1, VecMask3);
		vec_t  R0 = VecSub(Constant1110, V0);
		R0 = VecSub(R0, V1);

		V0 = VecShuffleR(q, q, 3,1,0,0);
		V1 = VecShuffleR(Q0, Q0, 3,2,1,2);
		V0 = VecMul(V0, V1);

		V1 = VecShuffleR(q, q, 3,3,3,3);
		vec_t V2 = VecShuffleR(Q0, Q0, 3,0,2,1);
		V1 = VecMul(V1, V2);

		vec_t R1 = VecAdd(V0, V1);
		vec_t R2 = VecSub(V0, V1);

		V0 = VecShuffleR(R1, R2, 1, 0, 2, 1);
		V0 = VecShuffleR(V0, V0, 1, 3, 2, 0);
		V1 = VecShuffleR(R1, R2, 2, 2, 0, 0);
		V1 = VecShuffleR(V1, V1, 2, 0, 2, 0);

		Q1 = VecShuffleR(R0, V0, 1, 0, 3, 0);
		Q1 = VecShuffleR(Q1, Q1, 1, 3, 2, 0);
		M.r[0] = Q1;
		Q1 = VecShuffleR(R0, V0, 3, 2, 3, 1);
		Q1 = VecShuffleR(Q1, Q1, 1, 3, 0, 2);
		M.r[1] = Q1;
		Q1 = VecShuffleR(V1, R0, 3, 2, 1, 0);
		M.r[2] = Q1;
		M.r[3] = VecIdentityR3;
		return M;
		#endif
	}

	static Matrix4 VECTORCALL FromQuaternion(const float* quaternion)
	{
		return FromQuaternion({ quaternion[0], quaternion[1], quaternion[2], quaternion[3] });
	}

	__forceinline static Matrix4 VECTORCALL Transpose(const Matrix4 M)
	{
		#ifdef AX_ARM
		float32x4x2_t P0 = vzipq_f32(M.r[0], M.r[2]);
		float32x4x2_t P1 = vzipq_f32(M.r[1], M.r[3]);
		float32x4x2_t T0 = vzipq_f32(P0.val[0], P1.val[0]);
		float32x4x2_t T1 = vzipq_f32(P0.val[1], P1.val[1]);
		Matrix4 mResult;
		mResult.r[0] = T0.val[0];
		mResult.r[1] = T0.val[1];
		mResult.r[2] = T1.val[0];
		mResult.r[3] = T1.val[1];
		#else
		const vec_t vTemp1 = VecShuffleR(M.r[0], M.r[1], 1, 0, 1, 0);
		const vec_t vTemp3 = VecShuffleR(M.r[0], M.r[1], 3, 2, 3, 2);
		const vec_t vTemp2 = VecShuffleR(M.r[2], M.r[3], 1, 0, 1, 0);
		const vec_t vTemp4 = VecShuffleR(M.r[2], M.r[3], 3, 2, 3, 2);
		Matrix4 mResult;
		mResult.r[0] = VecShuffleR(vTemp1, vTemp2, 2, 0, 2, 0);
		mResult.r[1] = VecShuffleR(vTemp1, vTemp2, 3, 1, 3, 1);
		mResult.r[2] = VecShuffleR(vTemp3, vTemp4, 2, 0, 2, 0);
		mResult.r[3] = VecShuffleR(vTemp3, vTemp4, 3, 1, 3, 1);
		#endif
		return mResult;
	}

	__forceinline static vec_t VECTORCALL Vector4Transform(vec_t V, const Matrix4& M)
	{
		return ::Vector4Transform(V, M.r);
	}

	__forceinline static Vector3f VECTORCALL Vector3Transform(const Vector3f V, const Matrix4& M)
	{
		Vector3f result;
		Vec3Store(&result.x, Vector3Transform(VecLoad(&V.x), M));
		return result;
	}

	__forceinline static vec_t VECTORCALL Vector3Transform(vec_t vec, const Matrix4& M)
	{
		vec_t vResult = VecSplatX(vec);
		vResult       = VecMul(vResult, M.r[0]);
		vec_t vTemp   = VecSplatY(vec);
		vTemp   = VecMul(vTemp, M.r[1]);
		vResult = VecAdd(vResult, vTemp);
		vTemp   = VecSplatZ(vec);
		vTemp   = VecMul(vTemp, M.r[2]);
		vResult = VecAdd(vResult, vTemp);
		vResult = VecAdd(vResult, M.r[3]);
		return vResult;
	}
};
 
__forceinline void InitializeMatrix4(Matrix4& mat, float s) 
{
	mat.r[0] = mat.r[1] = mat.r[2] = mat.r[3] = VecSet1(s);
}

__forceinline void VECTORCALL InitializeMatrix4(Matrix4& r, vec_t x, vec_t y, const vec_t& z, const vec_t& w)
{
	r.r[0] = x; r.r[1] = y; r.r[2] = z; r.r[3] = w;
}

AX_END_NAMESPACE 