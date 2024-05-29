
/*********************************************************************************
*    Description:                                                                *
*        Row major right handed vectorized Matrix3x3 and Matrix4x4 structures.   *
*        That is capable of creating and manipulating position, rotation, scale, *
*        view and projection matrices.                                           * 
*    Note:                                                                       *
*        There is frustum culling code at bottom most lines of this file.        *
*    Author:                                                                     *
*        Anilcan Gulkaya 2023 anilcangulkaya7@gmail.com github @benanil          *
*********************************************************************************/

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
    
    float* GetPtr()              { return &m[0][0]; }
    const float* GetPtr() const  { return &m[0][0]; }
    
    static Matrix3 Make(float x, float y, float z,
                        float a, float b, float c,
                        float u, float v, float s)
    {
    	Matrix3 M;
    	M.m[0][0] = x; M.m[0][1] = y; M.m[0][2] = z; 
    	M.m[1][0] = a; M.m[1][1] = b; M.m[1][2] = c; 
    	M.m[2][0] = u; M.m[2][1] = v; M.m[2][2] = s; 
    	return M;
    }
    
    static Matrix3 TBN(Vector3f normal, Vector3f tangent, Vector3f bitangent)
    {
        Matrix3 M;
        M.vec[0] = normal;
        M.vec[1] = tangent;
        M.vec[2] = bitangent;
        return M;
    }
    static Matrix3 Identity()
    {
        return Make(1.0f, 0.0f, 0.0f,
                    0.0f, 1.0f, 0.0f,
                    0.0f, 0.0f, 1.0f);
    }
    
    static Matrix3 LookAt(Vector3f direction, Vector3f up)
    {
        Matrix3 result;
        result.vec[2] = direction;
        Vector3f Right = Vector3f::Cross(up, result.vec[2]);
        result.vec[0] = Right * RSqrt(MAX(0.00001f, Vector3f::Dot(Right, Right)));
        result.vec[1] = Vector3f::Cross(result.vec[2], result.vec[0]);
        return result;
    }
    
    static Matrix3 Multiply(const Matrix3& a, const Matrix3& b)
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
    
    static float3 Multiply(const Matrix3& m, const float3& v) {
        return m.x * v.x + m.y * v.y + m.z * v.z;
    }
    
    static Matrix3 FromQuaternion(const Quaternion quat)
    {
        Matrix3 mat = {};
        MatrixFromQuaternion<3>(mat.GetPtr(), quat);
        return mat;
    }
    
    Quaternion ToQuaternion() const
    {
        Quaternion Orientation;
        QuaternionFromMatrix<3>((float*)&Orientation, &m[0][0]);
        return Orientation;
    }
};

static vec_t VECTORCALL Vector4Transform(vec_t v, const vec_t r[4])
{
    vec_t m0;
    m0 = VecMul(r[0], VecSplatX(v));
    m0 = VecFmaddLane(r[1], v, m0, 1);
    m0 = VecFmaddLane(r[2], v, m0, 2); 
    m0 = VecFmaddLane(r[3], v, m0, 3);
    return m0;
}

static vec_t VECTORCALL Vector3Transform(vec_t vec, const vec_t r[4])
{
    vec_t m0;
    m0 = VecMul(r[0], VecSplatX(vec));
    m0 = VecFmaddLane(r[1], vec, m0, 1);
    m0 = VecFmaddLane(r[2], vec, m0, 2); 
    m0 = VecAdd(r[3], m0);
    return m0;
}

struct alignas(16) Matrix4
{
    union
    {
        struct { vec_t r[4]; };
        struct { float m[4][4]; };
    };
    
    const vec_t& operator [] (int index) const { return r[index]; }
          vec_t& operator [] (int index)       { return r[index]; }
    
    vec_t    VECTORCALL operator *  (vec_t v)          { vec_t x; x = ::Vector4Transform(v, r); return x; };
    Matrix4  VECTORCALL operator *  (const Matrix4 M) { return Matrix4::Multiply(M, *this); };
    Matrix4& VECTORCALL operator *= (const Matrix4 M) { *this = Matrix4::Multiply(M, *this); return *this; };
    
          float* GetPtr()        { return &m[0][0]; }
    const float* GetPtr() const  { return &m[0][0]; }

    void SetPosition(Vector3f position) 
    {
        m[3][0] = position.x;
        m[3][1] = position.y;
        m[3][2] = position.z;
    }

    Vector3f GetPosition() { return ExtractPosition(*this); }

    Vector3f GetForward() const { Vector3f res; Vec3Store(&res.x, r[2]); return res; }
    Vector3f GetUp()      const { Vector3f res; Vec3Store(&res.x, r[1]); return res; }
    Vector3f GetRight()   const { Vector3f res; Vec3Store(&res.x, r[0]); return res; }
    
    static Matrix4 Make(float x, float y, float z, float w,
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
    
    static Matrix4 Identity()
    {
        Matrix4 M;
        M.r[0] = VecIdentityR0;
        M.r[1] = VecIdentityR1;
        M.r[2] = VecIdentityR2;
        M.r[3] = VecIdentityR3;
        return M;
    }
    
    static Matrix4 FromPosition(float x, float y, float z)
    {
        Matrix4 M;
        M.r[0] = VecIdentityR0;
        M.r[1] = VecIdentityR1;
        M.r[2] = VecIdentityR2;
        M.r[3] = VecSet(1.0f, z, y, x);
        return M;
    }
    
    static Matrix4 FromPosition(const float* vec3)
    {
    	return FromPosition(vec3[0], vec3[1], vec3[2]);
    }
    
    static Matrix4 FromPosition(Vector3f vec3)
    {
    	return FromPosition(vec3.x, vec3.y, vec3.z);
    }
    
    static Matrix4 CreateScale(float ScaleX, float ScaleY, float ScaleZ)
    {
        Matrix4 M;
        M.r[0] = VecSetR(ScaleX, 0.0f, 0.0f, 0.0f);
        M.r[1] = VecSetR(0.0f, ScaleY, 0.0f, 0.0f);
        M.r[2] = VecSetR(0.0f, 0.0f, ScaleZ, 0.0f);
        M.r[3] = VecIdentityR3;
        return M;
    }
    
    static Matrix4 CreateScale(Vector3f vec3)
    {
    	return CreateScale(vec3.x, vec3.y, vec3.z);
    }
    
    static Matrix4 CreateScale(float* vec3)
    {
    	return CreateScale(vec3[0], vec3[1], vec3[2]);
    }
    
    static Matrix4 CreateScale(float scale)
    {
    	return CreateScale(scale, scale, scale);
    }
    
    static Matrix4 CreateRotation(Vector3f right, Vector3f up, Vector3f forward)
    {
        Matrix4 m;
        m.r[0] = Vec3Load(&right.x);
        m.r[1] = Vec3Load(&up.x);
        m.r[2] = Vec3Load(&forward.x);
        m.r[3] = VecIdentityR3;
        return m;
    }
        
    // please assign normalized vectors, returns view matrix
    // creates view matrix
    static Matrix4 VECTORCALL LookAtRH(Vector3f eye, Vector3f center, const Vector3f& up)
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
    
    static Matrix4 PerspectiveFovRH(float fov, float width, float height, float zNear, float zFar)
    {
        float rad = Sin0pi(0.5f * fov); 
		AX_ASSUME(rad > 0.01f);
        float h = Sqrt(1.0f - (rad * rad)) / rad;
        float w = h * height / width; /// max(width , Height) / min(width , Height)?
        Matrix4 M = {};
        M.m[0][0] = w;
        M.m[1][1] = h;
        M.m[2][2] = -(zFar + zNear) / (zFar - zNear);
        M.m[2][3] = -1.0f;
        M.m[3][2] = -(2.0f * zFar * zNear) / (zFar - zNear);
        M.m[3][3] = 0.0f;
        return M;
    }
    
    static Matrix4 OrthoRH(float left, float right, float bottom, float top, float zNear, float zFar)
    {
        Matrix4 Result = {};
        Result.m[0][0] =  2.0f / (right - left);
        Result.m[1][1] =  2.0f / (top - bottom);
        Result.m[2][2] = -2.0f / (zFar - zNear);
        Result.m[3][0] = -(right + left) / (right - left);
        Result.m[3][1] = -(top + bottom) / (top - bottom);
        Result.m[3][2] = -(zFar + zNear) / (zFar - zNear);
        Result.m[3][3] = 1.0f;
        return Result;
    }
    
    static Matrix4 PositionRotationScale(Vector3f position, Quaternion rotation, const Vector3f& scale)
    {
        Matrix4 res;
        // Export rotation to matrix
        MatrixFromQuaternion<4>(res.GetPtr(), rotation);
        // Scale 3x3 matrix by given scale
        res.r[0] = VecMulf(res.r[0], scale.x);
        res.r[1] = VecMulf(res.r[1], scale.y);
        res.r[2] = VecMulf(res.r[2], scale.z);
        // Third row is position, x, y, z, 1.0f
        res.r[3] = VecLoad(&position.x);
        VecSetW(res.r[3], 1.0f);
        return res; 
    }
    
    static Matrix4 PositionRotationScale(const float* position, const float* rotation, const float* scale)
    {
        Matrix4 res = {};
        // Export rotation to matrix
        MatrixFromQuaternion<4>(res.GetPtr(), VecLoad(rotation));
        // Scale 3x3 matrix by given scale
        vec_t vecScale = VecLoad(scale);
        res.r[0] = VecMul(res.r[0], VecSplatX(vecScale));
        res.r[1] = VecMul(res.r[1], VecSplatY(vecScale));
        res.r[2] = VecMul(res.r[2], VecSplatZ(vecScale));
        // Third row is position, x, y, z, 1.0f
        res.r[3] = VecLoad(position);
        VecSetW(res.r[3], 1.0f);
        return res; 
    }
    
    static Vector3f VECTORCALL ExtractPosition(Matrix4 matrix)
    {
        Vector3f res;
        Vec3Store(&res.x, matrix.r[3]);
        return res;
    }
    
    static Quaternion VECTORCALL ExtractRotation(const Matrix4 M, bool rowNormalize = true) 
    {
        Quaternion res;
        QuaternionFromMatrix((float*)&res, &M.m[0][0]);
        return res;
    }
    
    static Vector3f VECTORCALL ExtractScale(const Matrix4 matrix) 
    {
    	return MakeVec3(Vec3Lenf(matrix.r[0]), Vec3Lenf(matrix.r[2]), Vec3Lenf(matrix.r[1]));
    }
        
    static vec_t VECTORCALL ExtractScaleV(const Matrix4 matrix) 
    {
    	return VecSetR(Vec3Lenf(matrix.r[0]), Vec3Lenf(matrix.r[2]), Vec3Lenf(matrix.r[1]), 0.0f);
    }

    static Matrix4 RotationX(float angleRadians) {
    	Matrix4 out_matrix = Identity();
    	float s, c;
    	SinCos(angleRadians, &s, &c);
    	out_matrix.m[1][1] = c;
    	out_matrix.m[1][2] = s;
    	out_matrix.m[2][1] = -s;
    	out_matrix.m[2][2] = c;
        // out_matrix.m[3][3] = 1.0f;
    	return out_matrix;
    }
    
    static Matrix4 RotationY(float angleRadians) {
        Matrix4 out_matrix = Identity();
        float s, c;
        SinCos(angleRadians, &s, &c);
        out_matrix.m[0][0] = c;
        out_matrix.m[0][2] = -s;
        out_matrix.m[2][0] = s;
        out_matrix.m[2][2] = c;
        // out_matrix.m[3][3] = 1.0f;
        return out_matrix;
    }
    
    static Matrix4 RotationZ(float angleRadians) {
        Matrix4 out_matrix = Identity();
        float s, c;
        SinCos(angleRadians, &s, &c);
        out_matrix.m[0][0] = c;
        out_matrix.m[0][1] = s;
        out_matrix.m[1][0] = -s;
        out_matrix.m[1][1] = c;
        // out_matrix.m[3][3] = 1.0f;
        return out_matrix;
    }
    
    static Matrix4 RotationFromEuler(Vector3f eulerRadians) 
    {
    	return FromQuaternion(QFromEuler(eulerRadians));
    }
    
    static Matrix4 RotationFromEuler(float x, float y, float z)
    {
      return FromQuaternion(QFromEuler(x, y, z));
    }
    
    static Matrix3 VECTORCALL ConvertToMatrix3(const Matrix4 M)
    {
        Matrix3 result;
        Vec3Store(&result.x.x, M.r[0]);
        Vec3Store(&result.y.x, M.r[1]);
        Vec3Store(&result.z.x, M.r[2]);
        return result;
    }
    // https://lxjk.github.io/2017/09/03/Fast-4x4-Matrix-Inverse-with-SSE-SIMD-Explained.html
    // for row major matrix
    // we use vec_t to represent 2x2 matrix as A = | A0  A1 |
    //                                             | A2  A3 |
    // 2x2 row major Matrix multiply A*B
    static vec_t VECTORCALL Mat2Mul(vec_t vec1, vec_t vec2)
    {
    	return VecAdd(VecMul(vec1, VecSwizzle(vec2, 0, 3, 0, 3)), 
    	              VecMul(VecSwizzle(vec1, 1, 0, 3, 2), VecSwizzle(vec2, 2, 1, 2, 1)));
    }
    // 2x2 row major Matrix adjugate multiply (A#)*B
    static vec_t VECTORCALL Mat2AdjMul(vec_t vec1, vec_t vec2)
    {
    	return VecSub(VecMul(VecSwizzle(vec1, 3, 3, 0, 0), vec2),
    	              VecMul(VecSwizzle(vec1, 1, 1, 2, 2), VecSwizzle(vec2, 2, 3, 0, 1)));
    }
    // 2x2 row major Matrix multiply adjugate A*(B#)
    static vec_t VECTORCALL Mat2MulAdj(vec_t vec1, vec_t vec2)
    {
    	return VecSub(VecMul(vec1, VecSwizzle(vec2, 3, 0, 3, 0)),
    	              VecMul(VecSwizzle(vec1, 1, 0, 3, 2), VecSwizzle(vec2, 2, 1, 2, 1)));
    }
    
    // this will not work on camera matrix this is for only transformation matricies
    static Matrix4 VECTORCALL InverseTransform(Matrix4 inM)
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
        sizeSqr = VecMul(out.r[0], out.r[0]);
        sizeSqr = VecFmadd(out.r[1], out.r[1], sizeSqr);
        sizeSqr = VecFmadd(out.r[2], out.r[2], sizeSqr);
        
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
        out.r[3] =       VecMul(out.r[0], VecSplatX(inM.r[3]));
        out.r[3] = VecFmaddLane(out.r[1], inM.r[3], out.r[3], 1);
        out.r[3] = VecFmaddLane(out.r[2], inM.r[3], out.r[3], 2);
        out.r[3] = VecSub(VecSetR(0.f, 0.f, 0.f, 1.f), out.r[3]);
        return out;
    }
    
    #if !defined(AX_ARM)
    static Matrix4 VECTORCALL Inverse(Matrix4 m)
    {
        vec_t A = VecShuffle_0101(m.r[0], m.r[1]);
        vec_t B = VecShuffle_2323(m.r[0], m.r[1]);
        vec_t C = VecShuffle_0101(m.r[2], m.r[3]);
        vec_t D = VecShuffle_2323(m.r[2], m.r[3]);
        
        vec_t detSub = VecSub(
        	VecMul(VecShuffle(m.r[0], m.r[2], 0, 2, 0, 2), VecShuffle(m.r[1], m.r[3], 1, 3, 1, 3)),
        	VecMul(VecShuffle(m.r[0], m.r[2], 1, 3, 1, 3), VecShuffle(m.r[1], m.r[3], 0, 2, 0, 2))
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
        
        detM = VecFmadd(detB, detC, detM);
        
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
    }
    #else
    static Matrix4 VECTORCALL Inverse(Matrix4 mat)
    {
        float32x4_t v0, v1, v2, v3,
                    t0, t1, t2, t3, t4, t5,
                    x0, x1, x2, x3, x4, x5, x6, x7, x8;
        float32x4x2_t a1;
        float32x2_t   lp, ko, hg, jn, im, fe, ae, bf, cg, dh;
        float32x4_t   x9 = VecSetR(-0.f,  0.f, -0.f,  0.f);
        
        x8 = vrev64q_f32(x9);
        /* l p k o, j n i m */
        a1  = vzipq_f32(mat.r[3], mat.r[2]);
        jn  = vget_high_f32(a1.val[0]);
        im  = vget_low_f32(a1.val[0]);
        lp  = vget_high_f32(a1.val[1]);
        ko  = vget_low_f32(a1.val[1]);
        hg  = vget_high_f32(mat.r[1]);
        
        x1  = vcombine_f32(vdup_lane_f32(lp, 0), lp);                   /* l p p p */
        x2  = vcombine_f32(vdup_lane_f32(ko, 0), ko);                   /* k o o o */
        x0  = vcombine_f32(vdup_lane_f32(lp, 1), vdup_lane_f32(hg, 1)); /* h h l l */
        x3  = vcombine_f32(vdup_lane_f32(ko, 1), vdup_lane_f32(hg, 0)); /* g g k k */
        
        t0 = vmlsq_f32(vmulq_f32(x3, x1), x2, x0);
        fe = vget_low_f32(mat.r[1]);
        x4 = vcombine_f32(vdup_lane_f32(jn, 0), jn);                   /* j n n n */
        x5 = vcombine_f32(vdup_lane_f32(jn, 1), vdup_lane_f32(fe, 1)); /* f f j j */
        
        t1 = vmlsq_f32(vmulq_f32(x5, x1), x4, x0);
        t2 = vmlsq_f32(vmulq_f32(x5, x2), x4, x3);
        
        x6 = vcombine_f32(vdup_lane_f32(im, 1), vdup_lane_f32(fe, 0)); /* e e i i */
        x7 = vcombine_f32(vdup_lane_f32(im, 0), im);                   /* i m m m */
        
        t3 = vmlsq_f32(vmulq_f32(x6, x1), x7, x0);
        t4 = vmlsq_f32(vmulq_f32(x6, x2), x7, x3);
        t5 = vmlsq_f32(vmulq_f32(x6, x4), x7, x5);
        
        /* h d f b, g c e a */
        a1 = vtrnq_f32(mat.r[0], mat.r[1]);
        x4 = vrev64q_f32(a1.val[0]); /* c g a e */
        x5 = vrev64q_f32(a1.val[1]); /* d h b f */
        
        ae = vget_low_f32(x4);
        cg = vget_high_f32(x4);
        bf = vget_low_f32(x5);
        dh = vget_high_f32(x5);
        
        x0 = vcombine_f32(ae, vdup_lane_f32(ae, 1)); /* a a a e */
        x1 = vcombine_f32(bf, vdup_lane_f32(bf, 1)); /* b b b f */
        x2 = vcombine_f32(cg, vdup_lane_f32(cg, 1)); /* c c c g */
        x3 = vcombine_f32(dh, vdup_lane_f32(dh, 1)); /* d d d h */
        
        v0 = VecXor(vmlaq_f32(vmlsq_f32(vmulq_f32(x1, t0), x2, t1), x3, t2), x8);
        v2 = VecXor(vmlaq_f32(vmlsq_f32(vmulq_f32(x0, t1), x1, t3), x3, t5), x8);
        v1 = VecXor(vmlaq_f32(vmlsq_f32(vmulq_f32(x0, t0), x2, t3), x3, t4), x9);
        v3 = VecXor(vmlaq_f32(vmlsq_f32(vmulq_f32(x0, t2), x1, t4), x2, t5), x9);
        /* determinant */
        x0 = vcombine_f32(vget_low_f32(vzipq_f32(v0, v1).val[0]),
                          vget_low_f32(vzipq_f32(v2, v3).val[0]));

        Matrix4 dest;
        x0 = VecDiv(VecOne(), VecHSum(vmulq_f32(x0, mat.r[0])));
        dest.r[0] = vmulq_f32(v0, x0);
        dest.r[1] = vmulq_f32(v1, x0);
        dest.r[2] = vmulq_f32(v2, x0);
        dest.r[3] = vmulq_f32(v3, x0);
        return dest;
    }
    #endif
    
    static Matrix4 VECTORCALL Multiply(Matrix4 in1, Matrix4 in2)
    {
        vec_t m0;
        m0 = VecMul(in1.r[0], VecSplatX(in2.r[0]));
        m0 = VecFmaddLane(in1.r[1], in2.r[0], m0, 1);
        m0 = VecFmaddLane(in1.r[2], in2.r[0], m0, 2); 
        m0 = VecFmaddLane(in1.r[3], in2.r[0], m0, 3); 
        in2.r[0] = m0;
        
        m0 = VecMul(in1.r[0], VecSplatX(in2.r[1]));
        m0 = VecFmaddLane(in1.r[1], in2.r[1], m0, 1); 
        m0 = VecFmaddLane(in1.r[2], in2.r[1], m0, 2); 
        m0 = VecFmaddLane(in1.r[3], in2.r[1], m0, 3); 
        in2.r[1] = m0;
        
        m0 = VecMul(in1.r[0], VecSplatX(in2.r[2]));
        m0 = VecFmaddLane(in1.r[1], in2.r[2], m0, 1); 
        m0 = VecFmaddLane(in1.r[2], in2.r[2], m0, 2); 
        m0 = VecFmaddLane(in1.r[3], in2.r[2], m0, 3); 
        in2.r[2] = m0;
        
        m0 = VecMul(in1.r[0], VecSplatX(in2.r[3]));
        m0 = VecFmaddLane(in1.r[1], in2.r[3], m0, 1); 
        m0 = VecFmaddLane(in1.r[2], in2.r[3], m0, 2); 
        m0 = VecFmaddLane(in1.r[3], in2.r[3], m0, 3); 
        in2.r[3] = m0;
        return in2;
    }
    
    static Matrix4 VECTORCALL FromQuaternion(Quaternion q)
    {
        #if defined(AX_ARM)
        Matrix4 mat = {};
        MatrixFromQuaternion<4>(mat.GetPtr(), q);
        mat.m[3][3] = 1.0f;
        return mat;
        #else
        Matrix4 M;
        const vec_t  Constant1110 = VecSetR(1.0f, 1.0f, 1.0f, 0.0f);
        
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
        return FromQuaternion(MakeQuat(quaternion[0], quaternion[1], quaternion[2], quaternion[3]));
    }
    
    static Matrix4 VECTORCALL Transpose(Matrix4 M)
    {
        Matrix4 mResult;
        #ifdef AX_ARM
        float32x4x2_t P0 = vzipq_f32(M.r[0], M.r[2]);
        float32x4x2_t P1 = vzipq_f32(M.r[1], M.r[3]);
        float32x4x2_t T0 = vzipq_f32(P0.val[0], P1.val[0]);
        float32x4x2_t T1 = vzipq_f32(P0.val[1], P1.val[1]);
        mResult.r[0] = T0.val[0];
        mResult.r[1] = T0.val[1];
        mResult.r[2] = T1.val[0];
        mResult.r[3] = T1.val[1];
        #else
        const vec_t vTemp1 = VecShuffleR(M.r[0], M.r[1], 1, 0, 1, 0);
        const vec_t vTemp3 = VecShuffleR(M.r[0], M.r[1], 3, 2, 3, 2);
        const vec_t vTemp2 = VecShuffleR(M.r[2], M.r[3], 1, 0, 1, 0);
        const vec_t vTemp4 = VecShuffleR(M.r[2], M.r[3], 3, 2, 3, 2);
        mResult.r[0] = VecShuffleR(vTemp1, vTemp2, 2, 0, 2, 0);
        mResult.r[1] = VecShuffleR(vTemp1, vTemp2, 3, 1, 3, 1);
        mResult.r[2] = VecShuffleR(vTemp3, vTemp4, 2, 0, 2, 0);
        mResult.r[3] = VecShuffleR(vTemp3, vTemp4, 3, 1, 3, 1);
        #endif
        return mResult;
    }
    
    static vec_t VECTORCALL Vector4Transform(vec_t V, const Matrix4& M)
    {
        return ::Vector4Transform(V, M.r);
    }
};
 
struct FrustumPlanes
{
    vec_t planes[6];
};

inline FrustumPlanes CreateFrustumPlanes(const Matrix4& viewProjection)
{
    FrustumPlanes result; // normalize each plane if you want to do sphere or cone intersection
    Matrix4 C = Matrix4::Transpose(viewProjection);
    result.planes[0] = VecAdd(C.r[3], C.r[0]); // m_left_plane
    result.planes[1] = VecSub(C.r[3], C.r[0]); // m_right_plane
    result.planes[2] = VecAdd(C.r[3], C.r[1]); // m_bottom_plane
    result.planes[3] = VecSub(C.r[3], C.r[1]); // m_top_plane
    result.planes[4] = C.r[2];                 // m_near_plane  
    // result.planes[5] = VecSub(C.r[3], C.r[2]); // m_far_plane
    return result;
}

__forceinline vec_t VECTORCALL MaxPointAlongNormal(vec_t min, vec_t max, vec_t n) 
{
    return VecSelect(min, max, VecCmpGe(n, VecZero()));
}

inline bool VECTORCALL CheckAABBCulled(vec_t minAABB, vec_t maxAABB, const FrustumPlanes& frustum, const Matrix4& matrix)
{
    vec_t min = Vector3Transform(minAABB, matrix.r);
    vec_t max = Vector3Transform(maxAABB, matrix.r);
    
    for (uint i = 0u; i < 5u; ++i) // make < 6 if you want far plane 
    {
        vec_t p = MaxPointAlongNormal(min, max, frustum.planes[i]);
        if (VecDotf(frustum.planes[i], p) < 0.0f)
        {
            return false;
        }
    }
    return true;
}

inline bool isPointCulled(const FrustumPlanes& frustum, const Vector3f& _point, const Matrix4& matrix)
{
    vec_t point = Vector3Transform(VecLoad(&_point.x), matrix.r);
    
    for (uint i = 0u; i < 5u; ++i)
    {
        if (VecDotf(frustum.planes[i], point) < 0.0f)
            return false;
    }
    return true;
}

AX_END_NAMESPACE 