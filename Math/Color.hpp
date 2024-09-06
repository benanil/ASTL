#pragma once

#include "Vector.hpp"
#include "SIMDVectorMath.hpp"


pureconst uint PackColorToUint(uint8 r, uint8 g, uint8 b, uint8 a) {
    return r | (uint(g) << 8) | (uint(b) << 16) | (uint(a) << 24);
}

pureconst uint MakeRGBAGrayScale(uint8 gray) {
    return uint(gray) * 0x01010101u;
}

pureconst uint PackColorToUint(float r, float g, float b) {
    return (uint)(r * 255.0f) | ((uint)(g * 255.0f) << 8) | ((uint)(b * 255.0f) << 16);
}

pureconst uint PackColor3PtrToUint(float* c) {
    return (uint)(*c * 255.0f) | ((uint)(c[1] * 255.0f) << 8) | ((uint)(c[2] * 255.0f) << 16);
}

pureconst uint PackColor4PtrToUint(float* c) {
    return (uint)(*c * 255.0f) | ((uint)(c[1] * 255.0f) << 8) | ((uint)(c[2] * 255.0f) << 16) | ((uint)(c[3] * 255.0f) << 24);
}

forceinline void UnpackColor3Uint(unsigned color, float* colorf) {
    const float toFloat = 1.0f / 255.0f;
    colorf[0] = float(color >> 0  & 0xFF) * toFloat;
    colorf[1] = float(color >> 8  & 0xFF) * toFloat;
    colorf[2] = float(color >> 16 & 0xFF) * toFloat;
}

forceinline void UnpackColor4Uint(unsigned color, float* colorf) {
    const float toFloat = 1.0f / 255.0f;
    colorf[0] = float(color >> 0  & 0xFF) * toFloat;
    colorf[1] = float(color >> 8  & 0xFF) * toFloat;
    colorf[2] = float(color >> 16 & 0xFF) * toFloat;
    colorf[3] = float(color >> 24) * toFloat;
}

purefn uint MultiplyU32Colors(uint a, uint b)
{
    uint result = 0u;
    result |= ((a & 0xffu) * (b & 0xffu)) >> 8u;
    result |= ((((a >> 8u) & 0xffu) * ((b >> 8u) & 0xffu)) >> 8u) << 8u;
    result |= ((((a >> 16u) & 0xffu) * ((b >> 16u) & 0xffu)) >> 8u) << 16u;
    return result;
}

pureconst Vector3f HUEToRGB(float h) {
    float r = Clamp01(Abs(h * 6.0f - 3.0f) - 1.0f);
    float g = Clamp01(2.0f - Abs(h * 6.0f - 2.0f));
    float b = Clamp01(2.0f - Abs(h * 6.0f - 4.0f));
    return { r, g, b };
}

// converts hue to rgb color
pureconst uint32 HUEToRGBU32(float h) {
    Vector3f v3 = HUEToRGB(h);
    uint32 res = PackColor3PtrToUint(v3.arr);
    return res | 0xFF000000u; // make the alpha 255
}

pureconst Vector3f RGBToHSV(Vector3f rgb)
{
    float r = rgb.x, g = rgb.y, b = rgb.z;
    float K = 0.0f;
    if (g < b) {
        Swap(g, b);
        K = -1.0f;
    }

    if (r < g) {
        Swap(r, g);
        K = -2.0f / 6.0f - K;
    }
    const float chroma = r - (g < b ? g : b);
    return {
        Abs(K + (g - b) / (6.0f * chroma + 1e-20f)),
        chroma / (r + 1e-20f),
        r
    };
}

forceinline void HSVToRGB(Vector3f hsv, float* dst)
{
    const Vector4x32f K = VecSetR(1.0f, 2.0f / 3.0f, 1.0f / 3.0f, 3.0f);
    Vector4x32f p = VecFabs(VecSub(VecMul(VecFract(VecAdd(VecSet1(hsv.x), K)), VecSet1(6.0f)), VecSet1(3.0f)));
    Vector4x32f kx = VecSplatX(K);
    Vector4x32f rv = VecMul(VecLerp(kx, VecClamp01(VecSub(p, kx)), hsv.y), VecSet1(hsv.z));
    Vec3Store(dst, rv);
}
