#pragma once
#include "Math/Matrix.hpp"

AX_NAMESPACE

AX_ALIGNED(16) struct BVHNode
{
  union { struct { float3 aabbMin; uint leftFirst; }; __m128 minv; };
  union { struct { float3 aabbMax; uint triCount; };  __m128 maxv; };
};

#pragma pack(push)
struct RGB8 {
  unsigned char r, g, b;
};
#pragma pack(pop)

struct MeshInfo {
  uint numTriangles;
  uint triangleStart;
  ushort materialStart;
  ushort numMaterials;
  const char* path;
};

#pragma pack(push)
AX_ALIGNED(16) struct Tri {
  union { struct { float3 vertex0; float centeroidx; };  __m128 v0; };
  union { struct { float3 vertex1; float centeroidy; };  __m128 v1; };
  union { struct { float3 vertex2; float centeroidz; };  __m128 v2; };

  half uv0x, uv0y;
  half uv1x, uv1y;
  half uv2x, uv2y;
  short materialIndex;
  half normal0x, normal0y, normal0z;
  half normal1x, normal1y, normal1z;
  half normal2x, normal2y, normal2z;
};
#pragma pack(pop)

// internal struct use
struct MeshInstance {
  Matrix4 inverseTransform;
  ushort meshIndex;
  ushort materialStart; // each submesh can have material
};

typedef uint MeshInstanceHandle;

AX_END_NAMESPACE