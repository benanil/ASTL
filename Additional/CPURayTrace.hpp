#pragma once
#include "BVH.hpp"

AX_NAMESPACE

#if defined(AX_SUPPORT_SSE) && defined(AX_BVH)

struct HitRecord
{
	Vector3f normal;
	Vector2f uv;
	float distance;
	uint color;
	uint index;
};

struct Texture {
  int width, height, offset, padd;
};

struct TextureInfo {
  char* path;
  char* name;
  uint glTextureIcon;
};

struct MaterialInfo {
  char* name;
};

// this material will store default mesh/submesh
// values and we will be able to change properties for different mesh instances
#pragma pack(push)
struct Material {
  uint color;
  uint specularColor;
  ushort albedoTextureIndex;
  ushort specularTextureIndex;
  half shininess, roughness;
};
#pragma pack(pop)

struct BVHData
{
  uint* BVHIndices;
  Tri* Triangles;
  BVHNode* BVHNodes;
  Material* Materials;
  RGB8* TexturePixels;
  Texture* Textures;
  MeshInstance* MeshInstances;
  uint NumMeshInstances;
};

extern const float RayacastMissDistance = 1e30f;

void CPU_RayTraceInitialize(const BVHData& data);

HitRecord CPU_RayCast(RaySSE ray);

#endif //AX_SUPPORT_SSE

AX_END_NAMESPACE
