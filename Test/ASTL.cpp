
#include "../Profiler.hpp"
// #include "AdventOfCode2021.cpp"

#include "../String.hpp"
#include "../Array.hpp"
#include <stdio.h>

#include "../Renderer.hpp"
#include "../Window.hpp"

#define CGLTF_IMPLEMENTATION
#include "../cgltf.h"

ParsedGLTF gltf;
Shader     shader;
Mesh       mesh;
Texture    texture;

cgltf_data* cgltfData = NULL;

#if AX_USE_WINDOW

void AXInit()
{
    SetWindowName("Anilcan Test");
    SetWindowSize(1920, 1080);
    SetWindowPosition(0, 0);
    
    cgltf_options options;
    memset(&options, 0, sizeof(cgltf_options));
    cgltf_result result = cgltf_parse_file(&options, "Meshes/Cube.gltf", &cgltfData);

    if (result == cgltf_result_success)
        result = cgltf_load_buffers(&options, cgltfData, "Meshes/Cube.gltf");
    
    if (result == cgltf_result_success)
        result = cgltf_validate(cgltfData);

    printf("Result: %d\n", result);

    if (result != cgltf_result_success)
    {
      ASSERT(0);
    }
}

int AXStart()
{
    InitRenderer();
    gltf   = ParseGLTF("Meshes/Cube.gltf");

    // Iterate through primitives
    cgltf_mesh cmesh = cgltfData->meshes[0];
    for (size_t j = 0; j < cmesh.primitives_count; ++j)
    {
      cgltf_primitive* primitive = &cmesh.primitives[j];

      // Access vertex attributes (positions, normals, texcoords, etc.)
      cgltf_accessor* positionAccessor = primitive->attributes[0].data;
      cgltf_accessor* normalAccessor   = primitive->attributes[1].data;
      cgltf_accessor* texcoordAccessor = primitive->attributes[2].data;

      // Access vertex data (assuming they are in float format)
      float* positions = (float*)positionAccessor->buffer_view->buffer->data + positionAccessor->offset / sizeof(float);
      float* normals   = (float*)normalAccessor->buffer_view->buffer->data + normalAccessor->offset / sizeof(float);
      // float* texcoords = (float*)texcoordAccessor->buffer_view->buffer->data + texcoordAccessor->offset / sizeof(float);
      
      cgltf_accessor* indexAccessor = primitive->indices;
      unsigned short* indices       = (unsigned short*)indexAccessor->buffer_view->buffer->data + indexAccessor->offset / sizeof(unsigned short);
      printf("x");
    }

    shader  = ImportShader("Shaders/3DFirstVert.glsl", "Shaders/3DFirstFrag.glsl");
    mesh    = CreateMeshFromGLTF(&gltf.meshes[0]);
    texture = LoadTexture(gltf.images[0].path);

    BindShader(shader);
    return 0;
}

// do rendering and main loop here
void AXLoop()
{
    SetTexture(texture, 0);
    RenderMesh(mesh);
    Render();
}

void AXExit()
{
    DeleteShader(shader);
    DeleteMesh(mesh);
    FreeGLTF(gltf);
    DestroyRenderer();
}
#else
int main()
{
  return 0;
}
#endif

