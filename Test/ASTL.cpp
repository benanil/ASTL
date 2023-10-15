
#include "../Profiler.hpp"
// #include "AdventOfCode2021.cpp"

#include "../String.hpp"
#include "../Array.hpp"
#include <stdio.h>

#include "../Renderer.hpp"
#include "../Window.hpp"
#include "../glad.hpp"

ParsedGLTF gltf;
Shader     shader;
Mesh       mesh;
Texture    texture;
Texture androidRobotTexture;

const char* fragmentShaderSource = "#version 150 core\n\
		out vec4 color;\
		in vec2 texCoord;\
		uniform sampler2D tex;\
		void main(){\
			color = texture(tex, texCoord);\
		}";

static Shader fullScreenShader{0};

#if AX_USE_WINDOW
void AXInit()
{
    SetWindowName("Test Window");
    SetWindowSize(1920, 1080);
    SetWindowPosition(0, 0);
}

int AXStart()
{
    InitRenderer();
    androidRobotTexture = LoadTexture("Textures/forest.jpg");
    fullScreenShader = CreateFullScreenShader(fragmentShaderSource);
    gltf    = ParseGLTF("Meshes/Duck.gltf");
    shader  = ImportShader("Shaders/3DFirstVert.glsl", "Shaders/3DFirstFrag.glsl");
    mesh    = CreateMeshFromGLTF(&gltf.meshes[0]);
    texture = LoadTexture(gltf.images[0].path);
    return 0;
}

// do rendering and main loop here
void AXLoop()
{
    ToggleDepthTest(false);

    // works like a skybox
    RenderFullScreen(fullScreenShader, androidRobotTexture.handle);

    ToggleDepthTest(true);

    BindShader(shader);
    SetTexture(texture, 0);
    RenderMesh(mesh);
    Render();
}

void AXExit()
{
    DeleteShader(shader);
    DeleteMesh(mesh);
    FreeGLTF(gltf);
    DeleteTexture(androidRobotTexture);
    DeleteShader(fullScreenShader);
    DestroyRenderer();
}
#else
int main()
{
  return 0;
}
#endif

// cgltf_options options;
// memset(&options, 0, sizeof(cgltf_options));
// cgltf_result result = cgltf_parse_file(&options, "Meshes/scene.gltf", &cgltfData);
// 
// if (result == cgltf_result_success)
// result = cgltf_load_buffers(&options, cgltfData, "Meshes/scene.gltf");
// 
// if (result == cgltf_result_success)
// result = cgltf_validate(cgltfData);
// 
// printf("Result: %d\n", result);
// 
// if (result != cgltf_result_success)
// {
//   ASSERT(0);
// }

// // Iterate through primitives
// cgltf_mesh cmesh = cgltfData->meshes[0];
// for (size_t j = 0; j < cmesh.primitives_count; ++j)
// {
//   cgltf_primitive* primitive = &cmesh.primitives[j];
// 
//   // Access vertex attributes (positions, normals, texcoords, etc.)
//   cgltf_accessor* positionAccessor = primitive->attributes[0].data;
//   cgltf_accessor* normalAccessor = primitive->attributes[1].data;
//   cgltf_accessor* texcoordAccessor = primitive->attributes[2].data;
// 
//   // Access vertex data (assuming they are in float format)
//   float* positions = (float*)positionAccessor->buffer_view->buffer->data + positionAccessor->offset / sizeof(float);
//   float* normals = (float*)normalAccessor->buffer_view->buffer->data + normalAccessor->offset / sizeof(float);
//   // float* texcoords = (float*)texcoordAccessor->buffer_view->buffer->data + texcoordAccessor->offset / sizeof(float);
// 
//   cgltf_accessor* indexAccessor = primitive->indices;
//   unsigned short* indices = (unsigned short*)indexAccessor->buffer_view->buffer->data + indexAccessor->offset / sizeof(unsigned short);
//   printf("x");
// }