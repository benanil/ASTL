
// #include "../Profiler.hpp"
// #include "AdventOfCode2021.cpp"

// #include "../String.hpp"
// #include "../Array.hpp"
#include <stdio.h>
// #include "../IO.hpp"

//#define STB_IMAGE_IMPLEMENTATION
#include "../stb_image.h"

#include "../Math/Camera.hpp"
#include "../Renderer.hpp"
#include "../Window.hpp"
#include "../Common.hpp"
#include "../Math/Matrix.hpp"
#include <math.h>

const Vector2i windowStartSize{1920, 1080};

ParsedGLTF gltf;
Shader     shader;
Mesh*      meshes;
Texture*   textures;
Texture androidRobotTexture;
Camera camera(MakeVec2(1920, 1080));

const char* fragmentShaderSource =
"#version 150 core\n\
out vec4 color;\n\
in vec2 texCoord;\n\
uniform sampler2D tex;\n\
void main(){\n\
    vec4 texColor = texture(tex, texCoord);\n\
    // Apply gamma correction to the color channels\n\
    color = vec4(pow(texColor.rgb, vec3(1.0 / 2.18)), texColor.a);\n\
}";
static Shader fullScreenShader{0};

void AXInit()
{
    SetWindowName("Duck Window");
    SetWindowSize(windowStartSize.x, windowStartSize.y);
    SetWindowPosition(0, 0);
}

void WindowResizeCallback(int width, int height)
{
    camera.RecalculateProjection(width, height);
}

int AXStart()
{
    InitRenderer();
    androidRobotTexture = LoadTexture("Textures/forest.jpg", false);
    fullScreenShader = CreateFullScreenShader(fragmentShaderSource);
    printf("parsing");
    gltf    = ParseGLTF("Meshes/Duck.gltf");
    if (gltf.error != GLTFError_NONE)
    {
        printf("gltf mesh parse failed");
        return 0;
    }
    ASSERT(gltf.error == GLTFError_NONE);
    shader = ImportShader("Shaders/3DFirstVert.glsl", "Shaders/3DFirstFrag.glsl");
    meshes = (Mesh*)calloc(sizeof(Mesh) * gltf.numMeshes, 1);

    for (int i = 0; i < gltf.numMeshes; ++i)
    {
        meshes[i] = CreateMeshFromGLTF(&gltf.meshes[i].primitives[0]);
    }
    
    textures = (Texture*)calloc(sizeof(Texture) * gltf.numImages, 1);
    
    for (int i = 0; i < gltf.numImages; ++i)
    {
        textures[i] = LoadTexture(gltf.images[i].path, true);
    }
    
    const float distance = 15.14159265f; // this is distance from cube but I did use pi anyways  
    camera.position = MakeVec3(Sin(1.0f) * distance, 0.0f, Cos(1.0f) * distance);
    camera.RecalculateProjection(windowStartSize.x, windowStartSize.y);
    camera.RecalculateView();
    return 0;
}

extern int windowHeight_, windowWidth_;
static Vector3f meshPosition{0.0f, -1.0f, 0.0f};

// do rendering and main loop here
void AXLoop()
{
    ToggleDepthTest(false);
    // works like a skybox
    RenderFullScreen(fullScreenShader, androidRobotTexture.handle);
    ToggleDepthTest(true);

    BindShader(shader);
    
    if (GetMousePressed(MouseButton_Left)) meshPosition.x += 2.1f;
    if (GetMousePressed(MouseButton_Right)) meshPosition.x -= 2.1f;
    camera.Update();

    Matrix4 model = Matrix4::CreateScale(0.05f, 0.05f, 0.05f) * Matrix4::FromPosition(meshPosition);
    Matrix4 mvp = model * camera.view * camera.projection;

    SetModelViewProjection(&mvp.m[0][0]);
    SetModelMatrix(&model.m[0][0]);

    SetTexture(textures[0], 0); //SetTexture(textures[material.textures[0].index], 0);
    RenderMesh(meshes[0]);

    Render();
}

void AXExit()
{
    DeleteShader(shader);

    free(meshes);
    free(textures);

    FreeGLTF(gltf);
    DeleteTexture(androidRobotTexture);
    DeleteShader(fullScreenShader);
    DestroyRenderer();
}

#define AX_USE_WINDOW

#ifndef AX_USE_WINDOW

int main()
{
    AXStart();
    return 0;
}
#endif

// in update before RenderFunction
// for (int i = 0; i < gltf.numNodes; i++) 
// {
//     GLTFNode node = gltf.nodes[i];
//     // in this scene first mesh is shitty so pass that
//     if (node.type != 0) continue;
// 
//     model = Matrix4::CreateScale(1.0f, 1.0f, 1.0f) * Matrix4::FromQuaternion(node.rotation) * Matrix4::FromPosition(node.translation);
//     mvp = model * view * projection;
//     
//     SetModelViewProjection(&mvp.m[0][0]);
//     SetModelMatrix(&model.m[0][0]);
// 
//     GLTFMesh mesh = gltf.meshes[node.index];
//     for (int j = 0; j < mesh.numPrimitives; ++j)
//     {
//         GLTFMaterial material = gltf.materials[mesh.primitives[j].material];
//         SetTexture(textures[0], 0); //SetTexture(textures[material.textures[0].index], 0);
//         RenderMesh(meshes[node.index]);
//     }
// }

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