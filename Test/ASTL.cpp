
// #include "../Profiler.hpp"
// #include "AdventOfCode2021.cpp"

// #include "../String.hpp"
// #include "../Array.hpp"
#include <stdio.h>
// #include "../IO.hpp"

//#define STB_IMAGE_IMPLEMENTATION
#include "../stb_image.h"

#if AX_USE_WINDOW

#include "../Renderer.hpp"
#include "../Window.hpp"
#include "../Common.hpp"
#include "../Math/Matrix.hpp"
#include <math.h>

ParsedGLTF gltf;
Shader     shader;
Mesh*      meshes;
Texture*   textures;
Texture androidRobotTexture;

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
    SetWindowName("Test Window");
    SetWindowSize(1920, 1080);
    SetWindowPosition(0, 0);
}

int AXStart()
{
    InitRenderer();
    androidRobotTexture = LoadTexture("Textures/forest.jpg", false);
    fullScreenShader = CreateFullScreenShader(fragmentShaderSource);
    gltf    = ParseGLTF("Meshes/Sponza/scene.gltf");
    ASSERT(gltf.error == GLTFError_NONE);
    shader = ImportShader("Shaders/3DFirstVert.glsl", "Shaders/3DFirstFrag.glsl");
    meshes = (Mesh*)calloc(sizeof(Mesh) * gltf.numMeshes, 1);

    for (int i = 0; i < gltf.numMeshes; ++i)
    {
        meshes[i] = CreateMeshFromGLTF(&gltf.meshes[i].primitives[0]);
    }
    
    // textures = (Texture*)calloc(sizeof(Texture) * gltf.numImages, 1);
    // 
    // for (int i = 0; i < gltf.numImages; ++i)
    // {
    //     textures[i] = LoadTexture(gltf.images[i].path, true);
    // }
    return 0;
}

extern int windowHeight_, windowWidth_;

// do rendering and main loop here
void AXLoop()
{
    ToggleDepthTest(false);
    // works like a skybox
    RenderFullScreen(fullScreenShader, androidRobotTexture.handle);
    ToggleDepthTest(true);

    Matrix4 projection, view, model, mvp;

    if (true)//(gltf.numCameras == 0)
    {
        static float f = 1.0f; f += 0.001f;
        const float distance = 15.14159265f; // this is distance from cube but I did use pi anyways  

        float verticalFOV = 65.0f, nearClip = 0.01f, farClip = 500.0f;
        Vector3f position = MakeVec3(sinf(f) * distance, 0.0f, cosf(f) * distance);
        projection = Matrix4::PerspectiveFovRH(verticalFOV * DegToRad, windowWidth_, windowHeight_, nearClip, farClip);
        view  = Matrix4::LookAtRH(position, -Vector3f::Normalize(position), Vector3f::Up());
    }
    else
    {
        GLTFCamera camera   = gltf.cameras[gltf.nodes[3].index];
        Vector3f   position = MakeVec3(gltf.nodes[3].translation);
        Quaternion rotation = MakeQuat(gltf.nodes[3].rotation);

        projection = Matrix4::PerspectiveFovRH(camera.yFov, windowWidth_, windowHeight_, camera.zNear, camera.zFar);
        view = Matrix4::LookAtRH(position, rotation.GetForward(), rotation.GetUp());
    }
    
    BindShader(shader);

    model = Matrix4::FromPosition(1.f, 0.f, 0.f) * Matrix4::CreateScale(1.f, 1.f, 1.f) * Matrix4::FromPosition(0.0f, -1.0f, 0.0f);
    mvp = model * view * projection;

    SetModelViewProjection(&mvp.m[0][0]);
    SetModelMatrix(&model.m[0][0]);

    // SetTexture(textures[0], 0); //SetTexture(textures[material.textures[0].index], 0);
    RenderMesh(meshes[0]);

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
#else

#include <string>
#include <Windows.h>

int main()
{
    char compName[512]{0};
    unsigned long compNameLen = 512;
    int ret = GetComputerNameA(compName, &compNameLen);
    int er = GetLastError();
    printf("computer name: %s", compName);
    
    std::string startText = 
    "object IntroForm: TIntroForm\n"
    "Left = 189\n"
    "Top = 225\n"
    "AlphaBlendValue = 100\n"
    "BorderStyle = bsNone\n"
    "Caption = 'IntroForm'\n"
    "ClientHeight = 230\n"
    "ClientWidth = 661\n"
    "Color = clBtnFace\n"
    "Font.Charset = DEFAULT_CHARSET\n"
    "Font.Color = clWindowText\n"
    "Font.Height = -11\n"
    "Font.Name = 'MS Sans Serif'\n"
    "Font.Style = []\n"
    "FormStyle = fsStayOnTop\n"
    "OldCreateOrder = False\n"
    "Position = poDesktopCenter\n"
    "Visible = True\n"
    "OnClose = FormClose\n"
    "OnCreate = FormCreate\n"
    "PixelsPerInch = 96\n"
    "TextHeight = 13\n"
    "object Image1: TImage\n"
    "  Left = 0\n"
    "  Top = 0\n"
    "  Width = 666\n"
    "  Height = 234\n"
    "  Picture.Data = {";
    
    int x, y, comp;
    unsigned char* img = stbi_load("banner.png", &x, &y, &comp, 4);
    std::string imgBin;
    
    const char* binToCharMap = "0123456789ABCDEF";
    for (int i = 0; i < x * y * 4;)
    {
        if ((i & 31) == 0) imgBin += "\n      ";
        imgBin.push_back(binToCharMap[img[i] & 0xF]); img[i] >>= 4;
        imgBin.push_back(binToCharMap[img[i]]); i++;

        if ((i & 31) == 0) imgBin += "\n      ";
        imgBin.push_back(binToCharMap[img[i] & 0xF]); img[i] >>= 4;
        imgBin.push_back(binToCharMap[img[i]]); i++;

        if ((i & 31) == 0) imgBin += "\n      ";
        imgBin.push_back(binToCharMap[img[i] & 0xF]); img[i] >>= 4;
        imgBin.push_back(binToCharMap[img[i]]); i++;

        if ((i & 31) == 0) imgBin += "\n      ";
        imgBin.push_back(binToCharMap[img[i] & 0xF]); img[i] >>= 4;
        imgBin.push_back(binToCharMap[img[i]]); i++;
    }

    std::string endText =
      "}\nend\n"
      "  object Timer1: TTimer\n"
      "    Interval = 2000\n"
      "    OnTimer = Timer1Timer\n"
      "    Left = 616\n"
      "    Top = 16\n"
      "  end\n"
      "end\n";

    std::string result = startText + imgBin + endText;
    // WriteAllBytes("Test.dfm", result.data(), result.size());
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