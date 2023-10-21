
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

ParsedGLTF gltf;
Shader     shader;
Mesh       mesh;
Texture    texture;
Texture androidRobotTexture;

const char* fragmentShaderSource =
"#version 150 core\n\
out vec4 color;\
in vec2 texCoord;\
uniform sampler2D tex;\
void main(){\
    color = texture(tex, texCoord);\
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
    androidRobotTexture = LoadTexture("Textures/forest.jpg");
    fullScreenShader = CreateFullScreenShader(fragmentShaderSource);
    gltf    = ParseGLTF("Meshes/Duck.gltf");
    ASSERT(gltf.error == GLTFError_NONE);
    shader  = ImportShader("Shaders/3DFirstVert.glsl", "Shaders/3DFirstFrag.glsl");
    mesh    = CreateMeshFromGLTF(&gltf.meshes[0].primitives[0]);
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