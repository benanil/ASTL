#ifndef AX_RENDERER_H
#define AX_RENDERER_H

#include "GLTFParser.hpp"

struct Shader
{
    unsigned int handle;
};

struct Texture
{
    int width, height;
    unsigned int handle;
};

enum GraphicType : char
{
    GraphicType_Byte,
    GraphicType_UnsignedByte, 
    GraphicType_Short,
    GraphicType_UnsignedShort, 
    GraphicType_Int,
    GraphicType_UnsignedInt,
    GraphicType_Float,
};

struct Mesh
{
    int numVertex, numIndex;
    // unsigned because opengl accepts unsigned
    unsigned int vertexLayoutHandle;
    unsigned int indexHandle;
    unsigned int indexType;  // uint32, uint64.
    unsigned int vertexHandles[6]; // opengl handles for, POSITION, TexCoord...
    int attributes; // POSITION, TexCoord...
};

// you have to set vertex attributes yourself
Mesh CreateMesh(void* vertexBuffer, void* indexBuffer, int numVertex, int numIndex, int vertexSize);

Mesh CreateMeshFromGLTF(GLTFMesh* mesh);

Texture CreateTexture(int width, int height, void* data);

Texture LoadTexture(const char* path);

Shader LoadShader(const char* vertexSource, const char* fragmentSource);

Shader CreateFullScreenShader(const char* fragmentSource);

Shader ImportShader(const char* vertexSource, const char* fragmentSource);

void DeleteTexture(Texture texture);

void DeleteShader(Shader shader);

void DeleteMesh(Mesh mesh);

void HandleInput();

void RenderFullScreen(Shader fullScreenShader, unsigned int texture);

void BindShader(Shader shader);

void SetTexture(Texture texture, int index);

void RenderMesh(Mesh mesh);

void Render();

void InitRenderer();

void DestroyRenderer();

void UpdateRenderArea();

#endif //AX_RENDERER_H