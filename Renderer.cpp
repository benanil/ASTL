
//                                                                   //
// simple graphics interface that runs in varius different platforms //
//                                                                   //

#ifdef __ANDROID__
    #include <game-activity/native_app_glue/android_native_app_glue.h>
    #include <GLES3/gl32.h>
    #define STBI_NO_STDIO
    #define STBI_NEON
#else
    #define GLAD_GL_IMPLEMENTATION
    #include "glad.hpp"

    #define AX_LOG(...) printf(__VA_ARGS__)
#endif

#include "Common.hpp"
#define STBI_ASSERT(x) ASSERT(x)
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "Renderer.hpp"
#include "IO.hpp"
#include "Math/Matrix.hpp"
#include <math.h>

static uint32  emptyVao         = 0;
static Shader  fullScreenShader = {0};
static Texture androidRobotTexture{};

bool CheckAndLogGlError() 
{
    GLenum error = glGetError();
    if (error != GL_NO_ERROR) return false;
    switch (error) {
    case GL_INVALID_ENUM:                  AX_LOG("GL Error: GL_INVALID_ENUM\n"); break;
    case GL_INVALID_VALUE:                 AX_LOG("GL Error: GL_INVALID_VALUE\n"); break;
    case GL_INVALID_OPERATION:             AX_LOG("GL Error: GL_INVALID_OPERATION\n"); break;
    case GL_INVALID_FRAMEBUFFER_OPERATION: AX_LOG("GL Error: GL_INVALID_FRAMEBUFFER_OPERATION\n"); break;
    case GL_OUT_OF_MEMORY:                 AX_LOG("GL Error: GL_OUT_OF_MEMORY\n"); break;
    default: AX_LOG("Unknown GL error: %d\n", error); 
    }
    return true;
}

// for now only .jpg files supported
Texture CreateTexture(int width, int height, void* data)
{
    Texture texture;
    glGenTextures(1, &texture.handle);
    glBindTexture(GL_TEXTURE_2D, texture.handle);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    texture.width  = width;
    texture.height = height;
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
    CheckAndLogGlError();
    return texture;
}

Texture LoadTexture(const char* path)
{
    int width, height, channels;
    unsigned char* image = nullptr;
#ifdef __ANDROID__
    AAsset* asset = AAssetManager_open(g_android_app->activity->assetManager, path, 0);
    off_t size = AAsset_getLength(asset);
    unsigned char* buffer = (unsigned char*)malloc(size);
    AAsset_read(asset, buffer, size);
    image = stbi_load_from_memory(buffer, size, &width, &height, &channels, 3);
    ASSERT(image);
    stbi_image_free(image);
    free(buffer);
    AAsset_close(asset);
    return texture;

#else
    image = stbi_load(path, &width, &height, &channels, 3);
    ASSERT(image);
    Texture texture = CreateTexture(width, height, image);
    return texture;
#endif
}

static GLenum ToGLType(GraphicType type)
{
    return GL_BYTE + (GLenum)type;
}

static int TypeToSize(GraphicType type)
{
    // BYTE, UNSIGNED_BYTE, SHORT, UNSIGNED_SHORT, INT, UNSIGNED_INT, FLOAT           
    const int sizes[8]{ 1, 1, 2, 2, 4, 4, 4 };
    return sizes[(int)type - (int)GraphicType_Byte];
}

// you have to set vertex attributes yourself
Mesh CreateMesh(void* vertexBuffer, void* indexBuffer, int numVertex, int numIndex, int vertexSize)
{
    Mesh mesh;
    glGenBuffers(1, &mesh.vertexHandles[0]);
    glGenBuffers(1, &mesh.indexHandle);

    glGenVertexArrays(1, &mesh.vertexLayoutHandle);
    glBindVertexArray(mesh.vertexLayoutHandle);

    glBindBuffer(GL_ARRAY_BUFFER, mesh.vertexHandles[0]);
    glBufferData(GL_ARRAY_BUFFER, vertexSize * numVertex, vertexBuffer, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh.indexHandle);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, numIndex * sizeof(uint32), indexBuffer, GL_STATIC_DRAW);
    mesh.indexType = GL_UNSIGNED_INT;
    mesh.numIndex = numIndex; mesh.numVertex = numVertex; mesh.indexType = GL_UNSIGNED_INT;
    return mesh;
}

Mesh CreateMeshFromGLTF(GLTFMesh* gltf)
{
    Mesh mesh;
    mesh.indexType  = gltf->indexType;
    mesh.numIndex   = gltf->numIndices;
    mesh.numVertex  = gltf->numVertices;
    mesh.attributes = gltf->attributes;
    
    glGenBuffers(1, &mesh.indexHandle);
    glGenBuffers(gltf->numAttributes, mesh.vertexHandles);

    glGenVertexArrays(1, &mesh.vertexLayoutHandle);
    glBindVertexArray(mesh.vertexLayoutHandle);
    
    int attributes = mesh.attributes;
    int i = 0, v = 0;
    
    while (attributes)
    {
        glBindBuffer(GL_ARRAY_BUFFER, mesh.vertexHandles[v]);
        unsigned int type = (int)(gltf->attribTypes[i]);
        int size = TypeToSize((GraphicType)type) * gltf->attribNumComps[i];
        glBufferData(GL_ARRAY_BUFFER, size * mesh.numVertex, gltf->vertexAttribs[i], GL_STATIC_DRAW);
        glVertexAttribPointer(v, gltf->attribNumComps[i], GL_BYTE + gltf->attribTypes[i], GL_FALSE, 0, nullptr);
        glEnableVertexAttribArray(v);
        attributes &= ~1;
        int tz = TrailingZeroCount(attributes);
        attributes >>= tz;
        i += tz;
        v++;
    }

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh.indexHandle);
    int indexSize = gltf->numIndices * TypeToSize((GraphicType)(gltf->indexType - GL_BYTE));
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indexSize, gltf->indices, GL_STATIC_DRAW);

    return mesh;
}

void CheckShaderError(uint shader)
{
    GLint isCompiled = 0;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &isCompiled);
    if (isCompiled == GL_FALSE)
    {
        GLint maxLength = 0;
        glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &maxLength);
        char infoLog[1024]{};
        glGetShaderInfoLog(shader, maxLength, &maxLength, infoLog);
        AX_LOG("No compile fs %s", infoLog) ;
        glDeleteShader(shader);
        DestroyRenderer();
        ASSERT(0);
    }
}

Shader LoadShader(const char* vertexSource, const char* fragmentSource)
{
    // Vertex shader
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexSource, NULL);
    glCompileShader(vertexShader);
    CheckShaderError(vertexShader);
    // Check for compile time errors
    GLint success;
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success); ASSERT(success);
    // Fragment shader
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentSource, NULL);
    glCompileShader(fragmentShader);
    CheckShaderError(fragmentShader);
    // Check for compile time errors
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success); ASSERT(success);
    // Link shaders
    uint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    // Check for linking errors
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success); ASSERT(success);

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    glUseProgram(shaderProgram);
    return {shaderProgram};
}

Shader CreateFullScreenShader(const char* fragmentSource)
{
    const GLchar* vertexShaderSource = "#version 150 core\n\
    out vec2 texCoord;\
    void main(){\
    	float x = -1.0 + float((gl_VertexID & 1) << 2);\
    	float y = -1.0 + float((gl_VertexID & 2) << 1);\
    	texCoord.x = (x + 1.0) * 0.5;\
    	texCoord.y = (y + 1.0) * 0.5;\
    	texCoord.y = 1.0 - texCoord.y;\
        gl_Position = vec4(x, y, 0, 1);\
    }";
    return LoadShader(vertexShaderSource, fragmentSource);
}

Shader ImportShader(const char* vertexPath, const char* fragmentPath)
{
    char* vertexText = ReadAllFile(vertexPath);
    char* fragmentText = ReadAllFile(fragmentPath);
    Shader shader = LoadShader(vertexText, fragmentText);
    free(vertexText);
    free(fragmentText);
    return shader;
}

void DeleteTexture(Texture texture) { glDeleteTextures(1, &texture.handle); }
void DeleteShader(Shader shader)    { glDeleteProgram(shader.handle);       }

void DeleteMesh(Mesh mesh)
{
    glDeleteVertexArrays(1, &mesh.vertexLayoutHandle);
    glDeleteBuffers(PopCount(mesh.attributes), mesh.vertexHandles);
    glDeleteBuffers(1, &mesh.indexHandle);
}

void GLDebugMessageCallback(GLenum source, GLenum type, GLuint id, GLenum severity,
                            GLsizei length, const GLchar* msg, const void* data) {
    AX_LOG("OpenGL error: %s", msg);
    AX_LOG("\n");
}

void InitRenderer()
{
#ifdef __ANDROID__
    InitWindow();
#endif
    androidRobotTexture = LoadTexture("Textures/forest.jpg");

	const GLchar* fragmentShaderSource = "#version 150 core\n\
		out vec4 color;\
		in vec2 texCoord;\
		uniform sampler2D tex;\
		void main(){\
			color = texture(tex, texCoord);\
		}";
    
    fullScreenShader = CreateFullScreenShader(fragmentShaderSource);
    // create empty vao unfortunately this step is necessary for ogl 3.2
    glGenVertexArrays(1, &emptyVao);
    // setup any other gl related global states
    glClearColor(0.2f, 0.8f, 0.25f, 1.0f);
}

void DestroyRenderer()
{
    DeleteTexture(androidRobotTexture);
    DeleteShader(fullScreenShader);
}

static unsigned int currentShader = 0;

void RenderFullScreen(Shader fullScreenShader, unsigned int texture)
{
    glUseProgram(fullScreenShader.handle);
    glBindVertexArray(emptyVao);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, texture);
    glUniform1i(glGetUniformLocation(currentShader, "tex"), 0);
    glDrawArrays(GL_TRIANGLES, 0, 3);
}

void BindShader(Shader shader)
{
    glUseProgram(shader.handle);
    currentShader = shader.handle;
}

void SetTexture(Texture texture, int index)
{
    glActiveTexture(GL_TEXTURE0 + index);
    glBindTexture(GL_TEXTURE_2D, texture.handle);
}

extern int windowHeight_, windowWidth_;

void RenderMesh(Mesh mesh)
{
    glBindVertexArray(mesh.vertexLayoutHandle);

    static float f = 1.0f; f += 0.01f;
    const float distance = 3.14159265f; // this is distance from cube but I did use pi anyways  
    Vector3f position    = MakeVec3(sinf(f) * distance, 0.0f, cosf(f) * distance);
    float verticalFOV    = 65.0f, nearClip = 0.01f, farClip = 500.0f;

    Matrix4 projection = Matrix4::PerspectiveFovRH(verticalFOV * DegToRad, windowWidth_, windowHeight_, nearClip, farClip);
    Matrix4 model      = Matrix4::FromPosition(0.0f, 0.0f, 0.0f);
    Matrix4 view       = Matrix4::LookAtRH(position + Vector3f::Up(), -Vector3f::Normalize(position), Vector3f::Up());
    Matrix4 mvp        = model * view * projection;

    glEnable(GL_CULL_FACE);
    glFrontFace(GL_CCW);

    GLint mvpLoc   = glGetUniformLocation(currentShader, "mvp");
    GLint modelLoc = glGetUniformLocation(currentShader, "model");

    glUniformMatrix4fv(mvpLoc  , 1, false, &mvp.m[0][0]);
    glUniformMatrix4fv(modelLoc, 1, false, &model.m[0][0]);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh.indexHandle);
    glDrawElements(GL_TRIANGLES, mesh.numIndex, mesh.indexType, nullptr);
}

void Render()
{
    UpdateRenderArea();

    // RenderFullScreen(fullScreenShader, androidRobotTexture.handle);
    // Present the rendered image. This is an implicit glFlush.
#ifdef __ANDROID__
    EGLBoolean swapResult = eglSwapBuffers(display_, surface_);
    ASSERT(swapResult);
    glClear(GL_COLOR_BUFFER_BIT);
#endif
}
