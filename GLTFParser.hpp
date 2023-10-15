
#pragma once

enum GLTFAttribType
{
    GLTFAttribType_POSITION   = 1 << 0,
    GLTFAttribType_NORMAL     = 1 << 1,
    GLTFAttribType_TEXCOORD_0 = 1 << 2,
    GLTFAttribType_TANGENT    = 1 << 3,
    GLTFAttribType_TEXCOORD_1 = 1 << 4
};

enum GLTFErrorType
{
    GLTFError_NONE,
    GLTFError_UNKNOWN,
    GLTFError_UNKNOWN_ATTRIB,
    GLTFError_UNKNOWN_DESCRIPTOR,
    GLTFError_MAX
};

struct GLTFMaterial
{
    // original value multiplied by 100 to get real value: float scale = ((float)scale) / 100.0f; 
    typedef short float16;
    char* name;
    float16 emissiveFactor[3]; 
    
    struct Texture
    {
        float16 scale;    
        float16 strength; 
        char index;
        char texCoord;
    } textures[3]; // normalTexture, occlusionTexture, emissiveTexture

    struct MetallicRoughness
    {
        Texture baseColorTexture;
        Texture metallicRoughnessTexture;
        float16 baseColorFactor[4];
        float16 metallicFactor;
        float16 roughnessFactor;
    } metallicRoughness;

    bool doubleSided;
};

struct GLTFImage
{
    char* path;
};

struct GLTFNode
{
    int   type;  // 0 mesh or 1 camera
    int   index; // index of mesh or camera
    char* name;
    float translation[3];
    float rotation[4];
    float scale[3];
    int   numChildren;
    int*  children;
};

struct GLTFMesh
{
    char  attribIndices[8]; // internal use only, +2 for padding
    void* indices;
    void* vertexAttribs[6];
    int   attributes; // GLTFAttribType Position, Normal, TexCoord, Tangent, masks
    char  attribTypes[8]; // + 2 for padding
    char  attribNumComps[8]; // + 2 for padding
    char* name;

    int   indexType;
    int   numVertices;
    int   numIndices;

    char  numAttributes;
    // internal use only. after parsing this is useless
    char  indiceIndex; // indice index 
    char  material; // material index
    char  mode; // 4 is triangle
};

struct GLTFTexture
{
    int   sampler;
    int   source;
    char* name;
};

struct GLTFCamera
{
    float aspectRatio, yFov, zFar, zNear;
    int type; // 0 orthographic, 1 perspective
};

struct ParsedGLTF
{
    int numMeshes;
    int numNodes;
    int numMaterials;
    int numTextures;
    int numImages;

    void* stringAllocator;
    void* intAllocator;

    GLTFMesh*     meshes;
    GLTFNode*     nodes;
    GLTFMaterial* materials;
    GLTFTexture*  textures;
    GLTFImage*    images;

    GLTFErrorType error;
};

// if there is an error error will be minus GLTFErrorType
ParsedGLTF ParseGLTF(const char* path);

void FreeGLTF(ParsedGLTF& gltf);
