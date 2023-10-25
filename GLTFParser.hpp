
#pragma once

enum GLTFAttribType
{
    GLTFAttribType_POSITION   = 1 << 0,
    GLTFAttribType_TEXCOORD_0 = 1 << 1,
    GLTFAttribType_NORMAL     = 1 << 2,
    GLTFAttribType_TANGENT    = 1 << 3,
    GLTFAttribType_TEXCOORD_1 = 1 << 4
};

enum GLTFErrorType
{
    GLTFError_NONE,
    GLTFError_UNKNOWN,
    GLTFError_UNKNOWN_ATTRIB,
    GLTFError_UNKNOWN_MATERIAL_VAR,
    GLTFError_UNKNOWN_PBR_VAR,
    GLTFError_UNKNOWN_NODE_VAR,
    GLTFError_UNKNOWN_TEXTURE_VAR,
    GLTFError_UNKNOWN_ACCESSOR_VAR,
    GLTFError_UNKNOWN_BUFFER_VIEW_VAR,
    GLTFError_UNKNOWN_MESH_VAR,
    GLTFError_UNKNOWN_CAMERA_VAR,
    GLTFError_UNKNOWN_MESH_PRIMITIVE_VAR,
    GLTFError_BUFFER_PARSE_FAIL,
    GLTFError_BIN_NOT_EXIST,
    GLTFError_FILE_NOT_FOUND,
    GLTFError_UNKNOWN_DESCRIPTOR,
    GLTFError_MAX
};

typedef struct GLTFMaterial_
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
} GLTFMaterial;

typedef struct GLTFImage_
{
    char* path;
} GLTFImage;

typedef struct GLTFNode_
{
    int   type;  // 0 mesh or 1 camera
    int   index; // index of mesh or camera
    char* name;
    float translation[3];
    float rotation[4];
    float scale[3];
    int   numChildren;
    int*  children;
} GLTFNode;

typedef struct GLTFPrimitive_
{
    // pointers to binary file to lookup position, texture, normal..
    void* vertexAttribs[6]; //<-when we are parsing we use this as an indicator to accessor.
    void* indices;
    int   attributes; // GLTFAttribType Position, Normal, TexCoord, Tangent, masks
    int   indexType; // GL_UNSIGNED_INT, GL_UNSIGNED_BYTE.. 
    int   numIndices;
    int   numVertices;

    // internal use only. after parsing this is useless
    short indiceIndex; // indice index to accessor
    char  material;    // material index
    char  mode;        // 4 is triangle
} GLTFPrimitive;

typedef struct GLTFMesh_
{
    char* name;  
    GLTFPrimitive* primitives;
    unsigned long long numPrimitives;
} GLTFMesh;

typedef struct GLTFTexture_
{
    int   sampler;
    int   source;
    char* name;
} GLTFTexture;

typedef struct GLTFCamera_
{
    union {
        float aspectRatio, yFov;
        float xmag, ymag;
    };
    float zFar, zNear;
    int type; // 0 orthographic, 1 perspective
    char* name;
} GLTFCamera;

typedef struct GLTFSampler_
{
    char magFilter; // 0 = GL_NEAREST = 0x2600 = 9728, or 1 = 0x2601 = GL_LINEAR = 9729
    char minFilter; // 0 or 1 like above
    char wrapS; // 0 GL_REPEAT, 1 GL_CLAMP_TO_EDGE, 2 GL_CLAMP_TO_BORDER, 3 GL_MIRRORED_REPEAT
    char wrapT; //       10497               33071                 33069                 33648
} GLTFSampler;

typedef struct GLTFScene_
{
    char* name;
    int   numNodes;   
    int*  nodes;
} GLTFScene;

typedef struct ParsedGLTF_
{
    short numMeshes;
    short numNodes;
    short numMaterials;
    short numTextures;
    short numImages;
    short numSamplers;
    short numCameras;
    short numScenes;

    short defaultSceneIndex;
    GLTFErrorType error;

    void* stringAllocator;
    void* intAllocator;

    GLTFMesh*     meshes;
    GLTFNode*     nodes;
    GLTFMaterial* materials;
    GLTFTexture*  textures;
    GLTFImage*    images;
    GLTFSampler*  samplers;
    GLTFCamera*   cameras;
    GLTFScene*    scenes;
} ParsedGLTF;

// if there is an error error will be minus GLTFErrorType
ParsedGLTF ParseGLTF(const char* path);

void FreeGLTF(ParsedGLTF& gltf);
