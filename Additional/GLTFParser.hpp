
#pragma once

enum AAttribType_
{
    AAttribType_POSITION   = 1 << 0,
    AAttribType_TEXCOORD_0 = 1 << 1,
    AAttribType_NORMAL     = 1 << 2,
    AAttribType_TANGENT    = 1 << 3,
    AAttribType_TEXCOORD_1 = 1 << 4
};
typedef int AAttribType;

enum AErrorType_
{
    AError_NONE,
    AError_UNKNOWN,
    AError_UNKNOWN_ATTRIB,
    AError_UNKNOWN_MATERIAL_VAR,
    AError_UNKNOWN_PBR_VAR,
    AError_UNKNOWN_NODE_VAR,
    AError_UNKNOWN_TEXTURE_VAR,
    AError_UNKNOWN_ACCESSOR_VAR,
    AError_UNKNOWN_BUFFER_VIEW_VAR,
    AError_UNKNOWN_MESH_VAR,
    AError_UNKNOWN_CAMERA_VAR,
    AError_UNKNOWN_MESH_PRIMITIVE_VAR,
    AError_BUFFER_PARSE_FAIL,
    AError_BIN_NOT_EXIST,
    AError_FILE_NOT_FOUND,
    AError_UNKNOWN_DESCRIPTOR,
    AError_HASH_COLISSION,
    AError_NON_UTF8,
    AError_EXT_NOT_SUPPORTED, // scenes other than GLTF, OBJ or Fbx
    AError_MAX
};
typedef int AErrorType;

typedef struct AMaterial_
{
    // original value multiplied by 400 to get real value: float scale = ((float)scale) / 400.0f; 
    typedef short float16;
    
    struct Texture
    {
        float16 scale;    
        float16 strength; 
        short index; // index to texture path
        short texCoord;
    } textures[3]; // normalTexture, occlusionTexture, emissiveTexture

    // pbrMetallicRoughness in gltf. but baseColorFactor is below
    Texture baseColorTexture;
    Texture specularTexture;
    Texture metallicRoughnessTexture;
    float16 metallicFactor;
    float16 roughnessFactor;

    char* name;
    float16 emissiveFactor[3]; 
    float16 specularFactor;
    unsigned diffuseColor, specularColor, baseColorFactor;
    float alphaCutoff;
    bool doubleSided;
    
#ifdef __cplusplus
    const Texture& GetNormalTexture()    const { return textures[0]; }
    const Texture& GetOcclusionTexture() const { return textures[1]; }
    const Texture& GetEmissiveTexture()  const { return textures[2]; }
    Texture& GetNormalTexture()    { return textures[0]; }
    Texture& GetOcclusionTexture() { return textures[1]; }
    Texture& GetEmissiveTexture()  { return textures[2]; }
#endif

} AMaterial;

typedef struct AImage_
{
    char* path;
} AImage;

typedef struct ANode_
{
    // Warning! order is important, we copy memory in assetmanager.cpp from fbx transform
    float translation[3];
    float rotation[4];
    float scale[3];

    int   type;  // 0 mesh or 1 camera
    int   index; // index of mesh or camera
    char* name;
    int   numChildren;
    int*  children;
} ANode;

typedef struct APrimitive_
{
    // pointers to binary file to lookup position, texture, normal..
    void* indices; 
    void* vertices;
    
    int   attributes; // AAttribType Position, Normal, TexCoord, Tangent, masks
    int   indexType; // GL_UNSIGNED_INT, GL_UNSIGNED_BYTE.. 
    int   numIndices;
    int   numVertices;

    // internal use only. after parsing this is useless
    int vertexAttribs[6]; //<-when we are parsing we use this as an indicator to accessor.
    short indiceIndex; // indice index to accessor
    char  material;    // material index
    char  mode;        // 4 is triangle
} APrimitive;

typedef struct AMesh_
{
    char* name;  
    APrimitive* primitives;
    int numPrimitives;
} AMesh;

typedef struct ATexture_
{
    int   sampler;
    int   source;
    char* name;
} ATexture;

typedef struct ACamera_
{
    union {
        struct { float aspectRatio, yFov; };
        struct { float xmag, ymag; };
    };
    float zFar, zNear;
    int type; // 0 orthographic, 1 perspective
    char* name;
} ACamera;

typedef struct ASampler_
{
    char magFilter; // 0 = GL_NEAREST = 0x2600 = 9728, or 1 = 0x2601 = GL_LINEAR = 9729
    char minFilter; // 0 or 1 like above
    char wrapS; // 0 GL_REPEAT, 1 GL_CLAMP_TO_EDGE, 2 GL_CLAMP_TO_BORDER, 3 GL_MIRRORED_REPEAT
    char wrapT; //       10497               33071                 33069                 33648
} ASampler;

static_assert(sizeof(ASampler) == sizeof(int), "size must be 4");

typedef struct AScene_
{
    char* name;
    int   numNodes;   
    int*  nodes;
} AScene;

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
    
    AErrorType error;

    void* stringAllocator;
    void* intAllocator;
    void* allVertices;
    void* allIndices;

    AMesh*     meshes;
    ANode*     nodes;
    AMaterial* materials;
    ATexture*  textures;
    AImage*    images;
    ASampler*  samplers;
    ACamera*   cameras;
    AScene*    scenes;
} ParsedGLTF;

typedef struct ParsedObj_
{
    short numMeshes;
    short numMaterials;
    short numImages;
    
    AErrorType error;

    void* allVertices;
    void* allIndices;

    char* materialText; // we hold material names in this.

    AMesh*     meshes;
    AMaterial* materials;
    AImage*    images;
} ParsedObj;

// if there is an error error will be minus GLTFErrorType
// out scene should not be null
int ParseGLTF(const char* path, ParsedGLTF* scene, float scale);

int ParseObj(const char* path, ParsedObj* scene);


// Free

void FreeParsedGLTF(ParsedGLTF* gltf);

const char* ParsedSceneGetError(AErrorType error);
