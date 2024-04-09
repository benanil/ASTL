
#ifndef ASTL_GLTF_PARSER
#define ASTL_GLTF_PARSER

// Warning! order is important, dependent in functions: ParseMeshes, ParseAttributes, ParseGLTF... 
enum AAttribType_
{
    AAttribType_POSITION   = 1 << 0,
    AAttribType_TEXCOORD_0 = 1 << 1,
    AAttribType_NORMAL     = 1 << 2,
    AAttribType_TANGENT    = 1 << 3,
    AAttribType_TEXCOORD_1 = 1 << 4,
    AAttribType_JOINTS     = 1 << 5,
    AAttribType_WEIGHTS    = 1 << 6,
    AAttribType_Count      = 7 ,
    AAttribType_MAKE32BIT  = 1 << 31
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

enum AMaterialAlphaMode_
{
    AMaterialAlphaMode_Opaque, AMaterialAlphaMode_Blend, AMaterialAlphaMode_Mask
};
typedef int AMaterialAlphaMode;

typedef struct AMaterial_
{
    // original value multiplied by 400 to get real value: float scale = ((float)scale) / 400.0f; 
    typedef short float16;
    
    struct Texture
    {
        float16 scale;    
        float16 strength; 
        short index; // index to texture path. -1 if is not exist
        short texCoord;
    } textures[3]; // normalTexture, occlusionTexture, emissiveTexture

    // pbrMetallicRoughness in gltf. but baseColorFactor is below
    Texture baseColorTexture;
    Texture specularTexture;
    Texture metallicRoughnessTexture;
    unsigned short metallicFactor;  // (float)(metallicFactor) / UINT16_MAX to get value
    unsigned short roughnessFactor; // (float)(roughnessFactor) / UINT16_MAX to get value

    char* name;
    float16 emissiveFactor[3];
    float16 specularFactor;
    unsigned diffuseColor, specularColor, baseColorFactor;
    float alphaCutoff;
    bool doubleSided;
    AMaterialAlphaMode alphaMode;

#ifdef __cplusplus
    const Texture& GetNormalTexture()    const { return textures[0]; }
    const Texture& GetOcclusionTexture() const { return textures[1]; }
    const Texture& GetEmissiveTexture()  const { return textures[2]; }
    Texture& GetNormalTexture()    { return textures[0]; }
    Texture& GetOcclusionTexture() { return textures[1]; }
    Texture& GetEmissiveTexture()  { return textures[2]; }
    
    static inline float16 MakeFloat16(float x) { return (float16)(x * 400.0f); }

#endif

} AMaterial;

typedef struct AImage_
{
    char* path;
} AImage;

typedef struct ANode_
{
    // Warning! order is important, we copy memory in assetmanager.cpp from fbx transform(pos,rot, scale)
    float translation[3];
    float rotation[4];
    float scale[3];

    int   type;  // 0 mesh or 1 camera
    int   index; // index of mesh or camera, -1 if node doesn't have mesh or camera
    int   skin; 
    int   numChildren;
    char* name;
    int*  children;
} ANode;

typedef struct AMorphTarget_
{
    AAttribType attributes;
    // accessor indices of: position, texcoord, tangent respectively
    unsigned short indexes[4]; 
} AMorphTarget;

typedef struct APrimitive_
{
    // pointers to binary file to lookup position, texture, normal..
    void* indices; 
    void* vertices;
    
    unsigned attributes; // AAttribType Position, Normal, TexCoord, Tangent, masks
    int indexType; // GraphicType_UnsignedInt, GraphicType_UnsignedShort.. 
    int numIndices;
    int numVertices;
    int indexOffset;
    short jointType;   // GraphicType_UnsignedInt, GraphicType_UnsignedShort.. 
    short jointCount;  // per vertex bone count (joint), 1-4
    short jointStride; // lets say index data is rgba16u  [r, g, b, a, .......] stride might be bigger than joint

    // weight count is equal to joint count
    short weightType;   // GraphicType_UnsignedInt, GraphicType_UnsignedShort.. 
    short weightStride; // lets say index data is rgba16u  [r, g, b, a, .......] stride might be bigger than joint
    
    // internal use only. after parsing this is useless
    short indiceIndex; // indice index to accessor
    short material;    // material index
    short mode;        // 4 is triangle

    // when we are parsing we use this as an indicator to accessor.
    // after parsing, this will become vertex pointers AAttribType_Position, AAttribType_TexCoord...
    // positions = (Vector3f*)vertexAttribs[0];
    // texCoords = (Vector2f*)vertexAttribs[1]; // note that tangent is vec4
    // ...
    void* vertexAttribs[AAttribType_Count]; 
    // AABB min and max
    float min[4];
    float max[4];

    AMorphTarget* morphTargets; // num morph targets is equal to mesh.num numMorphWeights
} APrimitive;

typedef struct AMesh_
{
    char* name;  
    APrimitive* primitives;
    int numPrimitives;
    int numMorphWeights;
    float* morphWeights;
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


typedef struct GLTFBuffer_
{
    void* uri;
    int byteLength;
} GLTFBuffer;

typedef struct ASkin_
{
    int skeleton; // < index of the root node
    int numJoints;
    float *inverseBindMatrices; // matrix4x4 * numJoints
    int   *joints; // < indices of bone nodes 
    char  *name;
} ASkin;

enum AAnimTargetPath_ {
    AAnimTargetPath_Translation, 
    AAnimTargetPath_Rotation, 
    AAnimTargetPath_Scale,
    AAnimTargetPath_Weight,
    AAnimTargetPath_Make32Bit = 1 << 31
};
typedef int AAnimTargetPath;

enum ASamplerInterpolation_ {
    ASamplerInterpolation_Linear, 
    ASamplerInterpolation_Step, 
    ASamplerInterpolation_CubicSpline
};
typedef int ASamplerInterpolation;

typedef struct AAnimChannel_
{
    int sampler;
    int targetNode;
    AAnimTargetPath targetPath; 
    int pad;
} AAnimChannel;

typedef struct AAnimSampler_
{
    float* input;
    float* output; // quaternion or vector array
    int count;
    int numComponent; // 3,4 vec3 or vec4
    ASamplerInterpolation interpolation;
} AAnimSampler;

typedef struct AAnimation_
{
    int numSamplers;
    int numChannels;
    float duration; // total duration
    float speed; // 1.0 by default, most of the time 1.0
    AAnimChannel* channels;
    AAnimSampler* samplers;
    char* name;
} AAnimation;

typedef struct SceneBundle_
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
    short numBuffers;
    short numAnimations;
    short numSkins;

    AErrorType error;

    void* stringAllocator;
    void* intAllocator;
    void* allVertices;
    void* allIndices;

    int totalVertices;
    int totalIndices;
    float scale;

    GLTFBuffer* buffers;

    AMesh *     meshes;
    ANode*      nodes;
    AMaterial*  materials;
    ATexture*   textures;
    AImage*     images;
    ASampler*   samplers;
    ACamera*    cameras;
    AScene*     scenes;
    AAnimation* animations;
    ASkin*      skins;
} SceneBundle;

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

// outScene should not be null
int ParseGLTF(const char* path, SceneBundle* outScene, float scale);

int ParseObj(const char* path, ParsedObj* scene);


void FreeGLTFBuffers(SceneBundle* gltf);

void FreeSceneBundle(SceneBundle* gltf);

const char* ParsedSceneGetError(AErrorType error);

#endif // ASTL_GLTF_PARSER