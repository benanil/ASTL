
/*****************************************************************
*                                                                *
*    Purpose:                                                    *
*        Simple and efficient parser for GLTF format             *
*        allows you to import 3d mesh, material and scene        *
*    Author:                                                     *
*        Anilcan Gulkaya 2023 anilcangulkaya7@gmail.com          *
*    Restrictions:                                               *
*        No extension support.                                   *
*    License:                                                    *
*        No License whatsoever do whatever you want.             *
*                                                                *
*****************************************************************/

#include "GLTFParser.hpp"

#include "../IO.hpp"
#include "../Array.hpp"
#include "../Math/Matrix.hpp"

#define __private static
#define __public 

struct GLTFAccessor
{
    int bufferView;
    int componentType; // GraphicType
    int count;
    int byteOffset;
    int type; // 1 = SCALAR, 2 = VEC2, 3 = VEC3, 4 = VEC4, mat4
};

struct GLTFBufferView
{
    int buffer;
    int byteOffset;
    int byteLength;
    int target;
    int byteStride;
};

typedef FixedSizeGrowableAllocator<char> AStringAllocator;
typedef FixedSizeGrowableAllocator<int>  AIntAllocator;

#define StrCMP16(_str, _otr) (sizeof(_otr) <= 9 ? StrCmp8(_str, _otr, sizeof(_otr)) : \
                                                  StrCmp16(_str, _otr, sizeof(_otr))) 
 
bool StrCmp8(const char* RESTRICT a, const char* b, uint64_t n)
{
    uint64_t al, bl;
    uint64_t mask = ~0ull >> (64 - ((n-1) * 8));
    al = UnalignedLoad64(a);
    bl = UnalignedLoad64(b);
    return ((al ^ bl) & mask) == 0;
}

bool StrCmp16(const char* RESTRICT a, const char* b, uint64_t bSize)
{
    bool equal = StrCmp8(a, b, 9);
    equal &= StrCmp8(a + 8, b + 8, bSize - 8);
    return equal;
}

inline const char* SkipUntill(const char* curr, char character)
{
    AX_NO_UNROLL while (*curr != character) curr++;
    return curr;
}

inline const char* SkipAfter(const char* curr, char character)
{
    AX_NO_UNROLL while (*curr++ != character);
    return curr;
}

__private const char* CopyStringInQuotes(char*& str, const char* curr, AStringAllocator& stringAllocator)
{
    while (*curr != '"') curr++; // find quote
    curr++; // skip "
    // get length in quotes
    const char* quote = curr;
    while (*quote != '"') quote++;
    int len = (int)(quote - curr);
    char* alloc = stringAllocator.AllocateUninitialized(len + 16);
    str = alloc;
    SmallMemCpy(alloc, curr, len);
    alloc += len;
    curr  += len;
    *alloc = '\0'; // null terminate
    return ++curr;// skip quote
}

__private const char* GetStringInQuotes(char* str, const char* curr)
{
    ++curr; // skip quote " 
    AX_NO_UNROLL while (*curr != '"') *str++ = *curr++;
    *str = '\0'; // null terminate
    return ++curr;
}

__private const char* HashStringInQuotes(uint64_t* hash, const char* curr)
{
    ++curr; // skip quote " 
    uint64_t h = 0ull;
    uint64_t shift = 0;
    AX_NO_UNROLL while (*curr != '"' && shift < 64) { h |= uint64_t(*curr++) << shift; shift += 8; }
    *hash = h;
    return ++curr;
}

// most 8 characters
constexpr uint64_t AHashString8(const char* curr)
{
    uint64_t h = 0ull;
    uint64_t shift = 0;

    while (*curr != '\0') 
    { 
        h |= uint64_t(*curr++) << shift; 
        shift += 8; 
    }
    return h;
}

__private const char* SkipToNextNode(const char* curr, char open, char close)
{
    int balance = 1;
    // skip to first open brackets
    AX_NO_UNROLL while (*curr++ != open);

    while (*curr && balance > 0)
    {
        balance += *curr == open;
        balance -= *curr++ == close;
    }
    return curr;
}

__private const char* ParseFloat16(const char*& curr, short& flt)
{
    flt = (short)(ParseFloat(curr) * 400.0f);
    return curr;
}

__private const char* ParseAccessors(const char* curr, Array<GLTFAccessor>& accessorArray)
{
    GLTFAccessor accessor{};
    curr += 10; // skip accessors"
    // read each accessor
    while (true)
    {
        // search for name
        while (*curr && *curr != '"')
        {
            if (*curr == '}') // next accessor
            {
                accessorArray.Add(accessor);
                MemsetZero(&accessor, sizeof(GLTFAccessor));
            }

            if (*curr == ']') // end all accessors
            return ++curr;
            curr++;
        }
        ASSERTR(*curr != '\0' && "parsing accessors failed probably you forget to close brackets!", return (const char*)AError_CloseBrackets);
        curr++;
        if (StrCMP16(curr, "bufferView"))         accessor.bufferView = ParsePositiveNumber(curr);
        else if (StrCMP16(curr, "byteOffset"))    accessor.byteOffset = ParsePositiveNumber(curr);
        else if (StrCMP16(curr, "componentType")) accessor.componentType = ParsePositiveNumber(curr) - 0x1400; // GL_BYTE 
        else if (StrCMP16(curr, "count"))         accessor.count = ParsePositiveNumber(curr);
        else if (StrCMP16(curr, "name"))
        {
            curr += sizeof("name'"); // we don't need accessor's name
            int numQuotes = 0;
            // skip two quotes
            while (numQuotes < 2)
                numQuotes += *curr++ == '"';
        }
        else if (StrCMP16(curr, "type"))
        {
            curr += sizeof("type'"); // skip type
            curr = SkipUntill(curr, '"');
            uint64_t hash;
            curr = HashStringInQuotes(&hash, curr);
            
            switch (hash)
            {   case AHashString8("SCALAR"): accessor.type = 1; break;
                case AHashString8("VEC2"):   accessor.type = 2; break;
                case AHashString8("VEC3"):   accessor.type = 3; break;
                case AHashString8("VEC4"):   accessor.type = 4; break;
                case AHashString8("MAT4"):   accessor.type = 16;break;
                default: ASSERT(0 && "Unknown accessor type");
            };
        }
        else if (StrCMP16(curr, "min")) curr = SkipToNextNode(curr, '[', ']'); // skip min and max
        else if (StrCMP16(curr, "max")) curr = SkipToNextNode(curr, '[', ']');
        else if (StrCMP16(curr, "normalized")) curr = SkipAfter(curr, '"');
        else
        {
            ASSERT(0 && "unknown accessor var");
            return (const char*)AError_UNKNOWN_ACCESSOR_VAR;
        }
    }
}

__private const char* ParseBufferViews(const char* curr, Array<GLTFBufferView>& bufferViews)
{
    GLTFBufferView bufferView{};
    curr += sizeof("bufferViews'");

    // read each buffer view
    while (true)
    {
        // search for name
        while (*curr && *curr != '"')
        {
            if (*curr == '}') // next buffer view
            {
                bufferViews.Add(bufferView);
                MemsetZero(&bufferView, sizeof(GLTFBufferView));
            }
            if (*curr++ == ']') return curr; // end all buffer views
        }
        ASSERTR(*curr != '0' && "buffer view parse failed, probably you forgot to close brackets!", return (const char*)AError_CloseBrackets);

        uint64_t hash;
        curr = HashStringInQuotes(&hash, curr);

        switch (hash)
        {
            case AHashString8("buffer"):    bufferView.buffer = ParsePositiveNumber(++curr); break; 
            case AHashString8("byteOffs"):  bufferView.byteOffset = ParsePositiveNumber(++curr); break; 
            case AHashString8("byteLeng"):  bufferView.byteLength = ParsePositiveNumber(++curr); break; 
            case AHashString8("byteStri"):  bufferView.byteStride = ParsePositiveNumber(++curr); break; 
            case AHashString8("target"):    bufferView.target = ParsePositiveNumber(++curr); break;
            case AHashString8("name"):  {
                int numQuote = 0;
                while (numQuote < 2)
                    numQuote += *curr++ == '"';
                break;
            }
            default: {
                ASSERT(0 && "UNKNOWN buffer view value!");
                return (const char*)AError_UNKNOWN_BUFFER_VIEW_VAR;
            }
        };
    }
}

static void DecodeBase64(char *dst, const char *src, size_t src_length)
{
    struct Base64Table
    {
        uint8_t table[256] = { 0 };
        constexpr Base64Table()
        {
            for (char c = 'A'; c <= 'Z'; c++) table[c] = (uint8_t)(c - 'A');
            for (char c = 'a'; c <= 'z'; c++) table[c] = (uint8_t)(26 + (c - 'a'));
            for (char c = '0'; c <= '9'; c++) table[c] = (uint8_t)(52 + (c - '0'));
            table['+'] = 62;
            table['/'] = 63;
        }
    };
    
    static constexpr Base64Table table{};
    
    for (uint64_t i = 0; i + 4 <= src_length; i += 4) {
        uint32_t a = table.table[src[i + 0]];
        uint32_t b = table.table[src[i + 1]];
        uint32_t c = table.table[src[i + 2]];
        uint32_t d = table.table[src[i + 3]];
           
        dst[0] = (char)(a << 2 | b >> 4);
        dst[1] = (char)(b << 4 | c >> 2);
        dst[2] = (char)(c << 6 | d);
        dst += 3;
    }
}

__private const char* ParseBuffers(const char* curr, const char* path, Array<GLTFBuffer>& bufferArray)
{
    GLTFBuffer buffer{};
    curr += sizeof("buffers'"); // skip buffers"
    char binFilePath[256]={0};
    char* endOfWorkDir = binFilePath;
    int binPathlen = StringLength(path);
    SmallMemCpy(binFilePath, path, binPathlen);
    endOfWorkDir = binFilePath + binPathlen;
    
    // remove bla.gltf
    while (binFilePath[binPathlen - 1] != '/' && binFilePath[binPathlen - 1] != '\\')
        binFilePath[--binPathlen] = '\0', endOfWorkDir--;

    // read each buffer
    while (true)
    {
        // search for name
        while (*curr && *curr != '"')
        {
            if (*curr == '}') // next buffer
            {
                bufferArray.Add(buffer);
                MemsetZero(&buffer, sizeof(GLTFBuffer));
                MemsetZero(endOfWorkDir, sizeof(binFilePath) - (size_t)(endOfWorkDir - binFilePath));
            }

            if (*curr++ == ']') return curr; // end all buffers
        }
        ASSERTR(*curr && "parsing buffers failed, probably you forgot to close braces", return (const char*)AError_CloseBrackets);
        curr++;
        if (StrCMP16(curr, "uri")) // is uri
        {
            curr += sizeof("uri'"); // skip uri": 
            while (*curr != '"') curr++;
            if (StartsWith(curr, "\"data:"))
            {
                curr = SkipAfter(curr, ',');
                uint64_t base64Size = 0;
                while (curr[base64Size] != '\"') {
                    base64Size++;
                }
                buffer.uri = (void*)new char[base64Size];
                DecodeBase64((char*)buffer.uri, curr, base64Size);
                curr += base64Size + 1;
            }
            else
            {
                curr = GetStringInQuotes(endOfWorkDir, curr);
                buffer.uri = ReadAllFile(binFilePath);
                ASSERT(buffer.uri && "uri is not exist");
                if (!buffer.uri) return (const char*)AError_BIN_NOT_EXIST;
            }
        }
        else if (StrCMP16(curr, "byteLength"))
        {
            buffer.byteLength = ParsePositiveNumber(++curr);
        }
        else
        {
            ASSERT(0 && "Unknown buffer variable! byteLength or uri excepted.");
            return (const char*)AError_BUFFER_PARSE_FAIL;
        }
    }
}

// write paths to path buffer, buffer is seperated by null terminators
__private const char* ParseImages(const char* curr, const char* path, Array<AImage>& images, AStringAllocator& stringAllocator)
{
    curr = SkipUntill(curr, '[');
    curr++;
    
    int pathLen = StringLength(path);
    while (path[pathLen-1] != '/') pathLen--;

    AImage image{};
    // read each buffer
    while (true)
    {
        // search for name
        while (*curr && *curr != '"')
            if (*curr++ == ']')
                return curr; // end all images
        
        ASSERTR(*curr != '\0' && "parse images failed probably you forgot to close brackets", return (const char*)AError_CloseBrackets);

        curr++;
        bool isUri = StrCMP16(curr, "uri");

        // mimeType and name is not supported
        if (isUri)
        {
            curr += 4;
            curr = SkipAfter(curr, '"');
            int uriSize = 0;
            while (curr[uriSize] != '"') uriSize++;

            image.path = stringAllocator.Allocate(uriSize + pathLen + 16);
            SmallMemCpy(image.path, path, pathLen);
            SmallMemCpy(image.path + pathLen, curr, uriSize);
            image.path[uriSize + pathLen] = '\0';
            curr += uriSize;
            images.PushBack(image);
        }
    }
    return nullptr;
}

__private const char* ParseTextures(const char* curr, Array<ATexture>& textures, AStringAllocator& stringAllocator)
{
    curr += sizeof("textures'");
    ATexture texture{};

    // read each buffer
    while (true)
    {
        // search for name
        while (*curr && *curr != '"')
        {
            if (*curr == '}') // next buffer view
            {
                textures.Add(texture);
                MemsetZero(&texture, sizeof(ATexture));
            }

            if (*curr == ']') // end all buffer views
            return ++curr;
            curr++;
        }
        ASSERTR(*curr != '\0' && "parse images failed probably you forgot to close brackets", return (const char*)AError_CloseBrackets);
        curr++;
        if (StrCMP16(curr, "sampler")) 
        {
            texture.sampler = ParsePositiveNumber(++curr);
        }
        else if (StrCMP16(curr, "source"))
        {
            texture.source = ParsePositiveNumber(++curr);
        }
        else if (StrCMP16(curr, "name"))
        {
            curr = CopyStringInQuotes(texture.name, curr + 5, stringAllocator);
        }
        else {
            ASSERT(0 && "Unknown buffer variable! sampler, source or name excepted.");
            return (const char*)AError_UNKNOWN_TEXTURE_VAR;
        }
    }
}

__private const char* ParseAttributes(const char* curr, APrimitive* primitive)
{
    curr += sizeof("attributes'");

    while (true)
    {
        while (*curr != '"')
        if (*curr++ == '}') return curr;
        
        curr++; // skip "
        unsigned maskBefore = primitive->attributes;
        if      (StrCMP16(curr, "POSITION"))   { primitive->attributes |= AAttribType_POSITION;   curr += sizeof("POSITION'");   }
        else if (StrCMP16(curr, "NORMAL"))     { primitive->attributes |= AAttribType_NORMAL;     curr += sizeof("NORMAL'");     }
        else if (StrCMP16(curr, "TEXCOORD_0")) { primitive->attributes |= AAttribType_TEXCOORD_0; curr += sizeof("TEXCOORD_0'"); }
        else if (StrCMP16(curr, "TANGENT"))    { primitive->attributes |= AAttribType_TANGENT;    curr += sizeof("TANGENT'");    }
        else if (StrCMP16(curr, "TEXCOORD_1")) { primitive->attributes |= AAttribType_TEXCOORD_1; curr += sizeof("TEXCOORD_1'"); }
        else if (StrCMP16(curr, "JOINTS_0"))   { primitive->attributes |= AAttribType_JOINTS;     curr += sizeof("JOINTS_0'");   }
        else if (StrCMP16(curr, "WEIGHTS_0"))  { primitive->attributes |= AAttribType_WEIGHTS;    curr += sizeof("WEIGHTS_0'");  }
        else if (StrCMP16(curr, "TEXCOORD_"))  { curr += sizeof("TEXCOORD_X'"); continue; } // < NO more than two texture coords
        else { ASSERT(0 && "attribute variable unknown!"); return (const char*)AError_UNKNOWN_ATTRIB; }

        // using bitmask will help us to order attributes correctly(sort) Position, Normal, TexCoord
        unsigned newIndex = TrailingZeroCount32(maskBefore ^ primitive->attributes);
        primitive->vertexAttribs[newIndex] = (void*)(uint64_t)ParsePositiveNumber(curr);
    }
}

__private const char* ParseMeshes(const char* curr, Array<AMesh>& meshes, AStringAllocator& stringAllocator)
{
    char text[64]{};
    curr += sizeof("meshes'"); // skip meshes" 
    AMesh mesh{};
    MemsetZero(&mesh, sizeof(AMesh));

    // parse all meshes
    while (true)
    {
        while (*curr != '"')
        {
            if (*curr == '}') 
            {
                meshes.Add(mesh);
                MemsetZero(&mesh, sizeof(AMesh));
            }
            if (*curr++ == ']') return curr; // end of meshes
        }
        curr = GetStringInQuotes(text, curr);
        
        if (StrCMP16(text, "name")) {
            curr = CopyStringInQuotes(mesh.name, curr, stringAllocator); 
            continue; 
        }
        else if (StrCMP16(text, "weights"))
        {
            curr = SkipAfter(curr, '[');
            const char* begin = curr;
            
            int numWeights = 1;
            while (*curr != ']')
                numWeights += *curr++ == ',';
            
            curr = begin;
            mesh.numMorphWeights = numWeights;
            mesh.morphWeights = new float[numWeights];
            for (int i = 0; i < numWeights; i++)
            {
                mesh.morphWeights[i] = ParseFloat(curr);
            }
            curr = SkipAfter(curr, ']');
            continue;
        }
        else if (!StrCMP16(text, "primitives")) { 
            ASSERT(0 && "only primitives, name and weights allowed"); 
            return (const char*)AError_UNKNOWN_MESH_VAR; 
        }

        APrimitive primitive{};  
        primitive.material = -1;
        // parse primitives
        while (true)
        {
            while (*curr != '"')
            {
                if (*curr == '}')
                {
                    SBPush(mesh.primitives, primitive);
                    MemsetZero(&primitive, sizeof(APrimitive));
                    mesh.numPrimitives++;
                    primitive.material = -1;
                }

                if (*curr++ == ']') goto end_primitives; // this is end of primitive list
            }
            curr++;
            
            if      (StrCMP16(curr, "attributes")) { curr = ParseAttributes(curr, &primitive); }
            else if (StrCMP16(curr, "indices"))    { primitive.indiceIndex = ParsePositiveNumber(curr); }
            else if (StrCMP16(curr, "mode"))       { primitive.mode        = ParsePositiveNumber(curr); }
            else if (StrCMP16(curr, "material"))   { primitive.material    = ParsePositiveNumber(curr); }
            else if (StrCMP16(curr, "targets"))    
            {
                curr += sizeof("targets'");
                AMorphTarget morphTarget = {};
                while (*curr)
                {
                    while (*curr != '"')
                    {
                        if (*curr == '}') {
                            SBPush(primitive.morphTargets, morphTarget);
                            MemsetZero(&morphTarget, sizeof(AMorphTarget));
                        }
                        if (*curr++ == ']') goto end_morphs;
                    }
                    curr++; // skip "

                    unsigned maskBefore = morphTarget.attributes;
                    if      (StrCMP16(curr, "POSITION"))   { morphTarget.attributes |= AAttribType_POSITION;   curr += sizeof("POSITION'"); }
                    else if (StrCMP16(curr, "TEXCOORD_0")) { morphTarget.attributes |= AAttribType_TEXCOORD_0; curr += sizeof("TEXCOORD_0'"); }
                    else if (StrCMP16(curr, "NORMAL"))     { morphTarget.attributes |= AAttribType_NORMAL;     curr += sizeof("NORMAL'"); }
                    else if (StrCMP16(curr, "TANGENT"))    { morphTarget.attributes |= AAttribType_TANGENT;    curr += sizeof("TANGENT'"); }
                    else if (StrCMP16(curr, "TEXCOORD_"))  { curr = SkipAfter(curr, '"'); continue; } // < NO more than one texture coords
                    else { ASSERT(0 && "attribute variable unknown!"); return (const char*)AError_UNKNOWN_ATTRIB; }
                 
                    // detect changed attribute.
                    unsigned addedAttribute = TrailingZeroCount32(maskBefore ^ morphTarget.attributes);
                    morphTarget.indexes[addedAttribute] = (short)ParsePositiveNumber(curr);
                }
                end_morphs:{}
            }
            else { ASSERT(0); return (const char*)AError_UNKNOWN_MESH_PRIMITIVE_VAR; }
        }
        end_primitives:{}
        curr++; // skip ]
    }
    return nullptr;
}

struct IntPtrPair { int numElements, *ptr; };

static IntPtrPair ParseIntArray(const char*& cr, AIntAllocator& intAllocator)
{
    const char* curr = cr;
    IntPtrPair result={};
    // find how many elements there are:
    while (!IsNumber(*curr)) curr++;
    
    const char* begin = curr;
    result.numElements = 1;
    while (true)
    {
        result.numElements += *curr == ',';
        if (*curr++ == ']') break;
    }
    curr = begin;
    result.ptr = intAllocator.AllocateUninitialized(result.numElements);
    result.numElements = 0;
    
    while (*curr != ']')
    {
        if (IsNumber(*curr))
        {
            result.ptr[result.numElements] = ParsePositiveNumber(curr);
            result.numElements++;
        }
        curr++;
    }
    cr = curr;
    return result;
}

__private const char* ParseNodes(const char* curr,
                                 Array<ANode>& nodes,
                                 AStringAllocator& stringAllocator,
                                 FixedSizeGrowableAllocator<int>& intAllocator, float scale)
{
    curr = SkipUntill(curr, '[');
    curr++;
    ANode node{};
    node.rotation[3] = 1.0f;
    node.scale[0] = node.scale[1] = node.scale[2] = scale; 
    node.index = -1;
    
    // read each node
    while (true)
    {
        // search for name
        while (*curr && *curr != '"')
        {
            if (*curr == '}')
            {
                nodes.Add(node);
                MemsetZero(&node, sizeof(ANode));
                node.rotation[3] = 1.0f;
                node.scale[0] = node.scale[1] = node.scale[2] = scale;
                node.index = -1;
            }
            if (*curr++ == ']') return curr; // end all nodes
        }
        ASSERTR(*curr != '\0' && "parsing nodes not possible, probably forgot to close brackets!", return (const char*)AError_CloseBrackets);
        curr++; // skips the "
        
        // mesh, name, children, matrix, translation, rotation, scale, skin
        if      (StrCMP16(curr, "mesh"))   { node.type = 0; node.index = ParsePositiveNumber(curr); continue; } // don't want to skip ] that's why continue
        else if (StrCMP16(curr, "camera")) { node.type = 1; node.index = ParsePositiveNumber(curr); continue; } // don't want to skip ] that's why continue
        else if (StrCMP16(curr, "children"))
        {
            IntPtrPair result = ParseIntArray(curr, intAllocator);
            node.numChildren = result.numElements;
            node.children = result.ptr;
        }
        else if (StrCMP16(curr, "matrix"))
        {
            Matrix4 m;
            float* matrix = m.GetPtr();
            
            for (int i = 0; i < 16; i++)
                matrix[i] = ParseFloat(curr);
            
            m = Matrix4::Transpose(m);
            node.translation[0] = matrix[12];
            node.translation[1] = matrix[13];
            node.translation[2] = matrix[14];
            QuaternionFromMatrix(node.rotation, matrix);

            vec_t v = VecMulf(Matrix4::ExtractScaleV(m), scale);
            Vec3Store(node.scale, v);
        }
        else if (StrCMP16(curr, "translation"))
        {
            node.translation[0] = ParseFloat(curr);
            node.translation[1] = ParseFloat(curr);
            node.translation[2] = ParseFloat(curr);
        }
        else if (StrCMP16(curr, "rotation"))
        {
            node.rotation[0] = ParseFloat(curr);
            node.rotation[1] = ParseFloat(curr);
            node.rotation[2] = ParseFloat(curr);
            node.rotation[3] = ParseFloat(curr);
        }
        else if (StrCMP16(curr, "scale"))
        {
            node.scale[0] = ParseFloat(curr) * scale;
            node.scale[1] = ParseFloat(curr) * scale;
            node.scale[2] = ParseFloat(curr) * scale;
        }
        else if (StrCMP16(curr, "name"))
        {
            curr = CopyStringInQuotes(node.name, curr + 5, stringAllocator);
            continue; 
        }
        else if (StrCMP16(curr, "skin"))
        {
            node.skin = ParsePositiveNumber(curr);
            continue; // continue because we don't want to skip ] and it is not exist
        }
        else
        {
            ASSERT(0 && "Unknown node variable");
            return (const char*)AError_UNKNOWN_NODE_VAR;
        }

        curr = SkipUntill(curr, ']');
        curr++;
    }
    return nullptr;
}

__private const char* ParseCameras(const char* curr, Array<ACamera>& cameras, AStringAllocator& stringAllocator)
{
    curr += sizeof("camera'");
    char text[64]{};
    ACamera camera{};
    // parse all meshes
    while (true)
    {
        while (*curr != '"')
        {
            if (*curr == '}') 
            {
                cameras.Add(camera);
                MemsetZero(&camera, sizeof(ACamera));
            }
            if (*curr++ == ']') return curr; // end of cameras
        }
        curr = GetStringInQuotes(text, curr);
        
        if (StrCMP16(text, "name")) {
            curr = CopyStringInQuotes(camera.name, curr, stringAllocator); 
            continue; 
        }
        if (StrCMP16(text, "type")) {
            curr = SkipUntill(curr, '"');
            curr++;
            camera.type = *curr == 'p'; // 0 orthographic 1 perspective 
            curr = SkipUntill(curr, '"');
            curr++;
            continue; 
        }
        else if (!StrCMP16(text, "orthographic") && !StrCMP16(text, "perspective")) { 
            ASSERT(0 && "unknown camera variable"); 
            return (const char*)AError_UNKNOWN_CAMERA_VAR; 
        }

        // parse primitives
        while (true)
        {
            while (*curr != '"')
            if (*curr++ == '}')  goto end_properties; // this is end of camera variables
            
            curr++;
            if      (StrCMP16(curr, "zfar"))        { camera.zFar        = ParseFloat(curr); }
            else if (StrCMP16(curr, "znear"))       { camera.zNear       = ParseFloat(curr); }
            else if (StrCMP16(curr, "aspectRatio")) { camera.aspectRatio = ParseFloat(curr); }
            else if (StrCMP16(curr, "yfov"))        { camera.yFov        = ParseFloat(curr); }
            else if (StrCMP16(curr, "xmag"))        { camera.xmag        = ParseFloat(curr); }
            else if (StrCMP16(curr, "ymag"))        { camera.ymag        = ParseFloat(curr); }
            else { ASSERT(0); return (const char*)AError_UNKNOWN_CAMERA_VAR; }
        }
        end_properties:{}
    }
    return nullptr;
}

__private const char* ParseScenes(const char* curr, Array<AScene>& scenes, 
                                  AStringAllocator& stringAllocator, FixedSizeGrowableAllocator<int>& intAllocator)
{
    curr = SkipUntill(curr, '[');
    curr++;
    AScene scene{};
    // read each node
    while (true)
    {
        // search for name
        while (*curr && *curr != '"')
        {
            if (*curr == '}')
            {
                scenes.Add(scene);
                MemsetZero(&scene, sizeof(AScene));
            }
            if (*curr++ == ']') return curr; // end all scenes
        }
        ASSERTR(*curr != '\0' && "parsing scenes not possible, probably forgot to close brackets!", return (const char*)AError_CloseBrackets);
        curr++; // skips the "
        
        if (StrCMP16(curr, "nodes"))
        {
            // find how many childs there are:
            while (!IsNumber(*curr)) curr++;
            const char* begin = curr;
            scene.numNodes = 1;
            while (true)
            {
                scene.numNodes += *curr == ',';
                if (*curr++ == ']') break;
            }
            curr = begin;
            scene.nodes = intAllocator.AllocateUninitialized(scene.numNodes);
            scene.numNodes = 0;

            while (*curr != ']')
            {
                if (IsNumber(*curr))
                {
                    scene.nodes[scene.numNodes] = ParsePositiveNumber(curr);
                    scene.numNodes++;
                }
                curr++;
            }
            curr++;// skip ]
        }
        else if (StrCMP16(curr, "name"))
        {
            curr = CopyStringInQuotes(scene.name, curr + 5, stringAllocator);
        }
    } 
}
inline char OGLWrapToWrap(int wrap)
{
    switch (wrap)
    {
        case 0x2901: return 0; // GL_REPEAT          10497
        case 0x812F: return 1; // GL_CLAMP_TO_EDGE   33071
        case 0x812D: return 2; // GL_CLAMP_TO_BORDER 33069
        case 0x8370: return 3; // GL_MIRRORED_REPEAT 33648
        default: ASSERT(0 && "wrong or undefined sampler type!"); return 0;
    }
}

__private const char* ParseSamplers(const char* curr, Array<ASampler>& samplers)
{
    curr = SkipUntill(curr, '[');
    curr++;
    
    ASampler sampler{};
    
    // read each node
    while (true)
    {
        // search for name
        while (*curr && *curr != '"')
        {
            if (*curr == '}')
            {
                samplers.Add(sampler);
                MemsetZero(&sampler, sizeof(ASampler));
            }
            if (*curr++ == ']') return curr; // end all nodes
        }
        ASSERTR(*curr != '\0' && "parsing nodes not possible, probably forgot to close brackets!", return (const char*)AError_CloseBrackets);
        curr++; // skips the "

        if      (StrCMP16(curr, "magFilter")) sampler.magFilter = (char)(ParsePositiveNumber(curr) - 0x2600); // GL_NEAREST 9728, GL_LINEAR 0x2601 9729
        else if (StrCMP16(curr, "minFilter")) sampler.minFilter = (char)(ParsePositiveNumber(curr) - 0x2600); // GL_NEAREST 9728, GL_LINEAR 0x2601 9729
        else if (StrCMP16(curr, "wrapS"))     sampler.wrapS = (char)OGLWrapToWrap(ParsePositiveNumber(curr));
        else if (StrCMP16(curr, "wrapT"))     sampler.wrapT = (char)OGLWrapToWrap(ParsePositiveNumber(curr));
        else { ASSERT(0 && "parse samplers failed!"); return (const char*)AError_UNKNOWN; }
    }
}

__private const char* ParseMaterialTexture(const char* curr, AMaterial::Texture& texture)
{
    curr = SkipUntill(curr, '{');
    curr++;
    texture.strength = AMaterial::MakeFloat16(1.0f);
    while (true)
    {
        while (*curr && *curr != '"')
        {
            if (*curr++ == '}') return curr;
        }
        ASSERTR(*curr && "parsing material failed, probably forgot to close brackets", return (const char*)AError_CloseBrackets);
        curr++;

        if (StrCMP16(curr, "scale")) {
            curr = ParseFloat16(curr, texture.scale);
        }
        else if (StrCMP16(curr, "index")) {
            texture.index = (char)ParsePositiveNumber(curr);
        }
        else if (StrCMP16(curr, "texCoord")) {
            texture.texCoord = (char)ParsePositiveNumber(curr);
        }
        else if (StrCMP16(curr, "strength")) {
            curr = ParseFloat16(curr, texture.strength);
        }
        else if (StrCMP16(curr, "extensions")) {
            curr = SkipToNextNode(curr, '{', '}'); // currently extensions are not supported 
        }
        else {
            ASSERT(0 && "unknown material texture value");
            return (const char*)AError_UNKNOWN_MATERIAL_VAR;
        }
    }
    return nullptr;
}

__private const char* ParseMaterials(const char* curr, Array<AMaterial>& materials, AStringAllocator& stringAllocator)
{
    // mesh, name, children, 
    // matrix, translation, rotation, scale
    curr = SkipUntill(curr, '[');
    curr++; // skip '['
    AMaterial material{};
    // read each material
    while (true)
    {
        // search for name
        while (*curr && *curr != '"')
        {
            if (*curr == '}')
            {
                materials.Add(material);
                MemsetZero(&material, sizeof(AMaterial));
                material.baseColorTexture.index = -1;
                material.metallicFactor = PackUnorm16(1.0f);
                material.roughnessFactor = PackUnorm16(1.0f);
            }
            if (*curr++ == ']') return curr; // end all nodes
        }
        ASSERTR(*curr && "parsing material failed, probably forgot to close brackets", return (const char*)AError_CloseBrackets);

        int texture = -1;
        curr++; // skips the "
        if (StrCMP16(curr, "name"))
        {
            curr = CopyStringInQuotes(material.name, curr + 5, stringAllocator);
        }
        else if (StrCMP16(curr, "doubleSided"))
        {
            curr += sizeof("doubleSided'"); // skip doubleSided"
            AX_NO_UNROLL while (!IsLower(*curr)) curr++;
            material.doubleSided = *curr == 't';
        }
        else if (StrCMP16(curr, "pbrMetallicRoug")) //pbrMetallicRoughhness
        {
            curr = SkipUntill(curr, '{');
    
            // parse until finish
            while (true)
            {
                // search for name
                while (*curr && *curr != '"')
                {
                    if (*curr++ == '}') { goto pbr_end; }
                }
                curr++; // skips the "
                
                if      (StrCMP16(curr, "baseColorTex"))  { curr = ParseMaterialTexture(curr, material.baseColorTexture); }
                else if (StrCMP16(curr, "metallicRough")) { curr = ParseMaterialTexture(curr, material.metallicRoughnessTexture); }
                else if (StrCMP16(curr, "baseColorFact"))
                {
                    float baseColorFactor[4] = { ParseFloat(curr), ParseFloat(curr), ParseFloat(curr), ParseFloat(curr)};
                    material.baseColorFactor = PackColorRGBAU32(baseColorFactor);
                    curr = SkipUntill(curr, ']');
                    curr++;
                }
                else if (StrCMP16(curr, "metallicFact"))
                {
                    material.metallicFactor = PackUnorm16(ParseFloat(curr));
                }
                else if (StrCMP16(curr, "roughnessFact"))
                {
                    material.roughnessFactor = PackUnorm16(ParseFloat(curr));
                }
                else
                {
                    ASSERT(0 && "unknown pbrMetallicRoughness value!");
                    return (char*)AError_UNKNOWN_PBR_VAR;
                }
            }
            pbr_end: {}
        }
        else if (StrCMP16(curr, "normalTexture"))    texture = 0;
        else if (StrCMP16(curr, "occlusionTextur"))  texture = 1;
        else if (StrCMP16(curr, "emissiveTexture"))  texture = 2;
        else if (StrCMP16(curr, "emissiveFactor")) 
        {
            curr = ParseFloat16(curr, material.emissiveFactor[0]); 
            curr = ParseFloat16(curr, material.emissiveFactor[1]);
            curr = ParseFloat16(curr, material.emissiveFactor[2]);
            curr = SkipUntill(curr, ']');
            curr++;
        }
        else if (StrCMP16(curr, "extensions"))
        {
            curr = SkipToNextNode(curr, '{', '}'); // currently extensions are not supported 
        }
        else if (StrCMP16(curr, "alphaMode"))
        {
            char text[16]={0};
            curr += sizeof("alphaMode'");
            curr = SkipUntill(curr, '"');
            curr = GetStringInQuotes(text, curr);
                 if (StrCMP16(text, "OPAQUE")) material.alphaMode = AMaterialAlphaMode_Opaque;
            else if (StrCMP16(text, "MASK"))   material.alphaMode = AMaterialAlphaMode_Mask;
            else if (StrCMP16(text, "BLEND"))  material.alphaMode = AMaterialAlphaMode_Blend;
        }
        else if (StrCMP16(curr, "alphaCutoff"))
        {
            material.alphaCutoff = ParseFloat(curr);
        }
        else if (StrCMP16(curr, "extras"))
        {
            curr = SkipAfter(curr, '{');
            int balance = 1;
            while (balance > 0) {
                balance += *curr == '{';
                balance -= *curr == '}';
                curr++;
            }
        }
        else
        {
            ASSERTR(0 && "undefined material variable!", return (const char*)AError_UNKNOWN_MATERIAL_VAR);
        }

        if (texture != -1)
        {
            curr = ParseMaterialTexture(curr, material.textures[texture]);
            if ((uint64_t)curr == AError_UNKNOWN_MATERIAL_VAR) return curr;
        }
    }
    return curr;
}

static const char* ParseSkins(const char* curr, Array<ASkin>& skins, AStringAllocator& stringAllocator, AIntAllocator& intAllocator)
{
    curr = SkipAfter(curr, '[');
    
    ASkin skin{};
    skin.skeleton = -1;

    // read each node
    while (true)
    {
        // search for name
        while (*curr && *curr != '"')
        {
            if (*curr == '}')
            {
                skins.Add(skin);
                MemsetZero(&skin, sizeof(ASkin));
                skin.skeleton = -1;
            }
            if (*curr++ == ']') return curr; // end all nodes
        }
        ASSERTR(*curr != '\0' && "parsing skins not possible, probably forgot to close brackets!", return (const char*)AError_CloseBrackets);
        curr++; // skips the "

        if (StrCMP16(curr, "inverseBindMatrices"))
        {
            // we will parse later, because we are not sure we are parsed accessors at this point
            skin.inverseBindMatrices = (float*)(size_t)ParsePositiveNumber(curr);
        }
        else if (StrCMP16(curr, "skeleton"))            skin.skeleton = ParsePositiveNumber(curr);
        else if (StrCMP16(curr, "name"))                curr = CopyStringInQuotes(skin.name, curr + 5, stringAllocator);
        else if (StrCMP16(curr, "joints"))
        {
            IntPtrPair result = ParseIntArray(curr, intAllocator);
            skin.numJoints = result.numElements;
            skin.joints = result.ptr;
            curr++; // skip ]
        }
    }
    return curr;
}

static const char* ParseAnimations(const char* curr, Array<AAnimation>& animations, 
                                   AStringAllocator& stringAllocator, AIntAllocator& intAllocator)
{
    curr = SkipAfter(curr, '[');
    Array<AAnimChannel> channels;
    Array<AAnimSampler> samplers;

    AAnimation animation{};
    animation.speed = 1.0f;
    // read each node
    while (true)
    {
        // search for name
        while (*curr && *curr != '"')
        {
            if (*curr == '}')
            {
                animation.numSamplers = samplers.Size();
                animation.numChannels = channels.Size();
                animation.samplers = samplers.TakeOwnership();
                animation.channels = channels.TakeOwnership();
                animations.Add(animation);
                MemsetZero(&animation, sizeof(AAnimation));
                animation.speed = 1.0f;
            }
            if (*curr++ == ']') 
                return curr; // end all nodes
        }
        ASSERTR(*curr != '\0' && "parsing animations not possible, probably forgot to close brackets!", return (const char*)AError_CloseBrackets);
        curr++; // skips the "

        if (StrCMP16(curr, "name"))
        {
            curr = CopyStringInQuotes(animation.name, curr + sizeof("name'"), stringAllocator);
        }
        else if (StrCMP16(curr, "channels"))
        {
            curr += sizeof("channels'");
            AAnimChannel channel;
            bool parsingTarget = false;
            while (true)
            {
                while (*curr && *curr != '"')
                {
                    if (*curr == ']') { curr++; /* skip ] */ goto end_parsing; }
                    if (*curr == '}') 
                    {
                        if (parsingTarget) parsingTarget = false;
                        else 
                        {
                            channels.Add(channel);
                            MemsetZero(&channel, sizeof(AAnimChannel));
                        }
                    }
                    curr++;
                }
                ASSERTR(*curr != '\0' && "parsing anim channels not possible, probably forgot to close brackets!", return (const char*)AError_CloseBrackets);

                uint64_t hash;
                curr = HashStringInQuotes(&hash, curr);

                switch (hash)
                {
                    case AHashString8("sampler"): channel.sampler = ParsePositiveNumber(curr);     break;
                    case AHashString8("node"):    channel.targetNode = ParsePositiveNumber(curr);  break;
                    case AHashString8("target"):  curr += sizeof("target'"); parsingTarget = true; break; 
                    case AHashString8("path"):
                    {
                        curr = SkipAfter(curr, '"');
                        switch (*curr) {
                            case 't': channel.targetPath = AAnimTargetPath_Translation; curr += sizeof("translation'"); break;
                            case 'r': channel.targetPath = AAnimTargetPath_Rotation;    curr += sizeof("rotation'"); break;
                            case 's': channel.targetPath = AAnimTargetPath_Scale;       curr += sizeof("scale'"); break;
                            case 'w': channel.targetPath = AAnimTargetPath_Weight;      curr += sizeof("weights'"); break;
                            default: ASSERT(0 && "Unknown animation path value");
                        };
                        break;
                    }
                    default: ASSERT(0 && "Unknown animation channel value");
                };
            }
        }
        else if (StrCMP16(curr, "samplers"))
        {
            curr += sizeof("samplers'");
            AAnimSampler sampler;
            while (true)
            {
                while (*curr && *curr != '"')
                {
                    if (*curr == ']') { curr++; /* skip ] */ goto end_parsing; }
                    if (*curr == '}') 
                    {
                        samplers.Add(sampler);
                        MemsetZero(&sampler, sizeof(AAnimSampler));
                    }
                    curr++;
                }
                ASSERTR(*curr != '\0' && "parsing anim channels not possible, probably forgot to close brackets!", return (const char*)AError_CloseBrackets);

                uint64_t hash;
                curr = HashStringInQuotes(&hash, curr);

                switch (hash)
                {
                    case AHashString8("input"):         sampler.input  = (float*)(size_t)ParsePositiveNumber(curr); break;
                    case AHashString8("output"):        sampler.output = (float*)(size_t)ParsePositiveNumber(curr); break;
                    case AHashString8("interpol"): // you've been searching from interpol hands up!!
                    {
                        curr += sizeof("interpolation") - sizeof("interpol");
                        curr = SkipAfter(curr, '"');
                        switch (*curr)
                        {
                            case 'L': sampler.interpolation = 0; curr += sizeof("Linear'");      break; // Linear
                            case 'S': sampler.interpolation = 1; curr += sizeof("Step'");        break; // Step
                            case 'C': sampler.interpolation = 2; curr += sizeof("CubicSpline'"); break; // CubicSpline
                            default: ASSERT(0 && "Unknown animation path value"); break;
                        };
                        break;
                    }
                    default: ASSERT(0 && "Unknown animation sampler value"); break;
                };
            }
        }
        end_parsing:{}
    }
    return curr;
}

__public int ParseGLTF(const char* path, SceneBundle* result, float scale)
{
    ASSERT(result && path);
    uint64_t sourceSize = 0;
    char* source = ReadAllFile(path, nullptr);
    MemsetZero(result, sizeof(SceneBundle));

    if (source == nullptr) { result->error = AError_FILE_NOT_FOUND; ASSERT(0); return 0; }

#if defined(DEBUG) || defined(_DEBUG)
    // ascii utf8 support check
    // if (IsUTF8ASCII(source, sourceSize) != 1) { result->error = AError_NON_UTF8; return; }
#endif
    Array<GLTFBufferView> bufferViews ;
    Array<GLTFBuffer>     buffers     ;
    Array<GLTFAccessor>   accessors   ;

    AStringAllocator stringAllocator(2048);
    FixedSizeGrowableAllocator<int> intAllocator(512);

    Array<AMesh>  meshes; Array<ANode>        nodes; Array<AMaterial> materials; Array<ATexture>  textures;
    Array<AImage> images; Array<ASampler>  samplers; Array<ACamera>     cameras; Array<AScene>    scenes;
    Array<ASkin>  skins; Array<AAnimation> animations;

    const char* curr = source;
    while (*curr)
    {
        // search for descriptor for example, accessors, materials, images, samplers
        AX_NO_UNROLL while (*curr && *curr != '"') curr++;
        
        if (*curr == '\0') break;

        curr++; // skips the "
        if      (StrCMP16(curr, "accessors"))    curr = ParseAccessors(curr, accessors);
        else if (StrCMP16(curr, "scenes"))       curr = ParseScenes(curr, scenes, stringAllocator, intAllocator);
        else if (StrCMP16(curr, "scene"))        result->defaultSceneIndex = ParsePositiveNumber(curr);
        else if (StrCMP16(curr, "bufferViews"))  curr = ParseBufferViews(curr, bufferViews);
        else if (StrCMP16(curr, "buffers"))      curr = ParseBuffers(curr, path, buffers);     
        else if (StrCMP16(curr, "images"))       curr = ParseImages(curr, path, images, stringAllocator);       
        else if (StrCMP16(curr, "textures"))     curr = ParseTextures(curr, textures, stringAllocator);   
        else if (StrCMP16(curr, "meshes"))       curr = ParseMeshes(curr, meshes, stringAllocator);
        else if (StrCMP16(curr, "materials"))    curr = ParseMaterials(curr, materials, stringAllocator);
        else if (StrCMP16(curr, "nodes"))        curr = ParseNodes(curr, nodes, stringAllocator, intAllocator, scale);
        else if (StrCMP16(curr, "samplers"))     curr = ParseSamplers(curr, samplers);    
        else if (StrCMP16(curr, "cameras"))      curr = ParseCameras(curr, cameras, stringAllocator); 
        else if (StrCMP16(curr, "skins"))        curr = ParseSkins(curr, skins, stringAllocator, intAllocator); 
        else if (StrCMP16(curr, "animations"))   curr = ParseAnimations(curr, animations, stringAllocator, intAllocator); 
        else if (StrCMP16(curr, "asset"))        curr = SkipToNextNode(curr, '{', '}'); // it just has text data that doesn't have anything to do with meshes, (author etc..) if you want you can add this feature :)
        else if (StrCMP16(curr, "extensionsUsed") || StrCMP16(curr, "extensionsRequ")) curr = SkipToNextNode(curr, '[', ']');
        else { ASSERT(0); curr = (const char*)AError_UNKNOWN_DESCRIPTOR; }

        if (curr < (const char*)AError_MAX) // is failed?
        {
            result->error = (AErrorType)(uint64_t)curr;
            FreeAllText(source);
            return 0;
        }
    }

    for (int m = 0; m < meshes.Size(); ++m)
    {
        // get number of vertex, getting first attribute count because all of the others are same
        AMesh mesh = meshes[m];
        mesh.numPrimitives = SBCount(mesh.primitives);
        for (int p = 0; p < mesh.numPrimitives; p++)
        {
            APrimitive& primitive = mesh.primitives[p];
            // get position attrib's count because all attributes same
            int numVertex = accessors[(int)(size_t)primitive.vertexAttribs[0]].count; 
            primitive.numVertices = numVertex;
        
            // get number of index
            GLTFAccessor accessor = accessors[primitive.indiceIndex];
            primitive.numIndices = accessor.count;

            accessor = accessors[primitive.indiceIndex];
            GLTFBufferView view = bufferViews[accessor.bufferView];
            int64_t offset = (int64_t)accessor.byteOffset + view.byteOffset;
            // copy indices
            primitive.indices = ((char*)buffers[view.buffer].uri) + offset;
            primitive.indexType = accessor.componentType;
            
            // get joint data that we need for creating vertices
            const uint32_t jointIndex = TrailingZeroCount32(AAttribType_JOINTS);
            accessor = accessors[(int)(size_t)primitive.vertexAttribs[jointIndex]];
            primitive.jointType   = (short)accessor.componentType;
            primitive.jointCount  = (short)accessor.type;
            primitive.jointStride = (short)bufferViews[accessor.bufferView].byteStride;

            // get weight data that we need for creating vertices
            const uint32_t weightIndex = TrailingZeroCount32(AAttribType_WEIGHTS);
            accessor = accessors[(int)(size_t)primitive.vertexAttribs[weightIndex]];
            primitive.weightType   = (short)accessor.componentType;
            primitive.weightStride = (short)bufferViews[accessor.bufferView].byteStride;

            // position, normal, texcoord are different buffers, 
            // we are unifying all attributes to Vertex* buffer here
            // even though attrib definition in gltf is not ordered, this code will order it, because we traversing set bits
            // for example it converts from this: TexCoord, Normal, Position to Position, Normal, TexCoord
            unsigned attributes = primitive.attributes;
            for (int j = 0; attributes > 0 && j < AAttribType_Count; j += NextSetBit(&attributes))
            {
                accessor     = accessors[(int)(size_t)primitive.vertexAttribs[j]];
                view         = bufferViews[accessor.bufferView];
                offset       = int64_t(accessor.byteOffset) + view.byteOffset;
                
                primitive.vertexAttribs[j] = (char*)buffers[view.buffer].uri + offset;
            }
        }
    }

    for (int s = 0; s < skins.Size(); s++)
    {
        ASkin& skin = skins[s];
        size_t skinIndex = (size_t)skin.inverseBindMatrices;
        GLTFAccessor   accessor  = accessors[(int)skinIndex];
        GLTFBufferView view      = bufferViews[accessor.bufferView];
        int64_t        offset    = int64_t(accessor.byteOffset) + view.byteOffset;
        skin.inverseBindMatrices = (float*)((char*)buffers[view.buffer].uri + offset);
    }

    result->numAnimations = (short)animations.Size();
    for (int a = 0; a < animations.Size(); a++)
    {
        AAnimation& animation = animations[a];
        animation.duration = 0.0f;

        for (int s = 0; s < animation.numSamplers; s++)
        {
            AAnimSampler& sampler = animation.samplers[s];
            size_t inputIndex = (size_t)sampler.input;
            GLTFAccessor   accessor  = accessors[(int)inputIndex];
            GLTFBufferView view      = bufferViews[accessor.bufferView];
            int64_t        offset    = int64_t(accessor.byteOffset) + view.byteOffset;
            
            sampler.input = (float*)((char*)buffers[view.buffer].uri + offset);
            sampler.count = accessor.count;
            
            size_t outputIndex = (size_t)sampler.output;
            accessor = accessors[(int)outputIndex];
            view     = bufferViews[accessor.bufferView];
            offset   = int64_t(accessor.byteOffset) + view.byteOffset;
            
            sampler.output = (float*)((char*)buffers[view.buffer].uri + offset);
            sampler.count = MIN(sampler.count, accessor.count);
            sampler.numComponent = accessor.type;
            
            animation.duration = MAX(animation.duration, sampler.input[sampler.count - 1]);
        }
    }

    // calculate num vertices and indices
    {
        int totalVertexCount = 0;
        int totalIndexCount = 0;
        // Calculate total vertices, and indices
        for (int m = 0; m < meshes.Size(); ++m)
        {
            AMesh mesh = meshes[m];
            for (int p = 0; p < mesh.numPrimitives; p++)
            {
                APrimitive& primitive = mesh.primitives[p];
                totalIndexCount  += primitive.numIndices;
                totalVertexCount += primitive.numVertices;
            }
        }

        result->totalIndices  = totalIndexCount;
        result->totalVertices = totalVertexCount;
    }

    result->stringAllocator = stringAllocator.TakeOwnership();
    result->intAllocator    = intAllocator.TakeOwnership();

    result->numMeshes     = meshes.Size();     result->meshes     = meshes.TakeOwnership();
    result->numNodes      = nodes.Size();      result->nodes      = nodes.TakeOwnership();
    result->numMaterials  = materials.Size();  result->materials  = materials.TakeOwnership();
    result->numTextures   = textures.Size();   result->textures   = textures.TakeOwnership(); 
    result->numImages     = images.Size();     result->images     = images.TakeOwnership();
    result->numSamplers   = samplers.Size();   result->samplers   = samplers.TakeOwnership();
    result->numCameras    = cameras.Size();    result->cameras    = cameras.TakeOwnership();
    result->numScenes     = scenes.Size();     result->scenes     = scenes.TakeOwnership();
    result->numBuffers    = buffers.Size();    result->buffers    = buffers.TakeOwnership();
    result->numAnimations = animations.Size(); result->animations = animations.TakeOwnership();
    result->numSkins      = skins.Size();      result->skins      = skins.TakeOwnership();
    result->scale = scale;
    result->error = AError_NONE;
    FreeAllText(source);
    return 1;
}

__public void FreeGLTFBuffers(SceneBundle* gltf)
{
    for (int i = 0; i < gltf->numBuffers; i++)
    {
        FreeAllText((char*)gltf->buffers[i].uri);
        gltf->buffers[i].uri = nullptr;
    }
    delete[] gltf->buffers;
    gltf->numBuffers = 0;
    gltf->buffers = nullptr;
}

__public void FreeSceneBundle(SceneBundle* gltf)
{
    struct CharFragment { 
        CharFragment* next; char* ptr; int64_t   size;
    };

    struct IntFragment { 
        IntFragment* next; int* ptr; int64_t   size;
    };

    for (int i = 0; i < gltf->numBuffers; i++)
    {
        FreeAllText((char*)gltf->buffers[i].uri);
    }

    if (gltf->stringAllocator)
    {
        // free allocators
        CharFragment* base = (CharFragment*)gltf->stringAllocator;

        while (base)
        {
            delete[] base->ptr;
            CharFragment* oldBase = base;
            base = base->next;
            delete oldBase;
        }
    }
    
    if (gltf->intAllocator)
    {
        IntFragment* ibase = (IntFragment*)gltf->intAllocator;

        while (ibase )
        {
            delete[] ibase->ptr;
            IntFragment* oldBase = ibase ;
            ibase  = ibase->next;
            delete oldBase;
        }
    }

    // also controls if arrays null or not
    for (int i = 0; i < gltf->numMeshes; i++)
        SBFree(gltf->meshes[i].primitives);

    if (gltf->meshes)      delete[] gltf->meshes;
    if (gltf->nodes)       delete[] gltf->nodes;
    if (gltf->materials)   delete[] gltf->materials;
    if (gltf->textures)    delete[] gltf->textures;
    if (gltf->images)      delete[] gltf->images;
    if (gltf->samplers)    delete[] gltf->samplers;
    if (gltf->cameras)     delete[] gltf->cameras;
    if (gltf->scenes)      delete[] gltf->scenes;
    if (gltf->skins)       delete[] gltf->skins;
    if (gltf->animations)
    {
        for (int i = 0; i < gltf->numAnimations; i++)
        {
            delete[] gltf->animations[i].samplers;
            delete[] gltf->animations[i].channels;
        }
        delete[] gltf->animations;
    }
    if (gltf->allVertices) FreeAligned(gltf->allVertices);
    if (gltf->allIndices)  FreeAligned(gltf->allIndices);
    
    MemsetZero(gltf, sizeof(SceneBundle));
}

const char* ParsedSceneGetError(AErrorType error)
{
    const char* SceneParseErrorToStr[] = {"NONE", 
                                          "UNKNOWN",
                                          "UNKNOWN_ATTRIB",
                                          "UNKNOWN_MATERIAL_VAR",
                                          "UNKNOWN_PBR_VAR",
                                          "UNKNOWN_NODE_VAR",
                                          "UNKNOWN_TEXTURE_VAR",
                                          "UNKNOWN_ACCESSOR_VAR",
                                          "UNKNOWN_BUFFER_VIEW_VAR",
                                          "UNKNOWN_MESH_VAR",
                                          "UNKNOWN_CAMERA_VAR",
                                          "UNKNOWN_MESH_PRIMITIVE_VAR",
                                          "BUFFER_PARSE_FAIL",
                                          "BIN_NOT_EXIST",
                                          "FILE_NOT_FOUND",
                                          "UNKNOWN_DESCRIPTOR",
                                          "HASH_COLISSION",
                                          "NON_UTF8",
                                          "MAX" };
    return SceneParseErrorToStr[error];
}