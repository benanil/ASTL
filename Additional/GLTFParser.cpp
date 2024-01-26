
/*****************************************************************
*                                                                *
*    Purpose:                                                    *
*        Simple and efficient parser for GLTF format             *
*        allows you to import 3d mesh, material and scene        *
*    Author:                                                     *
*        Anilcan Gulkaya 2023 anilcangulkaya7@gmail.com          *
*    Restrictions:                                               *
*        No animation and extension support yet.                 *
*    License:                                                    *
*        No License whatsoever do Whatever you want.             *
*                                                                *
*****************************************************************/

#include "GLTFParser.hpp"

#include "../Memory.hpp"
#include "../IO.hpp"
#include "../Algorithms.hpp"
#include "../Array.hpp"
#include "../String.hpp"
#include "../Math/Math.hpp"

#define __private static
#define __public 

struct GLTFAccessor
{
    int bufferView;
    int componentType; // GraphicType
    int count;
    int byteOffset;
    int type; // 1 = SCALAR, 2 = VEC2, 3 = VEC3, 4 = VEC4 
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


#if 0 // AX_SUPPORT_SSE
// _otr must be const char*
// generates bitmask for comparison for example if _otr == "names" 3rd argument is 0x11111
#define StrCMP16(_str, _otr) StrCmp16(_str, _otr, (1 << (sizeof(_otr) - 1)) - 1)
// compares strings length less than 16
// Return 1 if strings are equal, 0 otherwise
inline int StrCmp16(const char* a, const char* b, int cmpMask)
{
    __m128i va = _mm_loadu_si128((const __m128i*)a);
    __m128i vb = _mm_loadu_si128((const __m128i*)b);
    int movemask = _mm_movemask_epi8(_mm_cmpeq_epi8(va, vb));
    return (movemask & cmpMask) == cmpMask;  
}
#else
#define StrCMP16(_str, _otr) StrCmp16(_str, _otr, 0)
// 1 means equal
inline int StrCmp16(const char* a, const char* b, int cmpMask)
{
    bool equal = true;
    while (*b)
        equal &= *a++ == *b++;
    return equal;
}
#endif

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

    while (*curr != '"')
    {
        *alloc++ = *curr++;
    }
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
        ASSERT(*curr != '\0' && "parsing accessors failed probably you forget to close brackets!");
        curr++;
        if (StrCMP16(curr, "bufferView"))    accessor.bufferView = ParsePositiveNumber(curr);
        else if (StrCMP16(curr, "byteOffset"))    accessor.byteOffset = ParsePositiveNumber(curr);
        else if (StrCMP16(curr, "componentType")) accessor.componentType = ParsePositiveNumber(curr) - 0x1400; // GL_BYTE 
        else if (StrCMP16(curr, "count"))         accessor.count = ParsePositiveNumber(curr);
        else if (StrCMP16(curr, "name"))
        {
            curr += 6; // skip name":     because we don't need accessor's name
            int numQuotes = 0;
            // skip two quotes
            while (numQuotes < 2)
                numQuotes += *curr++ == '"';
        }
        else if (StrCMP16(curr, "type"))
        {
            curr += 6; // skip type
            int quoteCnt = 0;
            while (quoteCnt < 2)
            {
                quoteCnt += *curr++ == '"';
            }
            accessor.type = curr[-2] != 'R' ? curr[-2] - '0' : 1; // if scalar 1 otherwise Vec2 or Vec3
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
    curr += 13; // skip bufferViews"

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
        ASSERT(*curr != '0' && "buffer view parse failed, probably you forgot to close brackets!");
        curr++;
        if      (StrCMP16(curr, "buffer"))     bufferView.buffer = ParsePositiveNumber(++curr); 
        else if (StrCMP16(curr, "byteOffset")) bufferView.byteOffset = ParsePositiveNumber(++curr); 
        else if (StrCMP16(curr, "byteLength")) bufferView.byteLength = ParsePositiveNumber(++curr); 
        else if (StrCMP16(curr, "byteStride")) bufferView.byteStride = ParsePositiveNumber(++curr); 
        else if (StrCMP16(curr, "target"))     bufferView.target = ParsePositiveNumber(++curr);         
        else if (StrCMP16(curr, "name")) {
            curr += 6; // we don't need name of the buffer.
            int numQuote = 0;
            while (numQuote < 2)
            numQuote += *curr++ == '"';
        }
        else {
            ASSERT(0 && "UNKNOWN buffer view value!");
            return (const char*)AError_UNKNOWN_BUFFER_VIEW_VAR;
        }
    }
}

__private const char* ParseBuffers(const char* curr, const char* path, Array<GLTFBuffer>& bufferArray)
{
    GLTFBuffer buffer{};
    curr += 9; // skip buffers"
    char binFilePath[256]={0};
    char* endOfWorkDir = binFilePath;
    int binPathlen = StringLength(path);
    SmallMemCpy(binFilePath, path, binPathlen);
    endOfWorkDir = binFilePath + binPathlen;
    
    // remove bla.gltf
    while (binFilePath[binPathlen - 1] != '/' && binFilePath[binPathlen - 1] != '\\')
        binFilePath[--binPathlen] = '\0', endOfWorkDir--;

    // get path of 
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
                MemsetZero(binFilePath, sizeof(binFilePath));
            }

            if (*curr++ == ']') return curr; // end all buffers
        }
        ASSERT(*curr && "parsing buffers failed, probably you forgot to close braces");
        curr++;
        if (StrCMP16(curr, "uri")) // is uri
        {
            curr += 5; // skip uri": 
            while (*curr != '"') curr++;
            curr = GetStringInQuotes(endOfWorkDir, curr);
            buffer.uri = ReadAllText(binFilePath);
            ASSERT(buffer.uri && "uri is not exist");
            if (!buffer.uri) return (const char*)AError_BIN_NOT_EXIST;
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
        
        ASSERT(*curr != '\0' && "parse images failed probably you forgot to close brackets");

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
        ASSERT(*curr != '\0' && "parse images failed probably you forgot to close brackets");
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
            curr = CopyStringInQuotes(texture.name, curr + 7, stringAllocator);
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
        int maskBefore = primitive->attributes;
        if      (StrCMP16(curr, "POSITION"))   primitive->attributes |= AAttribType_POSITION;
        else if (StrCMP16(curr, "NORMAL"))     primitive->attributes |= AAttribType_NORMAL;
        else if (StrCMP16(curr, "TEXCOORD_0")) primitive->attributes |= AAttribType_TEXCOORD_0;
        else if (StrCMP16(curr, "TANGENT"))    primitive->attributes |= AAttribType_TANGENT;
        else if (StrCMP16(curr, "TEXCOORD_1")) primitive->attributes |= AAttribType_TEXCOORD_1;
        else if (StrCMP16(curr, "TEXCOORD_")) { curr = SkipAfter(curr, '"'); continue; }
        // todo add joints and weights
        else { ASSERT(0 && "attribute variable unknown!"); return (const char*)AError_UNKNOWN_ATTRIB; }

        curr = SkipUntill(curr, '"'); curr++;// skip quote because attribute in double quotes
        // using bitmask will help us to order attributes correctly(sort) Position, Normal, TexCoord
        int newIndex = TrailingZeroCount(maskBefore ^ primitive->attributes);
        primitive->vertexAttribs[newIndex] = (void*)(uint64_t)ParsePositiveNumber(curr);
    }
}

__private const char* ParseMeshes(const char* curr, Array<AMesh>& meshes, AStringAllocator& stringAllocator)
{
    char text[64]{};
    curr += 7; // skip meshes" 
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
        else if (!StrCMP16(text, "primitives")) { 
            ASSERT(0 && "only primitives and name allowed"); 
            return (const char*)AError_UNKNOWN_MESH_VAR; 
        }

        APrimitive primitive{};  
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
                }

                if (*curr++ == ']') goto end_primitives; // this is end of primitive list
            }
            curr++;
            
            if      (StrCMP16(curr, "attributes")) { curr = ParseAttributes(curr, &primitive); }
            else if (StrCMP16(curr, "indices"))    { primitive.indiceIndex = ParsePositiveNumber(curr); }
            else if (StrCMP16(curr, "mode"))       { primitive.mode        = ParsePositiveNumber(curr); }
            else if (StrCMP16(curr, "material"))   { primitive.material    = ParsePositiveNumber(curr); }
            else { ASSERT(0); return (const char*)AError_UNKNOWN_MESH_PRIMITIVE_VAR; }
        }
        end_primitives:{}
        curr++; // skip ]
    }
    return nullptr;
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
            }
            if (*curr++ == ']') return curr; // end all nodes
        }
        ASSERT(*curr != '\0' && "parsing nodes not possible, probably forgot to close brackets!");
        curr++; // skips the "
        
        // mesh, name, children, matrix, translation, rotation, scale
        if      (StrCMP16(curr, "mesh"))   { node.type = 0; node.index = ParsePositiveNumber(curr); continue; } // don't want to skip ] that's why continue
        else if (StrCMP16(curr, "camera")) { node.type = 1; node.index = ParsePositiveNumber(curr); continue; } // don't want to skip ] that's why continue
        else if (StrCMP16(curr, "children"))
        {
            // find how many childs there are:
            while (!IsNumber(*curr)) curr++;
           
            const char* begin = curr;
            node.numChildren = 1;
            while (true)
            {
                node.numChildren += *curr == ',';
                if (*curr++ == ']') break;
            }
            curr = begin;
            node.children = intAllocator.AllocateUninitialized(node.numChildren);
            node.numChildren = 0;

            while (*curr != ']')
            {
                if (IsNumber(*curr))
                {
                    node.children[node.numChildren] = ParsePositiveNumber(curr);
                    node.numChildren++;
                }
                curr++;
            }
        }
        else if (StrCMP16(curr, "matrix"))
        {
            float matrix[16]{};
            for (int i = 0; i < 16; i++)
                matrix[i] = ParseFloat(curr);

            node.translation[0] = matrix[12];
            node.translation[1] = matrix[13];
            node.translation[2] = matrix[14];
            QuaternionFromMatrix(node.rotation, matrix);
#ifdef AX_SUPPORT_SSE
            __m128 v;
            v = _mm_load_ps(matrix + 0); node.scale[0] = _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(v, v, 0x7F))) * scale;
            v = _mm_load_ps(matrix + 4); node.scale[1] = _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(v, v, 0x7F))) * scale;
            v = _mm_load_ps(matrix + 8); node.scale[2] = _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(v, v, 0x7F))) * scale;
#else
            node.scale[0] = Sqrt(matrix[0] * matrix[0] + matrix[1] * matrix[1] + matrix[2] * matrix[2]);
            node.scale[1] = Sqrt(matrix[4] * matrix[4] + matrix[5] * matrix[5] + matrix[6] * matrix[6]);
            node.scale[2] = Sqrt(matrix[8] * matrix[8] + matrix[9] * matrix[9] + matrix[10] * matrix[10]);
#endif
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
            curr += 5;
            curr = CopyStringInQuotes(node.name, curr, stringAllocator);
            continue; // continue because we don't want to skip ] and it is not exist
        }
        else
        {
            ASSERT(0 && "Unknown image value uri is the only val!");
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
        ASSERT(*curr != '\0' && "parsing scenes not possible, probably forgot to close brackets!");
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
            curr += 5;
            curr = CopyStringInQuotes(scene.name, curr, stringAllocator);
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
        ASSERT(*curr != '\0' && "parsing nodes not possible, probably forgot to close brackets!");
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
    while (true)
    {
        while (*curr && *curr != '"')
        {
            if (*curr++ == '}') return curr;
        }
        ASSERT(*curr && "parsing material failed, probably forgot to close brackets");
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
            }
            if (*curr++ == ']') return curr; // end all nodes
        }
        ASSERT(*curr && "parsing material failed, probably forgot to close brackets");

        int texture = -1;
        curr++; // skips the "
        if (StrCMP16(curr, "name"))
        {
            curr = CopyStringInQuotes(material.name, curr + 6, stringAllocator);
        }
        else if (StrCMP16(curr, "doubleSided"))
        {
            curr += 12; // skip doubleSided"
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
                    curr = ParseFloat16(curr, material.metallicFactor);
                }
                else if (StrCMP16(curr, "roughnessFact"))
                {
                    curr = ParseFloat16(curr, material.roughnessFactor);
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
            curr = SkipUntill(curr, ','); // skip alpha mode
        }
        else if (StrCMP16(curr, "alphaCutoff"))
        {
            material.alphaCutoff = ParseFloat(curr);
        }
        else {
            ASSERT(0 && "undefined material variable!");
            return (const char*)AError_UNKNOWN_MATERIAL_VAR;
        }

        if (texture != -1)
        {
            curr = ParseMaterialTexture(curr, material.textures[texture]);
            if ((uint64_t)curr == AError_UNKNOWN_MATERIAL_VAR) return curr;
        }
    }
    return curr;
}

__public int ParseGLTF(const char* path, ParsedGLTF* result, float scale)
{
    ASSERT(result && path);
    long sourceSize = 0;
    char* source = ReadAllText(path, nullptr, &sourceSize);
    MemsetZero(result, sizeof(ParsedGLTF));

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

    const char* curr = source;
    while (*curr)
    {
        // search for descriptor for example, accessors, materials, images, samplers
        AX_NO_UNROLL while (*curr && *curr != '"') curr++;
        
        if (*curr == '\0') break;

        curr++; // skips the "
        if      (StrCMP16(curr, "accessors"))    curr = ParseAccessors(curr, accessors);
        else if (StrCMP16(curr, "scenes"))       curr = ParseScenes(curr, scenes, stringAllocator, intAllocator); // todo add scenes
        else if (StrCMP16(curr, "scene"))        result->defaultSceneIndex = ParsePositiveNumber(curr);
        else if (StrCMP16(curr, "bufferViews"))  curr = ParseBufferViews(curr, bufferViews);
        else if (StrCMP16(curr, "buffers"))      curr = ParseBuffers(curr, path, buffers);     
        else if (StrCMP16(curr, "images"))       curr = ParseImages(curr, path, images, stringAllocator);       
        else if (StrCMP16(curr, "textures"))     curr = ParseTextures(curr, textures, stringAllocator);   
        else if (StrCMP16(curr, "meshes"))       curr = ParseMeshes(curr, meshes, stringAllocator);
        else if (StrCMP16(curr, "materials"))    curr = ParseMaterials(curr, materials, stringAllocator);
        else if (StrCMP16(curr, "nodes"))        curr = ParseNodes(curr, nodes, stringAllocator, intAllocator, scale);
        else if (StrCMP16(curr, "samplers"))     curr = ParseSamplers(curr, samplers);    
        else if (StrCMP16(curr, "cameras"))      curr = ParseCameras(curr, cameras, stringAllocator); // todo cameras
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
            // get position attrib's count' because all attributes same
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
            
            // position, normal, texcoord are different buffers, 
            // we are unifying all attributes to Vertex* buffer here
            int attributes = primitive.attributes;
            
            // even though attrib definition in gltf is not ordered, this code will order it, because we traversing set bits
            // for example it converts from this: TexCoord, Normal, Position to Position, Normal, TexCoord
            const int supportedAttributes = 6;
            for (int j = 0; attributes > 0 && j < supportedAttributes; j += NextSetBit(&attributes))
            {
                accessor     = accessors[(int)(size_t)primitive.vertexAttribs[j]];
                view         = bufferViews[accessor.bufferView];
                offset       = int64_t(accessor.byteOffset) + view.byteOffset;
                
                primitive.vertexAttribs[j] = (char*)buffers[view.buffer].uri + offset;
            }
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

    result->numMeshes    = meshes.Size();     result->meshes    = meshes.TakeOwnership();   
    result->numNodes     = nodes.Size();      result->nodes     = nodes.TakeOwnership();    
    result->numMaterials = materials.Size();  result->materials = materials.TakeOwnership();
    result->numTextures  = textures.Size();   result->textures  = textures.TakeOwnership(); 
    result->numImages    = images.Size();     result->images    = images.TakeOwnership();   
    result->numSamplers  = samplers.Size();   result->samplers  = samplers.TakeOwnership();
    result->numCameras   = cameras.Size();    result->cameras   = cameras.TakeOwnership();
    result->numScenes    = scenes.Size();     result->scenes    = scenes.TakeOwnership();
    result->numBuffers   = buffers.Size();    result->buffers   = buffers.TakeOwnership(); 
    result->scale = scale;
    result->error = AError_NONE;
    FreeAllText(source);
    return 1;
}

__public void FreeParsedGLTF(ParsedGLTF* gltf)
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

    if (gltf->meshes)    delete[] gltf->meshes;
    if (gltf->nodes)     delete[] gltf->nodes;
    if (gltf->materials) delete[] gltf->materials;
    if (gltf->textures)  delete[] gltf->textures;
    if (gltf->images)    delete[] gltf->images;
    if (gltf->samplers)  delete[] gltf->samplers;
    if (gltf->cameras)   delete[] gltf->cameras;
    if (gltf->scenes)    delete[] gltf->scenes;
    if (gltf->allVertices) delete[] (float*)gltf->allVertices;
    if (gltf->allIndices)  FreeAligned(gltf->allIndices);
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