
/*****************************************************************
*                                                                *
*    Purpose:                                                    *
*        Simple and efficient parser for GLTF format             *
*        allows you to import 3d mesh, material and scene        *
*    Author:                                                     *
*        Anilcan Gulkaya 2023 anilcangulkaya7@gmail.com          *
*    Restrictions:                                               *
*        No animation and extension support yet.                 *
*    Warning:                                                    *
*        Complex scenes will couse crash work in progress.       *
*    License:                                                    *
*        No License whatsoever do WTF you want.                  *
*                                                                *
*****************************************************************/

#include "GLTFParser.hpp"

#include <memory.h>
#include <stdint.h>

#include "../Memory.hpp"
#include "../IO.hpp"
#include "../Algorithms.hpp"

#define __private static
#define __public 

// https://github.com/nothings/stb/blob/master/deprecated/stretchy_buffer.h
#define SBFree(a)         ((a) ? free(stb__sbraw(a)),0 : 0)
#define SBPush(a,v)       (stb__sbmaybegrow(a,1), (a)[stb__sbn(a)++] = (v))
#define SBCount(a)        ((a) ? stb__sbn(a) : 0)
#define SBAdd(a,n)        (stb__sbmaybegrow(a,n), stb__sbn(a)+=(n), &(a)[stb__sbn(a)-(n)])
#define SBLast(a)         ((a)[stb__sbn(a)-1])

#define stb__sbraw(a) ((int *) (void *) (a) - 2)
#define stb__sbm(a)   stb__sbraw(a)[0]
#define stb__sbn(a)   stb__sbraw(a)[1]

#define stb__sbneedgrow(a,n)  ((a)==0 || stb__sbn(a)+(n) >= stb__sbm(a))
#define stb__sbmaybegrow(a,n) (stb__sbneedgrow(a,(n)) ? stb__sbgrow(a,n) : 0)
#define stb__sbgrow(a,n)      (*((void **)&(a)) = stb__sbgrowf((a), (n), sizeof(*(a))))

#pragma warning (disable : 6011)

__forceinline int Max3GLTF(int a, int b, int c) 
{
    int res = a > b ? a : b;
    return res > c ? res : c;
}

static void* stb__sbgrowf(void* arr, int increment, int itemsize)
{
    int dbl_cur = arr ? CalculateArrayGrowth(stb__sbm(arr)) : 0;
    int min_needed = SBCount(arr) + increment;
    int m = Max3GLTF(min_needed, dbl_cur, 32);
    int* p;
    int size = (itemsize * m) + (sizeof(int) * 2);
    if (!arr) p = (int*)malloc(size);
    else      p = (int*)realloc((int*)arr-2, size);
    
    if (!arr) p[1] = 0;
    p[0] = m;
    return p + 2;
}

struct GLTFAccessor
{
    int bufferView;
    int componentType; // GLTFAccessor::ComponentType
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

struct GLTFBuffer
{
    void* uri;
    int byteLength;
};

typedef FixedSizeGrowableAllocator<char> GLTFStringAllocator;

#ifdef AX_SUPPORT_SSE
inline int StringLength(const char* s)
{
    __m128i* mem = reinterpret_cast<__m128i*>(const_cast<char*>(s));
    const __m128i zeros = _mm_setzero_si128();
    
    for (int result = 0; /**/; mem++, result += 16) 
    {
        const __m128i data = _mm_loadu_si128(mem);
        const uint8_t mode = _SIDD_UBYTE_OPS | _SIDD_CMP_EQUAL_EACH | _SIDD_LEAST_SIGNIFICANT;

        if (_mm_cmpistrc(data, zeros, mode)) 
        {
            int idx = _mm_cmpistri(data, zeros, mode);
            return result + idx;
        }
    }
}
#else
inline int StringLength(const char* s)
{
    const char* begin = s;
    while (*s) s++;
    return s - begin;
}
#endif

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

__private const char* CopyStringInQuotes(char*& str, const char* curr, GLTFStringAllocator& stringAllocator)
{
    while (*curr != '"') curr++; // find quote
    curr++; // skip "
    // get length in quotes
    const char* quote = curr;
    while (*quote != '"') quote++;
    int len = quote - curr;
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
    flt = (short)(ParseFloat(curr) * 1000.0f);
    return curr;
}

__private const char* ParseAccessors(const char* curr, GLTFAccessor*& accessorArray)
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
                SBPush(accessorArray, accessor);
                MemsetZero(&accessor, sizeof(GLTFAccessor));
            }

            if (*curr == ']') // end all accessors
            return ++curr;
            curr++;
        }
        ASSERT(*curr != '\0' && "parsing accessors failed probably you forget to close brackets!");
        curr++;
        if      (StrCMP16(curr, "bufferView"))    accessor.bufferView = ParsePositiveNumber(curr); 
        else if (StrCMP16(curr, "byteOffset"))    accessor.byteOffset = ParsePositiveNumber(curr);
        else if (StrCMP16(curr, "componentType")) accessor.componentType = ParsePositiveNumber(curr); 
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
        else 
        return (const char*)GLTFError_UNKNOWN_ACCESSOR_VAR;
    }
}

__private const char* ParseBufferViews(const char* curr, GLTFBufferView*& bufferViews)
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
                SBPush(bufferViews, bufferView);
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
            return (const char*)GLTFError_UNKNOWN_BUFFER_VIEW_VAR;
        }
    }
}

#if _WIN32
#ifndef _MINWINBASE_
extern "C"
{
#define GetCurrentDirectory(bufferLength, buffer) GetCurrentDirectoryA(bufferLength, buffer)
#ifndef _PROCESSENV_
    __declspec(dllimport) uint64_t GetCurrentDirectoryA(uint64_t nBufferLength, char* lpBuffer);
#endif
}
#endif
#else
inline uint64_t GetCurrentDirectory(char* buffer, uint64_t bufferSize)
{
    ASSERT(getcwd(buffer, bufferSize));
    return StringLength(buffer);
}
#endif // windows.h included

__private const char* ParseBuffers(const char* curr, const char* path, GLTFBuffer*& bufferArray)
{
    GLTFBuffer buffer{};
    curr += 9; // skip buffers"
    char binFilePath[256]{0};
    char* endOfWorkDir = binFilePath;
    int binPathlen = StringLength(path);
    
#ifdef _WIN32
    constexpr char ASTL_FILE_SEPERATOR = '\\';
#else
    constexpr char ASTL_FILE_SEPERATOR = '/';
#endif

    if (path[0] >= 'C' && path[0] <= 'G' && path[1] == ':') // path is absolute
    {
        SmallMemCpy(binFilePath, path, binPathlen); // copy path to binfilePath: cwd/meshAdress/bla.gltf
    }
    else // path is relative
    {
        // gets the current work directory into bin file path
        uint64_t cwdLen = GetCurrentDirectory(256, binFilePath);
        endOfWorkDir += cwdLen; // make cwd
        *endOfWorkDir++ = ASTL_FILE_SEPERATOR; // make cwd/
        SmallMemCpy(endOfWorkDir, path, binPathlen); // make cwd/meshAdress/bla.gltf
        binPathlen += cwdLen;
    } 
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
                SBPush(bufferArray, buffer);
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
            buffer.uri = ReadAllFile(binFilePath);
            ASSERT(buffer.uri && "uri is not exist");
            if (!buffer.uri) return (const char*)GLTFError_BIN_NOT_EXIST;
        }
        else if (StrCMP16(curr, "byteLength"))
        {
            buffer.byteLength = ParsePositiveNumber(++curr);
        }
        else
        {
            ASSERT(0 && "Unknown buffer variable! byteLength or uri excepted.");
            return (const char*)GLTFError_BUFFER_PARSE_FAIL;
        }
    }
}

// write paths to path buffer, buffer is seperated by null terminators
__private const char* ParseImages(const char* curr, GLTFImage*& images, GLTFStringAllocator& stringAllocator)
{
    curr = SkipUntill(curr, '[');
    curr++;
    
    GLTFImage image{};
    // read each buffer
    while (true)
    {
        // search for name
        while (*curr && *curr != '"')
        if (*curr++ == ']')
        return curr; // end all images
        
        ASSERT(*curr != '\0' && "parse images failed probably you forgot to close brackets");

        curr++;
        ASSERT(StrCMP16(curr, "uri") && "Unknown image value uri is the only val!");
        curr = CopyStringInQuotes(image.path, curr + 4, stringAllocator);
        SBPush(images, image);
    }
    return nullptr;
}

__private const char* ParseTextures(const char* curr, GLTFTexture*& textures, GLTFStringAllocator& stringAllocator)
{
    curr += sizeof("textures'");
    GLTFTexture texture{};

    // read each buffer
    while (true)
    {
        // search for name
        while (*curr && *curr != '"')
        {
            if (*curr == '}') // next buffer view
            {
                SBPush(textures, texture);
                MemsetZero(&texture, sizeof(GLTFTexture));
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
            return (const char*)GLTFError_UNKNOWN_TEXTURE_VAR;
        }
    }
}

__private const char* ParseAttributes(const char* curr, GLTFPrimitive* primitive)
{
    curr += sizeof("attributes'");

    while (true)
    {
        while (*curr != '"')
        if (*curr++ == '}') return curr;
        
        curr++; // skip "
        int maskBefore = primitive->attributes;
        if      (StrCMP16(curr, "POSITION"))   primitive->attributes |= GLTFAttribType_POSITION;
        else if (StrCMP16(curr, "NORMAL"))     primitive->attributes |= GLTFAttribType_NORMAL;
        else if (StrCMP16(curr, "TEXCOORD_0")) primitive->attributes |= GLTFAttribType_TEXCOORD_0;
        else if (StrCMP16(curr, "TANGENT"))    primitive->attributes |= GLTFAttribType_TANGENT;
        else if (StrCMP16(curr, "TEXCOORD_1")) primitive->attributes |= GLTFAttribType_TEXCOORD_1;
        else if (StrCMP16(curr, "TEXCOORD_")) { curr = SkipAfter(curr, '"'); continue; }
        // todo add joints and weights
        else { ASSERT(0 && "attribute variable unknown!"); return (const char*)GLTFError_UNKNOWN_ATTRIB; }

        curr = SkipUntill(curr, '"'); curr++;// skip quote because attribute in double quotes
        // using bitmask will help us to order attributes correctly(sort) Position, Normal, TexCoord
        int newIndex = TrailingZeroCount(maskBefore ^ primitive->attributes);
        primitive->vertexAttribs[newIndex] = (void*)(uint64_t)ParsePositiveNumber(curr);    
    }
}

__private const char* ParseMeshes(const char* curr, GLTFMesh*& meshes, GLTFStringAllocator& stringAllocator)
{
    char text[64]{};
    curr += 7; // skip meshes" 
    GLTFMesh mesh{};
    int numMeshes = 0;
    // parse all meshes
    while (true)
    {
        while (*curr != '"')
        {
            if (*curr == '}') 
            {
                SBPush(meshes, mesh);
                MemsetZero(&mesh, sizeof(GLTFMesh));
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
            return (const char*)GLTFError_UNKNOWN_MESH_VAR; 
        }

        GLTFPrimitive primitive{};  
        // parse primitives
        while (true)
        {
            while (*curr != '"')
            {
                if (*curr == '}')
                {
                    SBPush(mesh.primitives, primitive);
                    MemsetZero(&primitive, sizeof(GLTFPrimitive));
                    mesh.numPrimitives++;
                }

                if (*curr++ == ']') goto end_primitives; // this is end of primitive list
            }
            curr++;
            
            if      (StrCMP16(curr, "attributes")) { curr = ParseAttributes(curr, &primitive); }
            else if (StrCMP16(curr, "indices"))    { primitive.indiceIndex = ParsePositiveNumber(curr); }
            else if (StrCMP16(curr, "mode"))       { primitive.mode        = ParsePositiveNumber(curr); }
            else if (StrCMP16(curr, "material"))   { primitive.material    = ParsePositiveNumber(curr); }
            else { ASSERT(0); return (const char*)GLTFError_UNKNOWN_MESH_PRIMITIVE_VAR; }
        }
        end_primitives:{}
        curr++; // skip ]
    }
    return nullptr;
}

__forceinline float Sqrt(float a)
{
#ifdef AX_SUPPORT_SSE
    return _mm_cvtss_f32(_mm_sqrt_ps(_mm_set_ps1(a)));  
#elif defined(__clang__)
    return __builtin_sqrt(a);
#endif
}

inline void QuaternionFromMatrix(float* Orientation, const float* m)
{
    int i, j, k = 0;
    const int numCol = 4;

    float root, trace = m[0*numCol+0] + m[1 * numCol + 1] + m[2 * numCol + 2];

    if (trace > 0.0f)
    {
        root = Sqrt(trace + 1.0f);
        Orientation[3] = 0.5f * root;
        root = 0.5f / root;
        Orientation[0] = root * (m[1 * numCol + 2] - m[2 * numCol + 1]);
        Orientation[1] = root * (m[2 * numCol + 0] - m[0 * numCol + 2]);
        Orientation[2] = root * (m[0 * numCol + 1] - m[1 * numCol + 0]);
    }
    else
    {
        static const int Next[3] = { 1, 2, 0 };
        i = 0;
        i += m[1 * numCol + 1] > m[0 * numCol + 0]; // if (M.m[1][1] > M.m[0][0]) i = 1
        if (m[2 * numCol + 2] > m[i * numCol + i]) i = 2;
        j = Next[i];
        k = Next[j];

        root = Sqrt(m[i * numCol + i] - m[j * numCol + j] - m[k * numCol + k] + 1.0f);

        Orientation[i] = 0.5f * root;
        root = 0.5f / root;
        Orientation[j] = root * (m[i * numCol + j] + m[j * numCol + i]);
        Orientation[k] = root * (m[i * numCol + k] + m[k * numCol + i]);
        Orientation[3] = root * (m[j * numCol + k] - m[k*numCol+j]);
    } 
}

__private const char* ParseNodes(const char* curr,
                                 GLTFNode*& nodes,
                                 GLTFStringAllocator& stringAllocator,
                                 FixedSizeGrowableAllocator<int>& intAllocator)
{
    curr = SkipUntill(curr, '[');
    curr++;
    GLTFNode node{};
    node.rotation[3] = 1.0f;
    node.scale[0] = node.scale[1] = node.scale[2] = 1.0f; 
    // read each node
    while (true)
    {
        // search for name
        while (*curr && *curr != '"')
        {
            if (*curr == '}')
            {
                SBPush(nodes, node);
                MemsetZero(&node, sizeof(GLTFNode));
                node.rotation[3] = 1.0f;
                node.scale[0] = node.scale[1] = node.scale[2] = 1.0f;
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
            node.translation[0] = ParseFloat(curr);

            node.translation[0] = matrix[12];
            node.translation[1] = matrix[13];
            node.translation[2] = matrix[14];
            QuaternionFromMatrix(node.rotation, matrix);
#ifdef AX_SUPPORT_SSE
            __m128 v;
            v = _mm_load_ps(matrix + 0); node.scale[0] = _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(v, v, 0x71)));
            v = _mm_load_ps(matrix + 4); node.scale[1] = _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(v, v, 0x71)));
            v = _mm_load_ps(matrix + 8); node.scale[2] = _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(v, v, 0x71)));
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
            node.scale[0] = ParseFloat(curr);
            node.scale[1] = ParseFloat(curr);
            node.scale[2] = ParseFloat(curr);
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
            return (const char*)GLTFError_UNKNOWN_NODE_VAR;
        }

        curr = SkipUntill(curr, ']');
        curr++;
    }
    return nullptr;
}

__private const char* ParseCameras(const char* curr, GLTFCamera*& cameras, GLTFStringAllocator& stringAllocator)
{
    curr += sizeof("camera'");
    char text[64]{};
    GLTFCamera camera{};
    // parse all meshes
    while (true)
    {
        while (*curr != '"')
        {
            if (*curr == '}') 
            {
                SBPush(cameras, camera);
                MemsetZero(&camera, sizeof(GLTFCamera));
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
            return (const char*)GLTFError_UNKNOWN_CAMERA_VAR; 
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
            else { ASSERT(0); return (const char*)GLTFError_UNKNOWN_CAMERA_VAR; }
        }
        end_properties:{}
    }
    return nullptr;
}

__private const char* ParseScenes(const char* curr, GLTFScene*& scenes, 
                                  GLTFStringAllocator& stringAllocator, FixedSizeGrowableAllocator<int>& intAllocator)
{
    curr = SkipUntill(curr, '[');
    curr++;
    GLTFScene scene{};
    // read each node
    while (true)
    {
        // search for name
        while (*curr != '"')
        {
            if (*curr == '}')
            {
                SBPush(scenes, scene);
                MemsetZero(&scene, sizeof(GLTFScene));
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

__private const char* ParseSamplers(const char* curr, GLTFSampler*& samplers)
{
    curr = SkipUntill(curr, '[');
    curr++;
    
    GLTFSampler sampler{};
    
    // read each node
    while (true)
    {
        // search for name
        while (*curr && *curr != '"')
        {
            if (*curr == '}')
            {
                SBPush(samplers, sampler);
                MemsetZero(&sampler, sizeof(GLTFSampler));
            }
            if (*curr++ == ']') return curr; // end all nodes
        }
        ASSERT(*curr != '\0' && "parsing nodes not possible, probably forgot to close brackets!");
        curr++; // skips the "

        if      (StrCMP16(curr, "magFilter")) sampler.magFilter = (char)(ParsePositiveNumber(curr) - 0x2600); // GL_NEAREST 9728, GL_LINEAR 0x2601 9729
        else if (StrCMP16(curr, "minFilter")) sampler.minFilter = (char)(ParsePositiveNumber(curr) - 0x2600); // GL_NEAREST 9728, GL_LINEAR 0x2601 9729
        else if (StrCMP16(curr, "wrapS"))     sampler.wrapS = (char)OGLWrapToWrap(ParsePositiveNumber(curr));
        else if (StrCMP16(curr, "wrapT"))     sampler.wrapT = (char)OGLWrapToWrap(ParsePositiveNumber(curr));
        else { ASSERT(0 && "parse samplers failed!"); return (const char*)GLTFError_UNKNOWN; }
    }
}

__private const char* ParseMaterialTexture(const char* curr, GLTFMaterial::Texture& texture)
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
            return (const char*)GLTFError_UNKNOWN_MATERIAL_VAR;
        }
    }
    return nullptr;
}

__private const char* ParseMaterials(const char* curr, GLTFMaterial*& materials, GLTFStringAllocator& stringAllocator)
{
    // mesh, name, children, 
    // matrix, translation, rotation, scale
    curr = SkipUntill(curr, '[');
    curr++; // skip '['
    GLTFMaterial material{};
    // read each material
    while (true)
    {
        // search for name
        while (*curr && *curr != '"')
        {
            if (*curr == '}')
            {
                SBPush(materials, material);
                MemsetZero(&material, sizeof(GLTFMaterial));
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
                
                if      (StrCMP16(curr, "baseColorTex"))  { curr = ParseMaterialTexture(curr, material.metallicRoughness.baseColorTexture); }
                else if (StrCMP16(curr, "metallicRough")) { curr = ParseMaterialTexture(curr, material.metallicRoughness.metallicRoughnessTexture); }
                else if (StrCMP16(curr, "baseColorFact"))
                {
                    curr = ParseFloat16(curr, material.metallicRoughness.baseColorFactor[0]);
                    curr = ParseFloat16(curr, material.metallicRoughness.baseColorFactor[1]);
                    curr = ParseFloat16(curr, material.metallicRoughness.baseColorFactor[2]);
                    curr = ParseFloat16(curr, material.metallicRoughness.baseColorFactor[3]);
                    curr = SkipUntill(curr, ']');
                    curr++;
                }
                else if (StrCMP16(curr, "metallicFact"))
                {
                    curr = ParseFloat16(curr, material.metallicRoughness.metallicFactor);
                }
                else if (StrCMP16(curr, "roughnessFact"))
                {
                    curr = ParseFloat16(curr, material.metallicRoughness.roughnessFactor);
                }
                else
                {
                    ASSERT(0 && "unknown pbrMetallicRoughness value!");
                    return (char*)GLTFError_UNKNOWN_PBR_VAR;
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
        else {
            ASSERT(0 && "undefined material variable!");
            return (const char*)GLTFError_UNKNOWN_MATERIAL_VAR;
        }

        if (texture != -1)
        {
            curr = ParseMaterialTexture(curr, material.textures[texture]);
            if ((uint64_t)curr == GLTFError_UNKNOWN_MATERIAL_VAR) return curr;
        }
    }
    return curr;
}

__public ParsedGLTF ParseGLTF(const char* path)
{
    ParsedGLTF result{};
    char* source = ReadAllFile(path);
    if (source == nullptr) { result.error = GLTFError_FILE_NOT_FOUND; return result; }

    GLTFBufferView* bufferViews = nullptr;
    GLTFBuffer*     buffers     = nullptr;
    GLTFAccessor* accessors{};
    GLTFAccessor accessor{};

    GLTFStringAllocator stringAllocator(2048);
    FixedSizeGrowableAllocator<int> intAllocator(512);

    const char* curr = source;
    while (*curr)
    {
        // search for descriptor for example, accessors, materials, images, samplers
        AX_NO_UNROLL while (*curr && *curr != '"') curr++;
        
        if (*curr == '\0') break;

        curr++; // skips the "
        if      (StrCMP16(curr, "accessors"))    curr = ParseAccessors(curr, accessors);
        else if (StrCMP16(curr, "scenes"))       curr = ParseScenes(curr, result.scenes, stringAllocator, intAllocator); // todo add scenes
        else if (StrCMP16(curr, "scene"))        result.defaultSceneIndex = ParsePositiveNumber(curr);
        else if (StrCMP16(curr, "bufferViews"))  curr = ParseBufferViews(curr, bufferViews);
        else if (StrCMP16(curr, "buffers"))      curr = ParseBuffers(curr, path, buffers);     
        else if (StrCMP16(curr, "images"))       curr = ParseImages(curr, result.images, stringAllocator);       
        else if (StrCMP16(curr, "textures"))     curr = ParseTextures(curr, result.textures, stringAllocator);   
        else if (StrCMP16(curr, "meshes"))       curr = ParseMeshes(curr, result.meshes, stringAllocator);
        else if (StrCMP16(curr, "materials"))    curr = ParseMaterials(curr, result.materials, stringAllocator);
        else if (StrCMP16(curr, "nodes"))        curr = ParseNodes(curr, result.nodes, stringAllocator, intAllocator);
        else if (StrCMP16(curr, "samplers"))     curr = ParseSamplers(curr, result.samplers);    
        else if (StrCMP16(curr, "cameras"))      curr = ParseCameras(curr, result.cameras, stringAllocator); // todo cameras
        else if (StrCMP16(curr, "asset"))        curr = SkipToNextNode(curr, '{', '}'); // it just has text data that doesn't have anything to do with meshes, (author etc..) if you want you can add this feature :)
        else if (StrCMP16(curr, "extensionsUsed") || StrCMP16(curr, "extensionsRequ")) curr = SkipToNextNode(curr, '[', ']');
        else { ASSERT(0); curr = (const char*)GLTFError_UNKNOWN_DESCRIPTOR; }

        if (curr < (const char*)GLTFError_MAX) // is failed?
        {
            result.error = (GLTFErrorType)(uint64_t)curr;
            goto finish_parse;
        }
    }

    for (int m = 0, mlen = SBCount(result.meshes); m < mlen; ++m)
    {
        // get number of vertex, getting first attribute count because all of the others are same
        GLTFMesh mesh = result.meshes[m];
        mesh.numPrimitives = SBCount(mesh.primitives);
        for (uint64_t p = 0; p < mesh.numPrimitives; p++)
        {
            GLTFPrimitive& primitive = mesh.primitives[p];
            // get position attrib's count' because all attributes same
            int numVertex = accessors[(uint64_t)primitive.vertexAttribs[0]].count; 
            primitive.numVertices = numVertex;
        
            // get number of index
            accessor = accessors[primitive.indiceIndex];
            int numIndex = accessor.count;
            primitive.numIndices = numIndex;
        
            accessor = accessors[primitive.indiceIndex];
            GLTFBufferView view = bufferViews[accessor.bufferView];
            int64_t offset = (int64_t)accessor.byteOffset + view.byteOffset;
            // copy indices
            primitive.indices = ((char*)buffers[view.buffer].uri) + offset;
            primitive.indexType = accessor.componentType;

            int attributes = primitive.attributes;
            // even though attrib definition in gltf is not ordered, this code will order it, 
            // for example it converts from this: TexCoord, Normal, Position to Position, Normal, TexCoord
            int j = 0;
            while (attributes)
            {
                accessor    = accessors[(uint64_t)primitive.vertexAttribs[j]];
                view        = bufferViews[accessor.bufferView];
                offset      = int64_t(accessor.byteOffset) + view.byteOffset;
                primitive.vertexAttribs[j]  = (char*)buffers[view.buffer].uri + offset;
                // traverse set bits instead of traversing each bit
                attributes &= ~1;
                int tz = TrailingZeroCount(attributes);
                attributes >>= tz;
                j += tz;
            }    
        }
    }

    result.stringAllocator = stringAllocator.TakeOwnership();
    result.intAllocator    = intAllocator.TakeOwnership();
    result.numMeshes    = SBCount(result.meshes);   
    result.numNodes     = SBCount(result.nodes);    
    result.numMaterials = SBCount(result.materials);
    result.numTextures  = SBCount(result.textures); 
    result.numImages    = SBCount(result.images);   
    result.numSamplers  = SBCount(result.samplers);
    result.numCameras   = SBCount(result.cameras);
    result.numScenes    = SBCount(result.scenes);

    result.error = GLTFError_NONE;
    finish_parse:
    {
        SBFree(bufferViews);
        SBFree(buffers);
        SBFree(accessors);
        free(source);
    }
    return result;
}

__public void FreeGLTF(ParsedGLTF& gltf)
{
    struct CharFragment { 
        CharFragment* next; char* ptr; int64_t   size;
    };

    struct IntFragment { 
        IntFragment* next; int* ptr; int64_t   size;
    };

    if (gltf.stringAllocator)
    {
        // free allocators
        CharFragment* base = (CharFragment*)gltf.stringAllocator;

        while (base)
        {
            delete[] base->ptr;
            CharFragment* oldBase = base;
            base = base->next;
            delete oldBase;
        }
    }
    
    if (gltf.intAllocator)
    {
        IntFragment* ibase = (IntFragment*)gltf.intAllocator;

        while (ibase )
        {
            delete[] ibase->ptr;
            IntFragment* oldBase = ibase ;
            ibase  = ibase->next;
            delete oldBase;
        }
    }

    // also controls if arrays null or not
    SBFree(gltf.meshes);
    SBFree(gltf.nodes);
    SBFree(gltf.materials);
    SBFree(gltf.textures);
    SBFree(gltf.images);
    SBFree(gltf.samplers);
    SBFree(gltf.cameras);
    SBFree(gltf.scenes);
}