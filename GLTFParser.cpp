
//                                           // 
// simple fast and efficient parser for gltf //
// Anilcan Gulkaya 2023                      //
// No License whatsoever do WTF you want.    //
//                                           //

// Todo make this single file that people can include only one header and they are done.

#include "GLTFParser.hpp"
#include "String.hpp"
#include "Array.hpp"
#include "Algorithms.hpp"
#include "IO.hpp"

#ifdef AX_SUPPORT_SSE
#include "Math/SIMDCommon.hpp" // for SSEVectorLength
#endif

#define __private static
#define __public 


#ifndef NO_STRETCHY_BUFFER_SHORT_NAMES
#define sb_free   stb_sb_free
#define sb_push   stb_sb_push
#define sb_count  stb_sb_count
#define sb_add    stb_sb_add
#define sb_last   stb_sb_last
#endif

#define stb_sb_free(a)         ((a) ? free(stb__sbraw(a)),0 : 0)
#define stb_sb_push(a,v)       (stb__sbmaybegrow(a,1), (a)[stb__sbn(a)++] = (v))
#define stb_sb_count(a)        ((a) ? stb__sbn(a) : 0)
#define stb_sb_add(a,n)        (stb__sbmaybegrow(a,n), stb__sbn(a)+=(n), &(a)[stb__sbn(a)-(n)])
#define stb_sb_last(a)         ((a)[stb__sbn(a)-1])

#define stb__sbraw(a) ((int *) (void *) (a) - 2)
#define stb__sbm(a)   stb__sbraw(a)[0]
#define stb__sbn(a)   stb__sbraw(a)[1]

#define stb__sbneedgrow(a,n)  ((a)==0 || stb__sbn(a)+(n) >= stb__sbm(a))
#define stb__sbmaybegrow(a,n) (stb__sbneedgrow(a,(n)) ? stb__sbgrow(a,n) : 0)
#define stb__sbgrow(a,n)      (*((void **)&(a)) = stb__sbgrowf((a), (n), sizeof(*(a))))

#pragma warning (disable : 6011)

static void* stb__sbgrowf(void* arr, int increment, int itemsize)
{
    int dbl_cur = arr ? 2 * stb__sbm(arr) : 0;
    int min_needed = stb_sb_count(arr) + increment;
    int m = dbl_cur > min_needed ? dbl_cur : min_needed;
    int* p = (int*)realloc(arr ? stb__sbraw(arr) : 0, itemsize * m + sizeof(int) * 2);
    if (p) {
        if (!arr)
            p[1] = 0;
        p[0] = m;
        return p + 2;
    }
    else {
        // STRETCHY_BUFFER_OUT_OF_MEMORY
        return (void*)(2 * sizeof(int)); // try to force a NULL pointer exception later
    }
}

struct GLTFAccessor
{
    int bufferView;
    int componentType; // GLTFAccessor::ComponentType
    int count;
    int byteOffset;
    int type; // 1 = SCALAR, 2 = VEC2, 3 = VEC3, 4 = VEC4 
    
    enum ComponentType
    {
        ComponentType_BYTE           = 0x1400, // GL_BYTE            
        ComponentType_UNSIGNED_BYTE  = 0x1401, // GL_UNSIGNED_BYTE   
        ComponentType_SHORT          = 0x1402, // GL_SHORT           
        ComponentType_UNSIGNED_SHORT = 0x1403, // GL_UNSIGNED_SHORT  
        ComponentType_INT            = 0x1404, // GL_INT             
        ComponentType_UNSIGNED_INT   = 0x1405, // GL_UNSIGNED_INT    
        ComponentType_FLOAT          = 0x1406  // GL_FLOAT           
    };
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
typedef Array<GLTFAccessor> AccessorArray;
typedef Array<GLTFBufferView> BufferViewArray;
typedef Array<GLTFBuffer> BufferArray;

#ifdef AX_SUPPORT_SSE
// compares strings length less than 16
inline int StrCmp16(const char* a, const char* b)
{
    // Load the input strings into SSE registers
    __m128i xmm_a = _mm_loadu_si128((const __m128i*)a);
    __m128i xmm_b = _mm_loadu_si128((const __m128i*)b);
    // calculate size of string b
    const int mode = _SIDD_UBYTE_OPS | _SIDD_CMP_EQUAL_EACH | _SIDD_LEAST_SIGNIFICANT;
    const __m128i zeros = _mm_setzero_si128();
    int len = _mm_cmpistri(xmm_b, zeros, mode) - 1;
    
    // Compare the characters using SIMD instructions
    __m128i comparison = _mm_cmpeq_epi8(xmm_a, xmm_b);
    // Move the comparison result to the integer register
    int result     = _mm_movemask_epi8(comparison);
    int lenOneBits = (1 << len) - 1;
    // Check if the strings are equal by comparing the mask
    return (result & lenOneBits) == lenOneBits;  // Return 1 if strings are equal, 0 otherwise
}
#else
// 1 means equal
inline int StrCmp16(const char* a, const char* b)
{
    AX_ASSUME(len < 16);
    bool equal = true;
    while (*b)
        equal &= *a++ == *b++;
    return equal;
}
#endif

__private const char* GetStringInQuotes(char* str, const char* curr)
{
    ++curr; // skip quote " 
    while (*curr != '"') *str++ = *curr++;
    *str = '\0'; // null terminate
    return ++curr;
}

__private const char* SkipToNextNode(const char* curr, char open, char close)
{
    int balance = 1;
    // skip to first open brackets
    while (*curr && *curr != open)
        curr++;
    curr++; // skip open

    while (*curr && balance > 0)
    {
        balance += *curr == open;
        balance -= *curr++ == close;
    }
    curr += *curr != '\0';
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
    char* alloc = stringAllocator.AllocateUninitialized(len + 1);
    str = alloc;
    
    while (*curr != '"')
    {
        *alloc++ = *curr++;
    }
    *alloc = '\0'; // null terminate
    *curr++; // skip quote
    return curr;
}

__private const char* ParseFloat16(const char*& curr, short& flt)
{
    while (!IsNumber(*curr)) curr++;
    float f = ParseFloat(curr);
    flt = (short)(f * 1000.0f);
    return curr;
}

__private const char* ParseAccessors(const char* curr, GLTFAccessor*& accessorArray)
{
    printf("parsing accessors \n");
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
                sb_push(accessorArray, accessor);
                SmallMemSet(&accessor, 0, sizeof(GLTFAccessor));
            }

            if (*curr == ']') // end all accessors
                return ++curr;
            curr++;
        }

        curr++;
        if      (StrCmp16(curr, "bufferView"))    accessor.bufferView = ParsePositiveNumber(++curr); 
        else if (StrCmp16(curr, "byteOffset"))    accessor.byteOffset = ParsePositiveNumber(++curr);
        else if (StrCmp16(curr, "componentType")) accessor.componentType = ParsePositiveNumber(++curr); 
        else if (StrCmp16(curr, "count"))         accessor.count = ParsePositiveNumber(++curr); 
        else if (StrCmp16(curr, "name")) 
        {
            curr += 6; // skip name":
            int numQuotes = 0;
            // skip two quotes
            while (numQuotes < 2)
                numQuotes += *curr++ == '"';
        }
        else if (StrCmp16(curr, "type"))
        {
            curr += 6; // skip type
            int quoteCnt = 0;
            while (quoteCnt < 2)
            {
                quoteCnt += *curr++ == '"';
            }
            accessor.type = curr[-2] != 'R' ? curr[-2] - '0' : 1; // if scalar 1 otherwise Vec2 or Vec3
        }
        else if (!(*curr ^ 'm' | curr[1] ^ 'i' | curr[2] ^ 'n')) curr = SkipToNextNode(curr, '[', ']'); // skip min and max
        else if (!(*curr ^ 'm' | curr[1] ^ 'a' | curr[2] ^ 'x')) curr = SkipToNextNode(curr, '[', ']');
        else ASSERT(0);
    }
}

__private const char* ParseBufferViews(const char* curr, GLTFBufferView*& bufferViews)
{
    printf("parsing buffer views \n");
    GLTFBufferView bufferView;
    SmallMemSet(&bufferView, 0, sizeof(GLTFBufferView));
    curr += 13; // skip bufferViews"
    // read each buffer view
    while (true)
    {
        // search for name
        while (*curr && *curr != '"')
        {
            if (*curr == '}') // next buffer view
            {
                sb_push(bufferViews, bufferView);
                SmallMemSet(&bufferView, 0, sizeof(GLTFBufferView));
            }
            if (*curr++ == ']') return curr; // end all buffer views
        }

        curr++;
        if      (StrCmp16(curr, "buffer"))     bufferView.buffer = ParsePositiveNumber(++curr); 
        else if (StrCmp16(curr, "byteOffset")) bufferView.byteOffset = ParsePositiveNumber(++curr); 
        else if (StrCmp16(curr, "byteLength")) bufferView.byteLength = ParsePositiveNumber(++curr); 
        else if (StrCmp16(curr, "byteStride")) bufferView.byteStride = ParsePositiveNumber(++curr); 
        else if (StrCmp16(curr, "target"))     bufferView.target = ParsePositiveNumber(++curr);         
        else if (StrCmp16(curr, "name"))
        {
            curr += 6;
            int numQuote = 0;
            while (numQuote < 2)
              numQuote += *curr++ == '"';
        }
        else ASSERT(0 && "UNKNOWN buffer view value!"); 
    }
}

__private const char* ParseBuffers(const char* curr, const char* path, GLTFBuffer*& bufferArray)
{
    printf("parsing buffer buffers \n");
    GLTFBuffer buffer{};
    curr += 9; // skip buffers"
    char binFilePath[256]{0};
    char* endOfWorkDir = binFilePath;
    int binPathlen = StringLength(path);
    
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
                sb_push(bufferArray, buffer);
                SmallMemSet(&buffer, 0, sizeof(GLTFBuffer));
                SmallMemSet(binFilePath, 0, sizeof(binFilePath));
            }

            if (*curr++ == ']') return curr; // end all buffers
        }
        curr++;
        if (!(*curr ^ 'u' | curr[1] ^ 'r' | curr[2] ^ 'i')) // is uri
        {
            curr += 5; // skip uri": 
            while (*curr != '"') curr++;
            curr = GetStringInQuotes(endOfWorkDir, curr);
            buffer.uri = ReadAllFile(binFilePath);
            ASSERT(buffer.uri && "uri is not exist");
        }
        else if (!(*curr ^ 'b' | curr[1] ^ 'y' | curr[2] ^ 't' | curr[4] ^ 'L'))
        {
            buffer.byteLength = ParsePositiveNumber(++curr);
        }
        else 
            ASSERT(0 && "Unknown buffer variable! byteLength or uri excepted.");
    }
}

// write paths to path buffer, buffer is seperated by null terminators
__private const char* ParseImages(const char* curr, GLTFImage*& images, GLTFStringAllocator& stringAllocator)
{
    printf("parsing buffer images \n");
    while (*curr != '[') curr++;
    curr++;
    
    GLTFImage image{};
    // read each buffer
    while (true)
    {
        // search for name
        while (*curr && *curr != '"')
        {
            if (*curr++ == ']')
            {
                return curr; // end all images
            }
        }

        curr++; // skips the "
        if (StrCmp16(curr, "uri"))
        {
            curr = CopyStringInQuotes(image.path, curr + 4, stringAllocator);
            sb_push(images, image);
        }
        else
        {
            ASSERT(0 && "Unknown image value uri is the only val!");
        }
    }
    return nullptr;
}

__private const char* ParseTextures(const char* curr, GLTFTexture*& textures, GLTFStringAllocator& stringAllocator)
{
    printf("parsing buffer textures \n");
    curr += sizeof("textures\"");
    GLTFTexture texture{};
    // read each buffer
    while (true)
    {
        // search for name
        while (*curr && *curr != '"')
        {
            if (*curr == '}') // next buffer view
            {
                sb_push(textures, texture);
                SmallMemSet(&texture, 0, sizeof(GLTFTexture));
            }

            if (*curr == ']') // end all buffer views
                return ++curr;
            curr++;
        }
        curr++;
        if (!(*curr ^ 's' | curr[1] ^ 'a' | curr[2] ^ 'm')) // sampler?
        {
            texture.sampler = ParsePositiveNumber(++curr);
        }
        else if (!(*curr ^ 's' | curr[1] ^ 'o' | curr[2] ^ 'u' | curr[3] ^ 'r'))
        {
            texture.source = ParsePositiveNumber(++curr);
        }
        else if (!(*curr ^ 'n' | curr[1] ^ 'a' | curr[2] ^ 'm' | curr[3] ^ 'e'))
        {
            curr = CopyStringInQuotes(texture.name, curr + 7, stringAllocator);
        }
        else
            ASSERT(0 && "Unknown buffer variable! sampler, source or name excepted.");
    }
}

__private const char* ParseMeshes(const char* curr, GLTFMesh*& meshes, GLTFStringAllocator& stringAllocator)
{
    printf("parsing buffer meshes \n");
    char text[64]{};
    curr += 7; // skip meshes" 
    // parse all primitives
    while (true)
    {
        // skip untill we find "primitives
        while (*curr && *curr != '"') curr++;

        SmallMemSet(text, 0, sizeof(text));
        curr = GetStringInQuotes(text, curr);
        char* name = nullptr;

        // there is a possibility that name becomes before primitives
        if (StrCmp16(text, "name"))
        {
            curr = CopyStringInQuotes(name, curr, stringAllocator);

            // skip untill we find "primitives
            while (*curr && *curr != '"') curr++;
            SmallMemSet(text, 0, sizeof(text));
            curr = GetStringInQuotes(text, curr);
        }

        if (!StrCmp16(text, "primitives")) return nullptr;

        // skip untill we find second {
        // it will point "attributes" : {
        int curlyCount = 0;
        while (*curr && curlyCount != 2)
        {
            curlyCount += *curr++ == '{';
        }
        curr++; // skip curly

        // parse attribs
        GLTFMesh tmp{};
        sb_push(meshes, tmp);
        GLTFMesh& primitive = sb_last(meshes);

        while (true)
        {
            while (*curr && (IsWhitespace(*curr) || *curr == '\n')) curr++; // skip whitespace
        
            if (*curr == '}')
                break; // this is end of the attributes

            // skip until we get first attribute
            while (*curr && *curr != '"') curr++;
            SmallMemSet(text, 0, sizeof(text));
            curr = GetStringInQuotes(text, curr);
            curr += 1; // skip :
            
            int maskBefore = primitive.attributes;
            if      (StrCmp16(text, "POSITION"))   primitive.attributes |= GLTFAttribType_POSITION;
            else if (StrCmp16(text, "NORMAL"))     primitive.attributes |= GLTFAttribType_NORMAL;
            else if (StrCmp16(text, "TEXCOORD_0")) primitive.attributes |= GLTFAttribType_TEXCOORD_0;
            else if (StrCmp16(text, "TANGENT"))    primitive.attributes |= GLTFAttribType_TANGENT;
            else if (StrCmp16(text, "TEXCOORD_1")) primitive.attributes |= GLTFAttribType_TEXCOORD_1;
            else return (const char*)GLTFError_UNKNOWN_ATTRIB;
            // using bitmask will help us to order attributes correctly(sort) Position, Normal, TexCoord
            int newIndex = TrailingZeroCount(maskBefore ^ primitive.attributes);
            primitive.attribIndices[newIndex] = ParsePositiveNumber(curr);
        }
        primitive.numAttributes = PopCount(primitive.attributes);

        while (true)
        {
            if (*curr == '"')
            {
                SmallMemSet(text, 0, sizeof(text));
                curr = GetStringInQuotes(text, curr);
                curr += 1; // skip :
                if (StrCmp16(text, "indices"))  primitive.indiceIndex = ParsePositiveNumber(curr);
                else if (StrCmp16(text, "mode")) primitive.mode = ParsePositiveNumber(curr);
                else if (StrCmp16(text, "material")) primitive.material = ParsePositiveNumber(curr);
                else ASSERT(0);
            }
            else if (*curr == ']')
            {
                break;
            }
            curr++;
        }
        // search for name if cannot find that means we are at the end of the primitive
        while (true)
        {
            if (*curr == '"') // probably we found name
            {
                SmallMemSet(text, 0, sizeof(text));
                curr = GetStringInQuotes(text, curr);
                if (StrCmp16(text, "name"))
                {
                    curr = CopyStringInQuotes(name, curr, stringAllocator);
                }
                else ASSERT(0);
            }
            else if (*curr == '}')
            {
                break;
            }
            curr++;
        }
        primitive.name = name;
        
        // detect if we have another primitive if so continue parsing else we have parsed all of the primitives
        while (true)
        {
            if (*curr == ']') return curr; // end of meshes
            if (*curr == '{') break; // next primitive exist
            curr++;
        }
    }

    return nullptr;
}

__private const char* ParseNodes(const char* curr,
                                 GLTFNode*& nodes, 
                                 GLTFStringAllocator& stringAllocator,
                                 FixedSizeGrowableAllocator<int>& children)
{
    // mesh, name, children, 
    // matrix, translation, rotation, scale
    while (*curr != '[') curr++;
    curr++; 
    GLTFNode node{};
    // read each node
    while (true)
    {
        // search for name
        while (*curr && *curr != '"')
        {
            if (*curr == '}')
            {
                sb_push(nodes, node);
                SmallMemSet(&node, 0, sizeof(GLTFNode));
            }
            if (*curr++ == ']') return curr; // end all nodes
        }

        curr++; // skips the "
        if      (StrCmp16(curr, "mesh"))   { node.type = 0; node.index = ParsePositiveNumber(curr); }
        else if (StrCmp16(curr, "camera")) { node.type = 1; node.index = ParsePositiveNumber(curr);  }
        else if (StrCmp16(curr, "children"))
        {
            // find how many childs there are:
            while (!IsNumber(*curr)) curr++;
            const char* begin = curr;
            node.numChildren = 0;
            while (true)
            {
                node.numChildren += *curr == ',';
                if (*curr++ == ']') break;
            }
            curr = begin;
            node.children = children.AllocateUninitialized(node.numChildren);
            node.numChildren = 0;

            while (true)
            {
                if (IsNumber(*curr))
                {
                    node.children[node.numChildren++] = ParsePositiveNumber(curr);
                }
                else if (*curr++ == ']') break;
            }
        }
        else if (StrCmp16(curr, "matrix"))
        {
            float matrix[16]{};
            for (int i = 0; i < 16; i++)
            {
                node.translation[0] = ParseFloat(curr);
            }
            node.translation[0] = matrix[12];
            node.translation[1] = matrix[13];
            node.translation[2] = matrix[14];
            QuaternionFromMatrix(node.rotation, matrix);
#ifdef AX_SUPPORT_SSE
            node.scale[0] = SSEVectorLengthf(_mm_load_ps(matrix + 0));
            node.scale[1] = SSEVectorLengthf(_mm_load_ps(matrix + 4));
            node.scale[2] = SSEVectorLengthf(_mm_load_ps(matrix + 8));
#else
            node.scale[0] = Sqrt(matrix[0] * matrix[0] + matrix[1] * matrix[1] + matrix[2] * matrix[2]);
            node.scale[1] = Sqrt(matrix[4] * matrix[4] + matrix[5] * matrix[5] + matrix[6] * matrix[6]);
            node.scale[2] = Sqrt(matrix[8] * matrix[8] + matrix[9] * matrix[9] + matrix[10] * matrix[10]);
#endif
            while (*curr != ']') curr++;
            curr++;
        }
        else if (StrCmp16(curr, "translation"))
        {
            node.translation[0] = ParseFloat(curr);
            node.translation[1] = ParseFloat(curr);
            node.translation[2] = ParseFloat(curr);
            while (*curr != ']') curr++;
            curr++;
        }
        else if (StrCmp16(curr, "rotation"))
        {
            node.rotation[0] = ParseFloat(curr);
            node.rotation[1] = ParseFloat(curr);
            node.rotation[2] = ParseFloat(curr);
            node.rotation[2] = ParseFloat(curr);
            while (*curr != ']') curr++;
            curr++;
        }
        else if (StrCmp16(curr, "scale"))
        {
            node.scale[0] = ParseFloat(curr);
            node.scale[1] = ParseFloat(curr);
            node.scale[2] = ParseFloat(curr);
            while (*curr != ']') curr++;
            curr++;
        }
        else if (StrCmp16(curr, "name"))
        {
            curr += 5;
            curr = CopyStringInQuotes(node.name, curr, stringAllocator);
        }
        else
        {
            ASSERT(0 && "Unknown image value uri is the only val!");
        }
    }
    return nullptr;
}

__private const char* ParseMaterialTexture(const char* curr, GLTFMaterial::Texture& texture)
{
    printf("parsing materials \n");
    while (*curr != '{') curr++;
    curr++;

    while (true)
    {
        while (*curr != '"')
        {
            if (*curr++ == '}') return curr;
        }
        curr++;

        if (StrCmp16(curr, "scale"))
        {
            curr = ParseFloat16(curr, texture.scale);
        }
        else if (StrCmp16(curr, "index"))
        {
            while (!IsNumber(*curr)) curr++;
            texture.index = (char)ParsePositiveNumber(curr);
        }
        else if (StrCmp16(curr, "texCoord"))
        {
            while (!IsNumber(*curr)) curr++;
            texture.texCoord = (char)ParsePositiveNumber(curr);
        }
        else if (StrCmp16(curr, "strength"))
        {
            curr = ParseFloat16(curr, texture.strength);
        }
        else if (StrCmp16(curr, "extensions"))
        {
            curr = SkipToNextNode(curr, '{', '}'); // currently extensions are not supported 
        }
        else
        {
            ASSERT(0 && "unknown material texture value");
            break;
        }
    }

    return nullptr;
}

__private const char* ParseMaterials(const char* curr, GLTFMaterial*& materials, GLTFStringAllocator& stringAllocator)
{
    printf("parsing materials \n");
    // mesh, name, children, 
    // matrix, translation, rotation, scale
    while (*curr != '[') curr++;
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
                sb_push(materials, material);
                SmallMemSet(&material, 0, sizeof(GLTFMaterial));
            }
            if (*curr++ == ']') return curr; // end all nodes
        }

        int texture = -1;
        curr++; // skips the "
        if (StrCmp16(curr, "name"))
        {
            curr = CopyStringInQuotes(material.name, curr + 6, stringAllocator);
        }
        else if (StrCmp16(curr, "doubleSided"))
        {
            curr += 12; // skip doubleSided"
            while (!IsLower(*curr)) curr++;
            material.doubleSided = *curr == 't';
        }
        else if (StrCmp16(curr, "pbrMetallicRoug")) //pbrMetallicRoughhness
        {
            while (*curr != '{') curr++;
            // parse until finish
            while (true)
            {
                // search for name
                while (*curr && *curr != '"')
                {
                    if (*curr++ == '}') { goto pbr_end; }
                }
                curr++; // skips the "
                
                if      (StrCmp16(curr, "baseColorTex"))  { curr = ParseMaterialTexture(curr, material.metallicRoughness.baseColorTexture); }
                else if (StrCmp16(curr, "metallicRough")) { curr = ParseMaterialTexture(curr, material.metallicRoughness.metallicRoughnessTexture); }
                else if (StrCmp16(curr, "baseColorFact"))
                {
                    curr = ParseFloat16(curr, material.metallicRoughness.baseColorFactor[0]);
                    curr = ParseFloat16(curr, material.metallicRoughness.baseColorFactor[1]);
                    curr = ParseFloat16(curr, material.metallicRoughness.baseColorFactor[2]);
                    curr = ParseFloat16(curr, material.metallicRoughness.baseColorFactor[3]);
                    while (*curr != ']') curr++;
                    curr++;
                }
                else if (StrCmp16(curr, "metallicFact"))
                {
                    curr = ParseFloat16(curr, material.metallicRoughness.metallicFactor);
                }
                else if (StrCmp16(curr, "roughnessFact"))
                {
                    curr = ParseFloat16(curr, material.metallicRoughness.roughnessFactor);
                }
                else ASSERT(0 && "pbrMetallicRoughness value!");
            }
            pbr_end: {}
        }
        else if (StrCmp16(curr, "normalTexture"))    texture = 0;
        else if (StrCmp16(curr, "occlusionTextur"))  texture = 1;
        else if (StrCmp16(curr, "emissiveTexture"))  texture = 2;
        else if (StrCmp16(curr, "emissiveFactor")) 
        {
            curr = ParseFloat16(curr, material.emissiveFactor[0]); 
            curr = ParseFloat16(curr, material.emissiveFactor[1]);
            curr = ParseFloat16(curr, material.emissiveFactor[2]);
            while (*curr != ']') curr++;
            curr++;
        }
        else if (StrCmp16(curr, "extensions"))
        {
            curr = SkipToNextNode(curr, '{', '}'); // currently extensions are not supported 
        }
        else ASSERT(0 && "undefined material variable!");

        if (texture != -1)
        {
            curr = ParseMaterialTexture(curr, material.textures[texture]);
        }
    }
    printf("test");
    return curr;
}

struct Vertex
{
    float normal[3];
    float position[3];
};

inline int GLTFComponentTypeToBytes(int componentType)
{
    // BYTE, UNSIGNED_BYTE, SHORT, UNSIGNED_SHORT, INT, UNSIGNED_INT, FLOAT           
    const int sizes[8] { 1, 1, 2, 2, 4, 4, 4 };
    return sizes[componentType - (int)GLTFAccessor::ComponentType_BYTE];
}

__public ParsedGLTF ParseGLTF(const char* path)
{
    char* source = ReadAllFile(path);
    ScopedText scopedSource(source); // this will free after scope means defer(source)
    ParsedGLTF result;
    SmallMemSet(&result, 0, sizeof(ParsedGLTF));

    if (source == nullptr) return result;

    GLTFTexture*  textures  = nullptr;
    GLTFMesh*     meshes    = nullptr;
    GLTFNode*     nodes     = nullptr;
    GLTFMaterial* materials = nullptr;
    GLTFImage*    images    = nullptr;
    
    GLTFBufferView* bufferViews = nullptr;
    GLTFBuffer*     buffers     = nullptr;
    GLTFAccessor*   accessors   = nullptr;

    GLTFStringAllocator stringAllocator(2048);
    FixedSizeGrowableAllocator<int> childs(512);

    const char* curr = source;

    while (*curr)
    {
        // search for descriptor for example, accessors, materials, images, samplers
        while (*curr && *curr != '"') curr++;
        if (*curr == '\0') break;

        curr++; // skips the "
        if      (StrCmp16(curr, "asset"))     curr = SkipToNextNode(curr, '{', '}');            
        else if (StrCmp16(curr, "accessors")) curr = ParseAccessors(curr, accessors);
        else if (StrCmp16(curr, "scenes"))    curr = SkipToNextNode(curr, '[', ']'); // todo add scenes
        else if (StrCmp16(curr, "scene")) 
        {
          // skip  "scene" : 0,
          while (*curr != ',') curr++; // skip to next "            
          curr++;
        }
        else if (StrCmp16(curr, "bufferViews")) curr = ParseBufferViews(curr, bufferViews);
        else if (StrCmp16(curr, "buffers"))     curr = ParseBuffers(curr, path, buffers);     
        else if (StrCmp16(curr, "images"))      curr = ParseImages(curr, images, stringAllocator);       
        else if (StrCmp16(curr, "textures"))    curr = ParseTextures(curr, textures, stringAllocator);   
        else if (StrCmp16(curr, "meshes"))      curr = ParseMeshes(curr, meshes, stringAllocator);
        else if (StrCmp16(curr, "nodes"))       curr = ParseNodes(curr, nodes, stringAllocator, childs);
        else if (StrCmp16(curr, "samplers"))    curr = SkipToNextNode(curr, '[', ']'); // todo samplers
        else if (StrCmp16(curr, "cameras"))     curr = SkipToNextNode(curr, '[', ']'); // todo cameras
        else if (StrCmp16(curr, "materials")) 
            curr = ParseMaterials(curr, materials, stringAllocator);
        else {
            curr = (const char*)GLTFError_UNKNOWN_DESCRIPTOR;
        }
        // todo: samplers

        // is failed?
        if (curr < (const char*)GLTFError_MAX)
        {
            result.error = (GLTFErrorType)(uint64_t)curr;
            return result;
        }
    }

    GLTFAccessor accessor{};

    for (int i = 0, ilen = sb_count(meshes); i < ilen; ++i)
    {
        // get number of vertex, getting first attribute count because all of the others are same
        accessor = accessors[meshes[i].attribIndices[0]];
        int numVertex = accessor.count;
        meshes[i].numVertices = numVertex;
        
        // get number of index
        accessor = accessors[meshes[i].indiceIndex];
        int numIndex = accessor.count;
        meshes[i].numIndices = numIndex;
        
        accessor = accessors[meshes[i].indiceIndex];
        GLTFBufferView view = bufferViews[accessor.bufferView];
        int64_t offset = (int64_t)accessor.byteOffset + view.byteOffset;
        // copy indices
        meshes[i].indices = ((char*)buffers[view.buffer].uri) + offset;
        meshes[i].indexType = accessor.componentType;

        int attributes = meshes[i].attributes;
        // even though attrib definition in gltf is not ordered, this code will order it, 
        // for example it converts from this: TexCoord, Normal, Position to Position, Normal, TexCoord
        int j = 0;
        while (attributes)
        {
            accessor    = accessors[meshes[i].attribIndices[j]];
            view        = bufferViews[accessor.bufferView];
            offset      = int64_t(accessor.byteOffset) + view.byteOffset;
            meshes[i].vertexAttribs[j]  = (char*)buffers[view.buffer].uri + offset;
            meshes[i].attribTypes[j]    = (char)(accessor.componentType - 0x1400); // 0x1400 GL_BYTE
            meshes[i].attribNumComps[j] = accessor.type;
            // traverse set bits instead of traversing each bit
            attributes &= ~1;
            int tz = TrailingZeroCount(attributes);
            attributes >>= tz;
            j += tz;
        }
    }

    result.stringAllocator = stringAllocator.TakeOwnership();
    result.intAllocator    = childs.TakeOwnership();

    result.numMeshes    = sb_count(meshes);    result.meshes    = meshes;
    result.numNodes     = sb_count(nodes);     result.nodes     = nodes;
    result.numMaterials = sb_count(materials); result.materials = materials;
    result.numTextures  = sb_count(textures);  result.textures  = textures;
    result.numImages    = sb_count(images);    result.images    = images;
    result.error = GLTFError_NONE;
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

    if (!gltf.stringAllocator) return;

    // free allocators
    CharFragment* base = (CharFragment*)gltf.stringAllocator;

    while (base)
    {
        delete[] base->ptr;
        CharFragment* oldBase = base;
        base = base->next;
        delete oldBase;
    }
    
    IntFragment* ibase = (IntFragment*)gltf.intAllocator;

    while (ibase )
    {
        delete[] ibase->ptr;
        IntFragment* oldBase = ibase ;
        ibase  = ibase->next;
        delete oldBase;
    }

    if (gltf.meshes)    sb_free(gltf.meshes);
    if (gltf.nodes)     sb_free(gltf.nodes);
    if (gltf.materials) sb_free(gltf.materials);
    if (gltf.textures)  sb_free(gltf.textures) ;
}
