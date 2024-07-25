
/*****************************************************************
*                                                                *
*    Purpose:                                                    *
*        Simple and efficient parser for OBJ format              *
*        allows you to import 3d mesh and material               *
*    Author:                                                     *
*        Anilcan Gulkaya 2023 anilcangulkaya7@gmail.com          *
*    License:                                                    *
*        No License whatsoever do Whatever you want.             *
*                                                                *
*****************************************************************/

#include "GLTFParser.hpp"
#include "../Algorithms.hpp"
#include "../Common.hpp"
#include "../IO.hpp"
#include "../String.hpp"
#include "../Random.hpp" // WYHash, MurmurHash
#include "../HashMap.hpp"
#include "../HashSet.hpp"

inline void ChangeExtension(char* path, const char* newExt, size_t len)
{
    path[len - 1] = newExt[2]; path[len - 2] = newExt[1]; path[len - 3] = newExt[0];
}

inline float ParseFloatNonConst(char*& text)
{
    char* ptr = text;
    while (!IsNumber(*ptr) && *ptr != '-') ptr++;

    double sign = 1.0;
    if (*ptr == '-') sign = -1.0, ptr++;

    double num = 0.0;

    while (IsNumber(*ptr))
        num = 10.0 * num + (double)(*ptr++ - '0');

    if (*ptr == '.') ptr++;

    double fra = 0.0, div = 1.0;

    while (IsNumber(*ptr) && div < 1e8) // 1e8 is 1 and 8 zero 100000000
        fra = 10.0f * fra + (double)(*ptr++ - '0'), div *= 10.0f;

    while (IsNumber(*ptr)) ptr++;

    num += fra / div;
    text = ptr;
    return (float)(sign * num);
}

static void ParseFloat162(char*& curr, short& flt)
{
    flt = (short)(ParseFloatNonConst(curr) * 400.0f);
}

struct OBJVertex { float pos[3]; float texCoord[2]; float normal[3]; };


inline uint64_t HashVertex(int pos, int tex, int norm)
{
    const uint64_t secret[4] = { 0xa0761d6478bd642full, 0xe7037ed1a0b428dbull, 0x8ebc6af09c88c6e3ull, 0x589965cc75374cc3ull };
    uint64_t a = (uint64_t)pos | (uint64_t(tex) << 22) | (uint64_t(norm) << 39);
    return WYHash::mix(a ^ secret[a & 3], MurmurHash(a));
    // uint64_t a = (uint64_t(pos) << 32U) | uint64_t(tex);
    // uint64_t b = (uint64_t(norm) << 32U) | uint64_t(pos + tex * norm);
    // return WYHash::mix(secret[1], WYHash::mix(a ^ secret[1], b ^ tex));
}   

int ParseObj(const char* path, ParsedObj* scene)
{
    if (!FileExist(path)) {
        scene->error = AError_FILE_NOT_FOUND;
        return 0;
    }

    int pathLen = StringLength(path);
    FixedSizeGrowableAllocator<char> stringAllocator(512);
    
    char* pathDup = stringAllocator.AllocateUninitialized(pathLen + 8);
    MemsetZero(pathDup, pathLen + 8);
    SmallMemCpy(pathDup, path, pathLen);

    const uint64 sz = FileSize(path);
    ScopedText objText = ReadAllText(path); 
    
    ChangeExtension(pathDup, "mtl", pathLen);
    char* mtlPath = pathDup; // for readibility
    uint msz = FileExist(mtlPath) ? (uint)FileSize(mtlPath) : 0u;
    
    char* mtlText = nullptr;
    // read if material path exist 
    if (msz) {
        mtlText = ReadAllText(mtlPath); // todo memory leak
        scene->materialText = mtlText;
    }

    Array<float> positions, texCoords, normals; 
    Array<AMesh> meshes;

    // hashMap for material indices, materialMap[materialNameHash & 255] = material Index 
    StackHashMap<uint32_t, int, 128, NoHasher<uint32_t>> materialMap;
	
    // import materials
    char* curr = mtlText, *currEnd = curr + msz;
    Array<AMaterial> materials;
    Array<AImage> images; // index to start of mtltext
    
    while (curr && *curr && curr < currEnd)
    {
        if (*curr == '#')  while (*curr++ != '\n'); // skip line
        while (*curr == '\n' || IsWhitespace(*curr)) curr++;
		
        if (curr[0] == 'n' && curr[1] == 'e' && curr[2] == 'w') // no need to m, t, l 
        {
            curr += 7; // skip newmtl + space
            materials.AddUninitialized(1);
            MemsetZero(&materials.Back(), sizeof(AMaterial));

            // set default properties
            materials.Back().specularColor   = ~0u, materials.Back().diffuseColor = ~0u; // white
            materials.Back().specularFactor  = AMaterial::MakeFloat16(2.2f);
            materials.Back().roughnessFactor = AMaterial::MakeFloat16(0.6f);
            materials.Back().name = mtlText + PointerDistance(mtlText, curr);
            
            uint32_t hash = WangHash(uint(curr[0]) | uint(curr[1] << 8) | uint(curr[2] << 16));
            while(*curr != '\n' && !IsWhitespace(*curr))
                hash = *curr++ + (hash << 6) + (hash << 16) - hash; 
            *curr++ = '\0'; // null terminator
        
            if (materialMap.Contains(hash))
            { scene->error = AError_HASH_COLISSION; return 0;}

            materialMap[hash] = scene->numMaterials++;
        }
        else if (curr[0] == 'N' && curr[1] == 's') 
        {
            // shininess
            materials.Back().specularFactor = PackUnorm16(ParseFloatNonConst(curr));
        }
        else if (curr[0] == 'd') 
        {
            // roughness
            materials.Back().roughnessFactor = PackUnorm16(ParseFloatNonConst(curr));
        }
        else if (curr[0] == 'K' && curr[1] == 'd') 
        { 
            // diffuse color
            float colorf[3] = { ParseFloatNonConst(curr), ParseFloatNonConst(curr), ParseFloatNonConst(curr)};
            materials.Back().baseColorFactor = PackColor3PtrToUint(colorf);
        }
        else if (curr[0] == 'K' && curr[1] == 's') { // specular color
            float colorf[3] = { ParseFloatNonConst(curr), ParseFloatNonConst(curr), ParseFloatNonConst(curr)};
            materials.Back().specularColor = PackColor3PtrToUint(colorf);
        }
        else if (curr[0] == 'm') // map_bla
        {
            if (curr[4] == 'K' && curr[5] == 'd') {
                // set this pointer as texture path data
                materials.Back().baseColorTexture.index = images.Size();
                AImage image{ mtlText + PointerDistance(mtlText, (curr += 7))};
                images.Add(image);
                // find end point of path
                while (*curr != '\n' && *curr != '\r') curr++;
                *curr++ = '\0'; // addd null terminator
            }
            else if (curr[4] == 'K' && curr[5] == 's') {
                materials.Back().specularTexture.index = images.Size();
                AImage image{ mtlText + PointerDistance(mtlText, (curr += 7))};
                images.Add(image);

                while (*curr != '\n' && *curr != '\r') curr++; // find end point of path
                *curr++ = '\0'; // addd null terminator
            }
            else while (*curr != '\n') curr++;
        } // skip line. header is unknown | unused
        else while (*curr != '\n' && *curr) curr++;
    }

    curr = objText.text, currEnd = curr + sz;
    unsigned currentMaterial = 0;
    
    HashMap<uint64_t, uint32_t, NoHasher<uint64_t>> uniqueVertices;
    Array<uint32_t> indices(256);
    Array<OBJVertex> vertices;
    
    int attributes = 0; // pos, tex, normal bits.
    
    while (curr && *curr && curr < currEnd)
    {
        if (*curr == 'x')
            printf("break");

        while (*curr && (*curr == '\n' || IsWhitespace(*curr))) 
            curr++;

        // skip comment lines
        if (*curr == '#') 
            while (*curr && *curr != '\n')
                curr++;

        if (*curr == 'v')
        {
            attributes |= AAttribType_POSITION;

            while (*curr != '#' && curr[1] == ' ') { // vertex=position 
                curr += 2; 
                positions.Add(ParseFloatNonConst(curr));
                positions.Add(ParseFloatNonConst(curr));
                positions.Add(ParseFloatNonConst(curr));
                // skip line&whiteSpace
                while (*curr == '\n' || IsWhitespace(*curr)) curr++;
            }

            while (curr[1] == 't') {
                attributes |= AAttribType_TEXCOORD_0;

                curr += 2;
                texCoords.Add(ParseFloatNonConst(curr));
                texCoords.Add(ParseFloatNonConst(curr));
                // if it has 3 component for texcoord skip the third
                if (IsNumber(curr[1])) ParseFloatNonConst(curr);
                // skip line&whiteSpace
                while (*curr == '\n' || IsWhitespace(*curr)) curr++;
            }

            while (curr[1] == 'n') {
                attributes |= AAttribType_NORMAL;
                curr += 2; 
                normals.Add(ParseFloatNonConst(curr));
                normals.Add(ParseFloatNonConst(curr));
                normals.Add(ParseFloatNonConst(curr));
                // skip line&whiteSpace
                while (*curr == '\n' || IsWhitespace(*curr)) curr++;
            }
        } 
		
        int materialIndex = 0;
        if (curr[0] == 'u' && curr[1] == 's' && curr[2] == 'e') // && curr[3] == 'm' && curr[4] == 't' no need
        {
            curr += 7; // skip usemtl + space
            unsigned hash = WangHash(uint(curr[0]) | uint(curr[1] << 8) | uint(curr[2] << 16));
            while(*curr != '\n' && !IsWhitespace(*curr)) // create texture path hash 
                hash = *curr++ + (hash << 6) + (hash << 16) - hash;
            while (*curr == '\n' || IsWhitespace(*curr)) curr++;
            materialIndex = materialMap[hash]; // use material index for this group of triangles
        }
        
        int vertexStart  = vertices.Size();
        int indiceStart  = indices.Size();
        int currentIndex = indiceStart;

        while (curr[0] == 'f')
        {
            curr += 2;
      
            // iterate over vaces
            for(int i = 0; i < 3; ++i) // compiler please unroll :D
            {
                int positionIdx = 0, textureIdx = 0, normalIdx = 0;
                // since we know indices are not negative values we are parsing like this
                while (IsNumber(*curr))
                    positionIdx = 10 * positionIdx + (*curr++ - '0'); 
                curr++; // last ptr++ for jump '/'
                
                if (*curr != '/') // texture coord might not exist
                {
                    while (IsNumber(*curr)) 
                        textureIdx = 10 * textureIdx + (*curr++ - '0');
                    curr++; // last ptr++ for jump '/'
                }
                
                while (IsNumber(*curr))
                    normalIdx = 10 * normalIdx + (*curr++ - '0');
                curr++;
               
                positionIdx--, textureIdx--, normalIdx--;// obj index always starts from 1

                uint64_t hash = HashVertex(positionIdx, textureIdx, normalIdx);
                KeyValuePair<uint64_t, uint32_t>* find = uniqueVertices.Find(hash);

                if (find != uniqueVertices.end())
                {
                    indices.Add(find->value);
                    continue;
                }
                
                indices.Add(currentIndex);
                uniqueVertices.Insert(hash, currentIndex++);
                vertices.AddUninitialized(1);

                OBJVertex& vertex = vertices.Back();
                SmallMemCpy(vertex.pos, positions.Data() + (3 * positionIdx), sizeof(float) * 3);
                
                if (textureIdx != -1) // texture coord might not exist
                    SmallMemCpy(vertex.texCoord, texCoords.Data() + (2 * textureIdx) , sizeof(float) * 2);
                
                SmallMemCpy(vertex.normal, normals.Data() + (3 * normalIdx), sizeof(float) * 3);
            }

            if (IsNumber(*curr)) // is face has more than 3 element (for example color)
                while (*curr && !IsWhitespace(*curr) || *curr != '\n')
                    curr++;

            // skip line&whiteSpace
            while (*curr && (IsWhitespace(*curr) || *curr == '\n'))
                curr++;
        }
        
        // any vertex parsed?
        if (vertexStart != vertices.Size())
        {
            meshes.AddUninitialized(1);
        
            APrimitive primitive{};
            primitive.indices     = indices.Data() + indiceStart; // indicate linear indices
            primitive.vertices    = vertices.Data() + vertexStart; // we will point this vertex to vertices buffer
            primitive.attributes  = attributes; // pos, tex, normal bits.
            primitive.indexType   = 0x1405 - 0x1400; // GL_UNSIGNED_INT - GL_BYTE to get index so result is GraphicType_UnsignedInt
            primitive.numIndices  = indices.Size() - indiceStart;
            primitive.numVertices = vertices.Size() - vertexStart;
            primitive.material    = materialIndex;
            
            meshes.Back().numPrimitives = 1;
            meshes.Back().primitives = nullptr;

            SBPush(meshes.Back().primitives, primitive);
        }

        // 'o' is name, 'm' mtllib
        if (*curr == 'o' || *curr == 'm' || *curr == 's' || *curr == 'g')
        {
            while (*curr && *curr++ != '\n'); // skip line, header is unknown|unused
        }
    }

    scene->allVertices  = vertices.TakeOwnership();
    scene->allIndices   = indices.TakeOwnership();

    scene->numMeshes    = meshes.Size();     scene->meshes    = meshes.TakeOwnership();   
    scene->numMaterials = materials.Size();  scene->materials = materials.TakeOwnership();
    scene->numImages    = images.Size();     scene->images    = images.TakeOwnership();   
    
    scene->error = AError_NONE;
    return 1;
}

void FreeParsedObj(ParsedObj* obj)
{
    // also controls if arrays null or not
    for (int i = 0; i < obj->numMeshes; i++)
        SBFree(obj->meshes[i].primitives);

    if (obj->meshes)      delete[] obj->meshes;
    if (obj->materials)   delete[] obj->materials;
    if (obj->images)      delete[] obj->images;
    if (obj->allVertices) delete[] (OBJVertex*)obj->allVertices;
    if (obj->allIndices)  delete[] (uint32_t*)obj->allIndices;

    if (obj->materialText) FreeAllText(obj->materialText);
}