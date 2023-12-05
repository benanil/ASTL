
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
#include "../Random.hpp"
#include "../Array.hpp"

inline void ChangeExtension(char* path, const char* newExt, size_t len)
{
    path[len - 1] = newExt[2]; path[len - 2] = newExt[1]; path[len - 3] = newExt[0];
}

static const char* ParseFloat16(const char*& curr, short& flt)
{
    flt = (short)(ParseFloat(curr) * 1000.0f);
    return curr;
}

void ParseObj(const char* path, ParsedScene* scene)
{
    if (!FileExist(path)) {
        scene->error = AError_BUFFER_PARSE_FAIL;
        return;
    }

    int pathLen = StringLength(path);
    FixedSizeGrowableAllocator<char> stringAllocator(2048);
    
    char* pathDup = stringAllocator.Allocate(pathLen + 8);
    SmallMemCpy(pathDup, path, pathLen);

    const uint64 sz = FileSize(path);
    char* objText = ReadAllFile(path);
    
#if defined(DEBUG) || defined(_DEBUG)
    // ascii utf8 support check
    if (IsUTF8ASCII(objText, sz) != 1)  { scene->error = AError_NON_UTF8; return result; } 
#endif
    ChangeExtension(pathDup, "mtl", pathLen);
    char* mtlPath = pathDup; // for readibility
    uint msz = FileExist(mtlPath) ? (uint)FileSize(mtlPath) : 0u;
    
    char* mtlText = nullptr;
    // read if material path exist 
    if (msz) {
        char* buffer = stringAllocator.Allocate(msz + 8);
        mtlText = ReadAllFile(mtlPath, buffer);
#if defined(DEBUG) || defined(_DEBUG)
        if (IsUTF8ASCII(mtlText, msz) != 1) { scene->error = AError_NON_UTF8; return result; }
#else
    }

    Array<float> positions, indices, texCoords, normals; 
    Array<AMesh> meshes;

    // hashMap for material indices, materialMap[materialNameHash & 255] = material Index 
    unsigned char materialMap[512] = {0};
	
    // import materials
    char* curr = mtlText, *currEnd = curr + msz;
    Array<AMaterial> materials;
    Array<AImage> images; // index to start of mtltext
#if 0
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
            materials.Back().specularColor = ~0u, currMaterial->diffuseColor = ~0u; // white
            materials.Back().specularFactor = 2.2f * 400.0f;
            materials.Back().roughness = 0.6f * 400.0f;
            materials.Back().name = PointerDistance(mtlText, curr);
            
            unsigned hash = WangHash(uint(curr[0]) | uint(curr[1] << 8) | uint(curr[2] << 16));
            while(*curr != '\n' && !IsWhitespace(*curr))
                hash = *curr++ + (hash << 6) + (hash << 16) - hash; 
            *curr++ = '\0'; // null terminator
            if (materialMap[hash & 511])
            { scene->error = AError_HASH_COLISSION; return;}

            materialMap[hash & 511] = scene->numMaterials++;
        }
        else if (curr[0] == 'N' && curr[1] == 's') 
        {
            // shininess
            curr = ParseFloat16(curr, materials.Back().specularFactor); 
        }
        else if (curr[0] == 'd') 
        {
            // roughness
            curr = ParseFloat16(curr, materials.Back().roughnessFactor);
        }
        else if (curr[0] == 'K' && curr[1] == 'd') 
        { 
            // diffuse color
            float colorf[3] = {ParseFloat(curr), ParseFloat(curr), ParseFloat(curr)};
            materials.Back().baseColorFactor = PackColorRGBU32(colorf);
        }
        else if (curr[0] == 'K' && curr[1] == 's') { // specular color
            float colorf[3] = {ParseFloat(curr), ParseFloat(curr), ParseFloat(curr)};
            materials.Back().specularColor = PackColorRGBU32(colorf);
        }
        else if (curr[0] == 'm') // map_bla
        {
            if (curr[4] == 'K' && curr[5] == 'd') {
                // set this pointer as texture path data
                materials.Back().baseColorTexture.index = images.Size();
                AImage image{ mtlText + PointerDistance(mtlText, (curr += 7)};
                images.Add(image);
                // find end point of path
                while (*curr != '\n' && *curr != '\r') curr++;
                *curr++ = '\0'; // addd null terminator
            }
            else if (curr[4] == 'K' && curr[5] == 's') {
                materials.Back().specularTexture.index = images.Size();
                AImage image{ mtlText + PointerDistance(mtlText, (curr += 7)};
                images.Add(image);

                while (*curr != '\n' && *curr != '\r') curr++; // find end point of path
                *curr++ = '\0'; // addd null terminator
            }
            else while (*curr != '\n') curr++;
        } // skip line. header is unknown | unused
        else while (*curr != '\n' && *curr) curr++;
    }

    curr = objText, currEnd = curr + sz;
    unsigned currentMaterial = 0;

    struct Vertex { float pos[3]; float texCoord[2]; float normal[3]; };
    Array<Vertex> vertices;

    while (curr && *curr && curr < currEnd)
    {
        if (*curr == '#')  while (*curr++ != '\n');
        while (*curr == '\n' || IsWhitespace(*curr)) curr++;

        if (*curr == 'v')
        {
            while (curr[1] == ' ') { // vertex=position 
                curr += 2; 
                positions.Add(ParseFloat(curr)); 
                positions.Add(ParseFloat(curr)); 
                positions.Add(ParseFloat(curr));
                // skip line&whiteSpace
                while (*curr == '\n' || IsWhitespace(*curr)) curr++;
            }

            while (curr[1] == 't') {
                curr += 2;
                texCoords.Add(ParseFloat(curr)); 
                texCoords.Add(ParseFloat(curr)); 
                // skip line&whiteSpace
                while (*curr == '\n' || IsWhitespace(*curr)) curr++;
            }

            while (curr[1] == 'n') {
                curr += 2; 
                normals.Add(ParseFloat(curr));
                normals.Add(ParseFloat(curr)); 
                normals.Add(ParseFloat(curr)); 
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
            materialIndex = materialMap[hash & 511]; // use material index for this group of triangles
        }
        
        int vertexStart = vertices.Size();
        while (curr[0] == 'f')
        {
            curr += 2;
      
            for(int i = 0; i < 3; ++i) // compiler please unroll :D
            {
                int positionIdx = 0, textureIdx = 0, normalIdx = 0;
                // since we know indices are not negative values we are parsing like this
                while (IsNumber(*curr)) positionIdx = 10 * positionIdx + (*curr++ - '0'); curr++; // last ptr++ for jump '/'
                while (IsNumber(*curr)) textureIdx = 10 * textureIdx + (*curr++ - '0'); curr++;
                while (IsNumber(*curr)) normalIdx = 10 * normalIdx + (*curr++ - '0'); curr++;
                vertices.AddUninitialized(1);
                Vertex& vertex = vertices.Back();
                SmallMemCpy(vertex.pos     , positions.Data() + (3 * positionsIdx), sizeof(float) * 3);
                SmallMemCpy(vertex.texCoord, texCoords.Data() + (2 * textureIdx)  , sizeof(float) * 2);
                SmallMemCpy(vertex.normal  , normals.Data()   + (3 * positionsIdx), sizeof(float) * 3);
            }

            // skip line&whiteSpace
            while (IsWhitespace(*curr) || *curr == '\n') curr++;
        }
        meshes.AddUninitialized(1);

        APrimitive primitive{};
        primitive.indices     = (void*)0XDEADBEAF; // indicate linear indices
        primitive.vertices    = (void*)vertexStart; // we will point this vertex to vertices buffer
        primitive.attributes  = 1 | 2 | 4; // pos, tex, normal bits.
        primitive.indexType   = 0x1405 - 0x1400; // GL_UNSIGNED_INT - GL_BYTE to get index
        primitive.numIndices  = vertices.Size() - vertexStart;
        primitive.numVertices = primitive.numIndices;
        primitive.material    = materialIndex;

        SBPush(meshes.Back().primitives, primitive);

        if (*curr == 'o' || *curr == 'm' || *curr == 's')  while (*curr++ != '\n'); // skip line, header is unknown|unused
    }

    scene->numNodes = meshes.Size();
    scene->nodes = new ANode[meshes.Size()];
    MemsetZero(scene->nodes, sizeof(ANode) * scene->numNodes);

    for (int i = 0; i < scene->numNodes; i++)
    {
        scene->nodes[i].rotation[3] = 1.0;
        scene->nodes[i].scale[0] = scene->nodes[i].scale[1] = scene->nodes[i].scale[2] = 1.0;
        scene->nodes[i].index = i;
    }

    scene->numMeshes    = meshes.Size();     scene->meshes    = Meshes.TakeOwnership();   
    scene->numMaterials = materials.Size();  scene->materials = Materials.TakeOwnership();
    scene->numImages    = images.Size();     scene->images    = Images.TakeOwnership();   
#endif
}