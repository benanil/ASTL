
/********************************************************************************
*    Purpose: Reading from files and Writing to files                           *
*    Author : Anilcan Gulkaya 2023 anilcangulkaya7@gmail.com github @benanil    *
********************************************************************************/

#pragma once

// AFile  AFileOpen(const char* fileName, AOpenFlag flag);
// void   AFileRead(void* dst, uint64_t size, AFile file);
// void   AFileWrite(const void* src, uint64_t size, AFile file);
// void   AFileClose(AFile file);
// bool   AFileExist(AFile file);
// uint64 AFileSize(AFile file);

// char*  GetFileExtension(path, size);
// bool   FileHasExtension(path, size, extension);
// char*  PathGoBackwards (path, end, skipSeparator);
// int    CopyFilename(char* out, char* path, int end); // < returns number of characters in filename
// bool   FileExist(file);
// uint64 FileSize(file);
// bool   RenameFile(oldFile, newFile);
// bool   MoveFile(oldLocation, newLocation);
// bool   DeleteFile(file);
// char*  ReadAllFile(fileName, buffer, numCharacters);
// void   CopyFile(source, dst, buffer);
// void   VisitFolder(path, visitFn);
// void   GetCurrentDirectory(buffer, bufferSize);
// void   AbsolutePath(path, outBuffer, bufferSize);
// void   CombinePaths(dst, a, b);
// bool   HasAnySubdir(path);

#include "Memory.hpp"
#include <stdio.h>
#include <stdint.h>
#include <sys/stat.h>

#ifdef _WIN32
    #include <io.h>
    #include <direct.h> 
    #include "dirent.h"
    #define F_OK 0
#else 
    #include <unistd.h>
    #include <sys/types.h>
    #include <fcntl.h>
    #include <dirent.h>
    #define _mkdir mkdir
    #define _fileno fileno
    #define _filelengthi64 filelength
#endif

#ifdef __ANDROID__
#include <android/asset_manager.h>
#include <game-activity/native_app_glue/android_native_app_glue.h>

extern "C" android_app* g_android_app;
#endif

#ifdef _WIN32
constexpr char ASTL_FILE_SEPERATOR = '\\';
#else
constexpr char ASTL_FILE_SEPERATOR = '/';
#endif

// path must be null terminated string
inline const char* GetFileExtension(const char* path, int size)
{
    while (path[size-1] != '.' && size > 0)
        size--;
    return path + size;
}

inline void ChangeExtension(char* path, int pathLen, const char* newExt)
{
    int lastDot = pathLen - 1;
    while (path[lastDot - 1] != '.')
        lastDot--;

    int i = lastDot;
    for (; *newExt; i++)
        path[i] = *newExt++;
    // clean the right with zeros
    while (i < pathLen)
        path[i++] = '\0';
}

inline bool FileHasExtension(const char* path, int size, const char* extension)
{
    int extLen = 0;
    // go to end of extension
    while (extension[1] != 0)
        extension++, extLen++;
    
    if (size <= extLen) return false;
    
    for (int i = 0; i <= extLen; i++)
    {
        if (*extension-- != path[--size])
            return false;
    }
    return true;
}

// returns pointer to the end of the new path
inline char* PathGoBackwards(char* path, int end, bool skipSeparator)
{
    if (path[end-1] == '/' || path[end-1] == '\\')
        path[--end] = '\0';
        
    while (end >= 0 && (path[end] != '/' && path[end] != '\\')) {
        path[end--] = '\0'; // Null-terminate the string.
    }
    
    if (skipSeparator && end >= 0) {
        path[end--] = '\0'; // Null-terminate again to remove the separator.
    }

    return path + end + 1; // Return the new starting point of the path.
}

// returns length of filename
inline int CopyFilename(char* out, const char* path, int end)
{
    int numChars = 0;
    while (end >= 0 && (path[end] != '/' && path[end] != '\\')) {
        out[numChars++] = path[end--]; // Null-terminate the string.
    }
    
    // reverse
    for (int i = 0, j = numChars-1; i < j; i++, j--) {
        char t = out[i];
        out[i] = out[j];
        out[j] = t;
    }

    return numChars; // Return the new starting point of the path.
}

// these functions works fine with file and folders
inline bool FileExist(const char* file)
{
#ifdef __ANDROID__
    AAsset* asset = AAssetManager_open(g_android_app->activity->assetManager, file, 0);
    if (asset == nullptr) {
        return false;
    } else {
        AAsset_close(asset);
        return true;
    }
#elif defined(_WIN32)
    return _access(file, F_OK) == 0;
#else
    return access(file, F_OK) == 0;
#endif
}

inline uint64_t FileSize(const char* file)
{
#ifdef __ANDROID__
    AAsset* asset = AAssetManager_open(g_android_app->activity->assetManager, file, 0);
    if (asset != nullptr) {
        off64_t sz = AAsset_getLength64(asset);
        AAsset_close(asset);
        return sz;
    }
    return 0;
#elif defined(_WIN32)
    struct stat sb;
    if (stat(file, &sb) == 0) return 0;
    return sb.st_size;
#else
    if (fseek(file.file, 0, SEEK_END) != 0) 
        return 0; // Or handle the error as appropriate

    long fileSize = ftell(file.file);
    if (fileSize == -1) 
        return 0; // Or handle the error as appropriate
    return (uint64_t)fileSize;
#endif
}

inline bool RenameFile(const char* oldFile, const char* newFile)
{
    return rename(oldFile, newFile) != 0;
}

inline bool CreateFolder(const char* folderName) 
{
    return _mkdir(folderName
               #ifndef _WIN32
               , 0777
               #endif 
               ) == 0;
}

inline bool IsDirectory(const char* path)
{
  struct stat file_info; 
  return stat(path, &file_info) == 0 && (S_ISDIR(file_info.st_mode));
}

enum AOpenFlag_
{
    AOpenFlag_ReadBinary, 
    AOpenFlag_WriteBinary,
    AOpenFlag_ReadText, 
    AOpenFlag_WriteText
};
typedef int AOpenFlag;

#ifdef __ANDROID__
struct AFile {
    AAsset* asset;
};

inline AFile AFileOpen(const char* fileName, AOpenFlag flag) {
    return { AAssetManager_open(g_android_app->activity->assetManager, fileName, 0) };
}

inline void AFileRead(void* dst, uint64_t size, AFile file, int alignment = 0) {
    AAsset_read(file.asset, dst, size);
}

inline void AFileSeekBegin(AFile file) {
    AAsset_seek(file.asset, 0, SEEK_SET);
}

inline void AFileSeek(long offset, AFile file) {
    AAsset_seek(file.asset, offset, SEEK_CUR);
}

inline void AFileWrite(const void* src, uint64_t size, AFile file, int alignment = 0)
{ }

inline void AFileClose(AFile file) {
    AAsset_close(file.asset);
}

inline bool AFileExist(AFile file) {
    return file.asset != nullptr;
}

inline uint64_t AFileSize(AFile file) {
    return AAsset_getLength(file.asset);
}

#else

struct AFile {
    FILE* file;
};

inline AFile AFileOpen(const char* fileName, AOpenFlag flag)
{
    FILE* file;
#ifdef _MSC_VER
    const char* modes[4] = { "rb", "wb", "rb", "wb" };
    fopen_s(&file, fileName, modes[flag]);
#else
    const char* modes[4] = { "rb", "wb", "r", "w"};
    file = fopen(fileName, modes[flag]);
    #endif

    if (flag == AOpenFlag_WriteText) {
        // Write the UTF-8 BOM
        unsigned char bom[] = {0xEF, 0xBB, 0xBF};
        fwrite(bom, sizeof(char), sizeof(bom), file);
    }

    AFile afile;
    afile.file = file;
    return afile;
}

inline void AFileRead(void* dst, uint64_t size, AFile file, int alignment = 1) {
    fread(dst, alignment, size, file.file);
}

inline void AFileWrite(const void* src, uint64_t size, AFile file, int alignment = 1) { 
    fwrite(src, alignment, size, file.file);
}

inline void AFileSeekBegin(AFile file) {
    fseek(file.file, 0, SEEK_SET);
}

inline void AFileSeek(long offset, AFile file) {
    fseek(file.file, offset, SEEK_CUR);
}

inline void AFileClose(AFile file) {
    fclose(file.file);
}

inline bool AFileExist(AFile file) { 
    return file.file != nullptr;
}

inline uint64_t AFileSize(AFile file)
{
#ifdef _WIN32
    return _filelengthi64(_fileno(file.file));
#elif defined(__ANDROID__)
    if (file.asset == nullptr) return 0;
    return AAsset_getLength(file.asset);
#else
    struct stat sb; stat(file.file, &sb);
    return sb.st_size;
#endif  
}

#endif

inline char* ReadAllFile(const char* fileName, char* buffer = 0)
{
    AFile file = AFileOpen(fileName, AOpenFlag_ReadBinary);
    if (!AFileExist(file))
        return nullptr;

    uint64_t fileSize = AFileSize(file);

    if (buffer == nullptr) {
        buffer = new char[fileSize+1]{}; // +1 for null terminator
    }

    AFileRead(buffer, fileSize, file);
    AFileClose(file);
    return buffer;
}

// don't forget to free using FreeAllText
// fileName      : path of the file that we want to load
// buffer        : is pre allocated memory if exist. otherwise null
// numCharacters : if not null returns length of the imported string
// startText     : if its not null will be added to start of the buffer
// note: if you define it you are responsible of deleting the buffer
inline char* ReadAllText(const char* fileName, char* buffer = 0, uint64_t* numCharacters = 0, const char* startText = 0)
{
    int startTextLen = 0;
    if (startText) 
        while (startText[startTextLen]) startTextLen++;

    AFile file = AFileOpen(fileName, AOpenFlag_ReadBinary);
    
    if (!AFileExist(file)) {
        perror("Error opening the file");
        return nullptr; // Return an error code
    }
    
    // Determine the file size
    uint64_t file_size = AFileSize(file);

    // Allocate memory to store the entire file
    if (buffer == nullptr) 
        buffer = new char[file_size + 40 + startTextLen] {}; // +1 for null terminator
    
    if (file_size >= 3)
    {
        AFileRead(buffer, 3, file);
        bool isBOM = buffer[0] == '\xEF' && buffer[1] == '\xBB' && buffer[2] == '\xBF';
        if (!isBOM) {
            AFileSeekBegin(file);
        }
    }

    if (startText) 
        while (*startText) *buffer++ = *startText++;
    
    if (buffer == NULL) {
        AFileClose(file);
        return nullptr; // Return an error code
    }

    // Read the entire file into the buffer
    AFileRead(buffer, file_size, file);
    buffer[file_size] = '\0'; // Null-terminate the buffer
    AFileClose(file);
    if (numCharacters) 
        *numCharacters = (uint64_t)file_size + 1;
    return buffer - startTextLen;
}

inline void FreeAllText(char* text)
{
    delete[] text;
}

inline void WriteAllBytes(const char *filename, const char *bytes, unsigned long size) 
{
    // Open the file for writing in binary mode
    AFile file = AFileOpen(filename, AOpenFlag_WriteBinary);
    if (!AFileExist(file)) 
    {
        perror("Failed to open file for writing");
        return;
    }

    AFileWrite(bytes, size, file);

    AFileClose(file);
}

struct ScopedText
{
    char* text;
    ScopedText(char* txt) : text(txt) {}
    ~ScopedText() { delete[] text; }
    operator char*() const { return text; }
};

// buffer is pre allocated memory if exist. otherwise null. 
// note: if you define it you are responsible of deleting the buffer
inline void CopyFile(const char* source, const char* dst, char* buffer = 0)
{
    uint64_t sourceSize = 0;
    bool bufferProvided = buffer != 0;
    char* sourceFile = ReadAllText(source, buffer, &sourceSize);
    
    AFile dstFile = AFileOpen(dst, AOpenFlag_WriteBinary);
    AFileWrite(sourceFile, sourceSize, dstFile);
    AFileClose(dstFile);

    if (!bufferProvided) delete[] buffer;
}

#ifndef _WIN32
    #define GetCurrentDirectory(size, outPath) getcwd(outPath, size)
#endif

// input: ../Textures/Tree.png
// output: C:/Source/Repos/Textures/Tree.png
inline void AbsolutePath(const char* path, char* outBuffer, int bufferSize)
{
    GetCurrentDirectory(bufferSize, outBuffer);
    int currLength = 0;

    const char* curr = outBuffer;
    while (*curr++) currLength++; // strlen

    while (*path)
    {
        if (path[0] == '.' && path[1] == '.')
        {
            const char* before = outBuffer + currLength;
            const char* newEnd = PathGoBackwards(outBuffer, currLength, true);
            currLength -= (int)(before - newEnd);
            path += 3; // skip two dot and seperator
            outBuffer[currLength++] = ASTL_FILE_SEPERATOR;
        }
        else 
        {
            while (*path && *path != '\\' && *path != '/')
                outBuffer[currLength++] = *path++;
            outBuffer[currLength++] = *path++;
        }
    }
    outBuffer[currLength] = '\n';
}

inline void CombinePaths(char* dst, const char* a, const char* b)
{
    int parentLen = StringLength(a);
    SmallMemCpy(dst, a, parentLen);
    dst[parentLen] = ASTL_FILE_SEPERATOR;
    SmallMemCpy(dst + parentLen + 1, b, StringLength(b));
}

// data is user defined, whatever data you need, and nullable
typedef void(*FolderVisitFn)(const char* path, void* data);

// thanks to https://github.com/tronkko/dirent this is windows and linux compatible
inline bool VisitFolder(const char* path, FolderVisitFn visitFn, void* data)
{
    DIR* dir;
    dirent* ent;
    if ((dir = opendir(path)) != NULL) 
    {
        char combined[512];
        while ((ent = readdir(dir)) != NULL)  
        {
            MemsetZero(combined, 512);
            CombinePaths(combined, path, ent->d_name);
            if (ent->d_name[0] != '.') visitFn(combined, data);
        }
        closedir(dir);
        return true;
    }
    return false;
}

inline bool HasAnySubdir(const char* path)
{
    DIR* dir;
    dirent* ent;
    if ((dir = opendir(path)) != NULL) 
    {
        char combined[512];
        int pathLen = StringLength(path);
        while ((ent = readdir(dir)) != NULL)  
        {
            MemsetZero(combined, 512);
            ASSERTR(pathLen + StringLength(ent->d_name) < 511, continue);
            CombinePaths(combined, path, ent->d_name);
            if (ent->d_name[0] != '.' && IsDirectory(combined)) return true;
        }
        closedir(dir);
    }
    return false;
}
