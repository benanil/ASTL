#pragma once

// char*  GetFileExtension(path, size);
// bool   FileHasExtension(path, size, extension);
// char*  PathGoBackwards (path, end, skipSeparator);
// bool   FileExist(file);
// uint64 FileSize(file);
// bool   RenameFile(oldFile, newFile);
// bool   MoveFile(oldLocation, newLocation);
// bool   DeleteFile(file);
// char*  ReadAllFile(fileName, buffer, numCharacters);
// void   CopyFile(source, dst, buffer);
// void   VisitFolder(path, pathLen, visitFn);
// void   GetCurrentDirectory(buffer, bufferSize);
// void   AbsolutePath(path, outBuffer, bufferSize);

#include <stdio.h>
#include <sys/stat.h>
#include <memory.h> // malloc free
#include <direct.h> // mkdir

#include "Algorithms.hpp"
#include "IntFltTypesLimits.hpp"

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

// these functions works fine with file and folders
inline bool FileExist(const char* file)
{
    struct stat sb;
    return stat(file, &sb) != 0;
}

inline uint64_t FileSize(const char* file)
{
    struct stat sb;
    if (stat(file, &sb) == 0) return 0;
    return sb.st_size;
}

inline bool RenameFile(const char* oldFile, const char* newFile)
{
    return rename(oldFile, newFile) != 0;
}

inline bool CreateFolder(const char* folderName) {
  return mkdir(folderName
               #ifndef _WIN32
               , 0777
               #endif 
               ) == 0;
}

inline bool IsDirectory(const char* path)
{
  struct stat file_info; 
  return stat(path, &file_info) == 0 && (file_info.st_mode & _S_IFDIR);
}

// don't forget to free
// buffer is pre allocated memory if exist. otherwise null
// note: if you define it you are responsible of deleting the buffer
inline char* ReadAllFile(const char* fileName, char* buffer = 0, int* numCharacters = 0)
{
#ifdef __ANDROID__
    AAsset* asset = AAssetManager_open(g_android_app->activity->assetManager, fileName, 0);
    off_t size = AAsset_getLength(asset);

    // Allocate memory to store the entire file
    if (buffer == nullptr) buffer = (char*)malloc(size + 1); // +1 for null terminat
    if (numCharacters) *numCharacters = size + 1;

    AAsset_read(asset, buffer, size);
    AAsset_close(asset);
    return buffer;
#else
    // Open the file for reading
    FILE* file = fopen(fileName, "rb");

    if (file == NULL) {
        perror("Error opening the file");
        return nullptr; // Return an error code
    }

    // Determine the file size
    fseek(file, 0, SEEK_END);
    long file_size = ftell(file);
    rewind(file);

    // Allocate memory to store the entire file
    if (buffer == nullptr) buffer = (char*)malloc(file_size + 1); // +1 for null terminator

    if (buffer == NULL) {
        fclose(file);
        return nullptr; // Return an error code
    }

    // Read the entire file into the buffer
    fread(buffer, 1, file_size, file);
    buffer[file_size] = '\0'; // Null-terminate the buffer
    fclose(file);
    if (numCharacters) *numCharacters = file_size + 1;
    return buffer;
#endif
}

inline void FreeAllText(char* text)
{
    free(text);
}

struct ScopedText
{
    char* text;
    ScopedText(char* txt) : text(txt) {}
   ~ScopedText() { free(text); }
};

// buffer is pre allocated memory if exist. otherwise null. 
// note: if you define it you are responsible of deleting the buffer
inline void CopyFile(const char* source, const char* dst, char* buffer = 0)
{
    int sourceSize = 0;
    bool bufferProvided = buffer != 0;
    char* sourceFile = ReadAllFile(source, buffer, &sourceSize);
    FILE* dstFile = fopen(dst, "w");
    fwrite(sourceFile, 1, sourceSize, dstFile);
    fclose(dstFile);

    if (!bufferProvided) free(buffer);
}

typedef void(*FolderVisitFn)(char* buffer, int bufferLength, const char* fileName, bool isFolder, uint64_t fileSize);

#if _WIN32

// #define WIN32_LEAN_AND_MEAN 
// #define NOMINMAX
// #include <Windows.h>

// if windows.h is not included forward declare the functions instead of including windows.h header
#ifndef _MINWINBASE_

#define INVALID_HANDLE_VALUE ((void*)(long)-1)
#define MAX_PATH                 260
#define FILE_ATTRIBUTE_DIRECTORY 0x00000010  

typedef struct _FILETIME {
    uint64_t dwLowDateTime;
    uint64_t dwHighDateTime;
} FILETIME;

typedef struct _WIN32_FIND_DATAA {
    uint64_t dwFileAttributes;
    FILETIME ftCreationTime;
    FILETIME ftLastAccessTime;
    FILETIME ftLastWriteTime;
    uint64_t nFileSizeHigh;
    uint64_t nFileSizeLow;
    uint64_t dwReserved0;
    uint64_t dwReserved1;
    char     cFileName[MAX_PATH];
    char     cAlternateFileName[14];
    uint64_t dwFileType; // Obsolete. Do not use.
    uint64_t dwCreatorType; // Obsolete. Do not use
    uint16_t wFinderFlags; // Obsolete. Do not use
} WIN32_FIND_DATA, *PWIN32_FIND_DATA, *LPWIN32_FIND_DATA;

extern "C"
{
    __declspec(dllimport) void* __cdecl FindFirstFileA(const char* lpFileName, LPWIN32_FIND_DATA lpFindFileData);

    __declspec(dllimport) int __cdecl FindNextFileA(void* hFindFile, LPWIN32_FIND_DATA lpFindFileData);

    __declspec(dllimport) int __cdecl FindClose(void* hFindFile);

#define GetCurrentDirectory(bufferLength, buffer) GetCurrentDirectoryA(bufferLength, buffer)
#ifndef _PROCESSENV_
    __declspec(dllimport) uint64_t GetCurrentDirectoryA(uint64_t nBufferLength, char* lpBuffer);
#endif
}
#endif // windows.h included

// visit all files and folders in the path
inline void VisitFolder(char* path, int pathLen, FolderVisitFn visitFn)
{
    WIN32_FIND_DATA findFileData;
    // we should add \\* to end the of the buffer
    path[pathLen]     = '\\';
    path[pathLen + 1] = '*';
    pathLen++; // skip '\\'
    void* hFind = FindFirstFileA(path, &findFileData);

    if (hFind == INVALID_HANDLE_VALUE) {
        printf("Error opening directory: %s \n", path);
        return;
    }

    if (path[pathLen-1] == '*')
        path[--pathLen] = 0;

    do 
    {
        bool isFolder = findFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY;
        visitFn(path, pathLen, findFileData.cFileName, isFolder, findFileData.nFileSizeLow);
    } while (FindNextFileA(hFind, &findFileData) != 0);
 
    FindClose(hFind);
}

#else
#include <sys/types.h>
#include <dirent.h>
#include "String.hpp"

#define MAX_PATH 260

inline uint64_t GetCurrentDirectory(char* buffer, uint64_t bufferSize)
{
  ASSERT(getcwd(buffer, bufferSize));
  return StringLength(buffer);
}

inline void VisitFolder(char *path, int pathLen, FolderVisitFn visitFn) 
{
    struct dirent *entry;
    struct stat fileStat;
    DIR *directory = opendir(path);

    if (directory == NULL) {
        perror("Error opening directory");
        return;
    }

    while ((entry = readdir(directory)) != NULL) 
    {
        if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0) {
            continue; // Skip "." and ".." entries.
        }

        // Build the full path to the file or folder.
        char filePath[pathLen + strlen(entry->d_name) + 2]; // +2 for '/' and '\0'
        snprintf(filePath, sizeof(filePath), "%s/%s", path, entry->d_name);

        if (stat(filePath, &fileStat) == 0) {
            bool isFolder = S_ISDIR(fileStat.st_mode);
            off_t fileSize = fileStat.st_size;
            visitFn(filePath, strlen(filePath), entry->d_name, isFolder, fileSize);
        } else {
            perror("Error getting file stat");
        }
    }

    closedir(directory);
}
#endif // 

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