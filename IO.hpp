
/********************************************************************************
*    Purpose: Reading from files and Writing to files                           *
*    Author : Anilcan Gulkaya 2023 anilcangulkaya7@gmail.com github @benanil    *
********************************************************************************/

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

#include "Algorithms.hpp"
#include <stdio.h>
#include <sys/stat.h>

#if defined __WIN32__ || defined _WIN32 || defined _Windows
    #if !defined S_ISDIR
            #define S_ISDIR(m) (((m) & _S_IFDIR) == _S_IFDIR)
    #endif
#endif

#ifdef __ANDROID__
#include <android/asset_manager.h>
#include <game-activity/native_app_glue/android_native_app_glue.h>

extern android_app* g_android_app;
#else
#include <direct.h> // mkdir
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

struct ScopedFILE
{
    FILE* file;
    ScopedFILE(FILE* _file) : file(_file) {}
    ~ScopedFILE()
    {
        fclose(file); 
    }
};

#ifdef WIN32
#include <io.h>
#define F_OK 0
#define access _access
#else 
#include <unistd.h>
#endif

// these functions works fine with file and folders
inline bool FileExist(const char* file)
{
    return access(file, F_OK) == 0;
}

inline uint64_t FileSize(const char* file)
{
#ifndef _WIN32
    struct stat sb;
    if (stat(file, &sb) == 0) return 0;
    return sb.st_size;
#else
    ScopedFILE f = fopen(file, "rb");
    return _filelengthi64(fileno(f.file));
#endif
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
  return stat(path, &file_info) == 0 && (S_ISDIR(file_info.st_mode));
}

enum AOpenFlag_
{
    AOpenFlag_Read, 
    AOpenFlag_Write
};
typedef int AOpenFlag;

#ifdef __ANDROID__
struct AFile
{
    AAsset* asset;
};

inline AFile AFileOpen(const char* fileName, AOpenFlag flag)
{
    return { AAssetManager_open(g_android_app->activity->assetManager, fileName, 0) };
}

inline void AFileRead(void* dst, uint64_t size, AFile file)
{
    AAsset_read(file.asset, dst, size);
}

inline void AFileWrite(const void* src, uint64_t size, AFile file)
{ }

inline void AFileClose(AFile file)
{
    AAsset_close(file.asset);
}

inline bool AFileExist(AFile file)
{
    return file.asset != nullptr;
}

inline uint64_t AFileSize(AFile file)
{
    return AAsset_getLength(file.asset);
}

#else

struct AFile
{
    FILE* file;
};

inline AFile AFileOpen(const char* fileName, AOpenFlag flag)
{
    const char* modes[2] = { "rb", "wb" };
    return { fopen(fileName, modes[flag]) };
}

inline void AFileRead(void* dst, uint64_t size, AFile file)
{
    fread(dst, 1, size, file.file);
}

inline void AFileWrite(const void* src, uint64_t size, AFile file)
{ 
    fwrite(src, 1, size, file.file);
}

inline void AFileClose(AFile file)
{
    fclose(file.file);
}

inline bool AFileExist(AFile file)
{ 
    return file.file != nullptr;
}

inline uint64_t AFileSize(AFile file)
{
#ifdef _WIN32
    return _filelengthi64(fileno(file.file));
#else
    struct stat sb; stat(file.file, &sb);
    return sb.st_size;
#endif  
}

#endif

// don't forget to free using FreeAllText
// fileName      : path of the file that we want to load
// buffer        : is pre allocated memory if exist. otherwise null
// numCharacters : if not null returns length of the imported string
// startText     : if its not null will be added to start of the buffer
// note: if you define it you are responsible of deleting the buffer
inline char* ReadAllFile(const char* fileName, char* buffer = 0, long* numCharacters = 0, const char* startText = 0)
{
    int startTextLen = 0;
    if (startText) 
        while (startText[startTextLen]) startTextLen++;

    AFile file = AFileOpen(fileName, AOpenFlag_Read);
    
    if (!AFileExist(file)) {
        perror("Error opening the file");
        return nullptr; // Return an error code
    }
    
    // Determine the file size
    size_t file_size = AFileSize(file);

    // Allocate memory to store the entire file
    if (buffer == nullptr) 
        buffer = new char[file_size + 40 + startTextLen] {}; // +1 for null terminator
    
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
    if (numCharacters) *numCharacters = file_size + 1;
    return buffer - startTextLen;
}

inline void FreeAllText(char* text)
{
    delete[] text;
}

inline void WriteAllBytes(const char *filename, const char *bytes, unsigned long size) 
{
    // Open the file for writing in binary mode
    AFile file = AFileOpen(filename, AOpenFlag_Write);
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
};

// buffer is pre allocated memory if exist. otherwise null. 
// note: if you define it you are responsible of deleting the buffer
inline void CopyFile(const char* source, const char* dst, char* buffer = 0)
{
    long sourceSize = 0;
    bool bufferProvided = buffer != 0;
    char* sourceFile = ReadAllFile(source, buffer, &sourceSize);
    FILE* dstFile = fopen(dst, "w");
    fwrite(sourceFile, 1, sourceSize, dstFile);
    fclose(dstFile);

    if (!bufferProvided) delete[] buffer;
}

typedef void(*FolderVisitFn)(char* buffer, int bufferLength, const char* fileName, bool isFolder, unsigned long long fileSize);

#if _WIN32

// #define WIN32_LEAN_AND_MEAN 
// #define NOMINMAX
// #define VC_EXTRALEAN
// #include <Windows.h>

// if windows.h is not included forward declare the functions instead of including windows.h header
#ifndef _MINWINBASE_

#define INVALID_HANDLE_VALUE ((void*)(long)-1)
#define MAX_PATH                 260
#define FILE_ATTRIBUTE_DIRECTORY 0x00000010  

typedef struct _FILETIME {
    unsigned long long dwLowDateTime;
    unsigned long long dwHighDateTime;
} FILETIME;

typedef struct _WIN32_FIND_DATAA {
    unsigned long long dwFileAttributes;
    FILETIME ftCreationTime;
    FILETIME ftLastAccessTime;
    FILETIME ftLastWriteTime;
    unsigned long long nFileSizeHigh;
    unsigned long long nFileSizeLow;
    unsigned long long dwReserved0;
    unsigned long long dwReserved1;
    char     cFileName[MAX_PATH];
    char     cAlternateFileName[14];
    unsigned long long dwFileType; // Obsolete. Do not use.
    unsigned long long dwCreatorType; // Obsolete. Do not use
    unsigned short wFinderFlags; // Obsolete. Do not use
} WIN32_FIND_DATA, *PWIN32_FIND_DATA, *LPWIN32_FIND_DATA;

extern "C"
{
    __declspec(dllimport) void* __cdecl FindFirstFileA(const char* lpFileName, LPWIN32_FIND_DATA lpFindFileData);

    __declspec(dllimport) int __cdecl FindNextFileA(void* hFindFile, LPWIN32_FIND_DATA lpFindFileData);

    __declspec(dllimport) int __cdecl FindClose(void* hFindFile);

#define GetCurrentDirectory(bufferLength, buffer) GetCurrentDirectoryA(bufferLength, buffer)
#ifndef _PROCESSENV_
    __declspec(dllimport) unsigned long long GetCurrentDirectoryA(unsigned long long nBufferLength, char* lpBuffer);
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

inline unsigned long long GetCurrentDirectory(unsigned long long bufferSize, char* buffer)
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
        if (*entry->d_name == '.') {
            continue; // Skip "." and ".." entries.
        }

        // Build the full path to the file or folder.
        char filePath[256]; // +2 for '/' and '\0'
        snprintf(filePath, sizeof(filePath), "%s/%s", path, entry->d_name);

        if (stat(filePath, &fileStat) == 0) {
            bool isFolder = S_ISDIR(fileStat.st_mode);
            off_t fileSize = fileStat.st_size;
            int pathLen=0;
            while (filePath[pathLen]) pathLen++;
            visitFn(filePath, pathLen, entry->d_name, isFolder, fileSize);
        } else {
            perror("Error getting file stat");
        }
    }

    closedir(directory);
}
#endif // 

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