
#include "../Profiler.hpp"
#include "AdventOfCode2021.cpp"
#include "AdventOfCode2022.cpp"
#include <stdio.h>

#include <windows.h>
#include <string.h>
#include "../String.hpp"
#include "../Array.hpp"

bool HasExtension(const char* path, int size, const char* extension)
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

bool IsRealClass(char* curr)
{
    // return false if template
    if (curr[-1] == '<' || curr[-2] == 'm') // is enum class?
        return false;
    
    curr += 6; // skip class and a space
    if (*curr == 'T' && curr[1] == 'K') // skip TK_API
        curr += 6;

    while (IsWhitespace(*curr) || *curr == '\n')
        curr++;

    while (IsChar(*curr) || IsNumber(*curr))
        curr++;
        
    while (IsWhitespace(*curr) || *curr == '\n')
        curr++;

    if (*curr == ':' || *curr == '{')
        return true;

    return false;
}

char* FindNextClass(char* curr)
{
    while (*curr)
    {
        curr = FindCharInString(curr, 'c');
        if (curr == nullptr)
            return nullptr;
        
        if (StringEqual(curr, "class ", 5))
        {
            if (!IsRealClass(curr))
            {
                curr++;
                continue;
            }
            return curr;
        }
        curr++;
    }
    return nullptr;
}

char* PathGoBackwards(char* path, int end, bool skipSeparator)
{
    while (end >= 0 && (path[end-1] != '/' && path[end-1] != '\\')) {
        path[end--] = '\0'; // Null-terminate the string.
    }
    
    path[end--] = '\0'; // Null-terminate the string.

    if (skipSeparator && end >= 0) {
        path[end--] = '\0'; // Null-terminate again to remove the separator.
    }

    return path + end + 1; // Return the new starting point of the path.
}

#include <sys/stat.h>

bool FolderExists(const char* folder)
{
    struct stat sb;
    return stat(folder, &sb) == 0 && S_ISDIR(sb.st_mode);
}

int numberOfClasses = 0;
HashSet<uint64_t> hashMap{};
Array<char, MallocAllocator<char>> textBuffer{};

void SearchInFolder(char* buffer, int fileSize)
{
    printf("search directory: %s\n", buffer);

    WIN32_FIND_DATA findFileData;
    HANDLE hFind = FindFirstFile(buffer, &findFileData);

    if (hFind == INVALID_HANDLE_VALUE) {
        printf("Error opening directory\n");
        getchar();
        return;
    }

    if (buffer[fileSize-1] == '*')
        buffer[--fileSize] = 0;

    printf("searching \n");

    const int textMaxCharSize = 120 * 1000;

    if (textBuffer.Capacity() < textMaxCharSize)
    textBuffer.Resize(textMaxCharSize);

    do 
    {
        // is file?
        if (!(findFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) 
        {
            int filelength = StringLength(findFileData.cFileName);
            if (HasExtension(findFileData.cFileName, filelength, "hpp") || 
                HasExtension(findFileData.cFileName, filelength, ".h") )
            {
                // zero the path buffer
                SmallMemSet(buffer + fileSize, 0, MAX_PATH - fileSize);
                // append file to folder
                SmallMemCpy(buffer + fileSize, findFileData.cFileName, filelength);
            
                FILE* file = fopen(buffer, "r");
                printf("opening the file: %s\n", buffer);

                if (file == NULL) {
                    printf("cannot oppen file");
                    return;
                }

                // Determine the file size
                uint64_t fileSize = findFileData.nFileSizeLow;

                if (fileSize + 1 > textBuffer.Capacity())
                    textBuffer.Resize((int)fileSize);

                MemSet(textBuffer.Data(), 0, textBuffer.Capacity());
                fread(textBuffer.Data(), 1, fileSize, file);
                textBuffer[fileSize] = '\n';

                char* curr = textBuffer.Data();

                while (curr = FindNextClass(curr))
                {
                    curr += 6; // skip class and a space
                    if (*curr == 'T' && curr[1] == 'K') // skip TK_API
                        curr += 6;
                    
                    while (IsWhitespace(*curr))
                        curr++;

                    int nameLen = 0;
                    while (IsChar(*curr) || IsNumber(*curr))
                        nameLen++, curr++;

                    *curr = '\0'; // null terminate  
                    curr -= nameLen;
                    uint64_t hash = StringToHash(curr, nameLen);
                    printf("class name: %s, hash: %llu \n", curr, hash);

                    if (hashMap.Contains(hash))
                    {
                        printf("your class have collission: %s \n", curr);
                        getchar();
                        return;
                    }
                    hashMap.Insert(hash);
                    numberOfClasses++;
                    curr += nameLen + 1;
                }
                fclose(file);
            }
        }
    } while (FindNextFile(hFind, &findFileData) != 0);

    printf("Success no class collission founded in number of classes: %i\n", numberOfClasses);

    FindClose(hFind);
}

int main()
{   
    char buffer[MAX_PATH]{};
    GetCurrentDirectory(MAX_PATH, buffer);
    
    int fileSize = StringLength(buffer);
    buffer[fileSize]     = '\\';
    buffer[fileSize + 1] = '*';
    fileSize++; // skip '\\'

    SearchInFolder(buffer, fileSize);

    // try to search in editor folder
    // remove the last header files name from the path
    SmallMemSet(buffer + fileSize-1, 0, MAX_PATH - fileSize + 1);
    char* endPath = PathGoBackwards(buffer, fileSize, false);
    SmallMemCpy(endPath, "Editor\\", 8);
    
    printf("editor file: %s \n", buffer);

    // (for toolkit)
    if (FolderExists(buffer))
    {
        printf("editor file is exists\n");
        fileSize = StringLength(buffer);
        buffer[fileSize++] = '*';
        SearchInFolder(buffer, fileSize);
    }

    getchar();

    return 0;
    // AdventOfCodeTests(); // 2022
    BeginProfile();
    {
        TimeBlock("all");
        Day1();
        Day2();
        Day4();
    }
    EndAndPrintProfile();
    getchar();
    return 0;
}
  