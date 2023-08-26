
#include "../Profiler.hpp"
#include "AdventOfCode2021.cpp"
#include "AdventOfCode2022.cpp"

#include "../String.hpp"
#include "../Array.hpp"

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

int numberOfClasses = 0;
bool hashCheckFailed = false;
HashSet<uint64_t> hashMap{};
Array<char, MallocAllocator<char>> textBuffer{};

// search header files for class names and hash check
void SearchHeaderFiles(char* buffer, int pathSize, const char* fileName, bool isFolder, uint64_t fileSize)
{
    int fileLength = StringLength(fileName);
    
    if (!(FileHasExtension(fileName, fileLength, "hpp") || 
          FileHasExtension(fileName, fileLength, ".h")))
    {
        return;
    }
    // zero the path buffer
    SmallMemSet(buffer + pathSize, 0, MAX_PATH - pathSize);
    // append file to folder
    SmallMemCpy(buffer + pathSize, fileName, pathSize);
        
    FILE* file = fopen(buffer, "r");
    printf("opening file: %s\n", fileName);
    
    if (file == NULL) {
        printf("cannot oppen file\n");
        return;
    }
    
    if (fileSize + 1 > textBuffer.Capacity())
        textBuffer.Resize((int)fileSize);
    
    MemSet(textBuffer.Data(), 0, textBuffer.Capacity());
    // read header into text buffers
    fread(textBuffer.Data(), 1, fileSize, file);
    textBuffer[fileSize] = '\n';
    
    char* curr = textBuffer.Data();
    
    // search class names and hash them, then check for uniqueness
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
            hashCheckFailed = true;
            return;
        }
        hashMap.Insert(hash);
        numberOfClasses++;
        curr += nameLen + 1;
    }
    fclose(file);
}

void SearchInFolder(char* buffer, int fileSize)
{
    hashCheckFailed = false;
    printf("search directory: %s\n", buffer);

    const int textMaxCharSize = 120 * 1000;
    
    if (textBuffer.Capacity() < textMaxCharSize)
        textBuffer.Resize(textMaxCharSize);
    
    // from IO.hpp
    VisitFolder(buffer, fileSize, SearchHeaderFiles);
  
    if (!hashCheckFailed)
    {
        printf("Success no class collission founded in number of classes: %i\n", numberOfClasses);
    }
}

int main()
{   
    char absolute[MAX_PATH]{0};
    AbsolutePath("..\\..\\Games\\Full\\..", absolute, MAX_PATH);
    printf("absolute path: %s ", absolute);
    getchar();
    return 0;

    char buffer[MAX_PATH]{};
    GetCurrentDirectory(buffer, MAX_PATH);
    
    int fileSize = StringLength(buffer);
    SearchInFolder(buffer, fileSize);

    // try to search in editor folder
    // remove the last header files name from the path
    SmallMemSet(buffer + fileSize-1, 0, MAX_PATH - fileSize + 1);
    char* endPath = PathGoBackwards(buffer, fileSize, false);
    SmallMemCpy(endPath, "Editor", 8);
    
    printf("editor file: %s \n", buffer);

    // (for toolkit)
    if (FileExist(buffer))
    {
        printf("editor file is exists\n");
        SearchInFolder(buffer, StringLength(buffer));
    }

    getchar();
    return 0;
    AdventOfCodeTests(); // 2022
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
  