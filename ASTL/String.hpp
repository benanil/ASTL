
#pragma once

#include <stdlib.h>
#include <string.h> //!< todo remove standard headers
#include <memory.h> //!< memcpy and memset

#include "Common.hpp"
#include "Random.hpp"
#include "Algorithms.hpp"
#include "Memory.hpp"

template<typename AllocatorT = MallocAllocator<char>>
class String
{
public:
    int capacity = 0;
    int size     = 0;
    char* ptr    = nullptr;
    AllocatorT allocator{};
public:
    //todo allocator

    ~String() 
    {
        if (ptr != nullptr) 
        {
            free(ptr); 
            ptr = nullptr;  
            size = 0; 
            capacity = 0; 
        }
    }

    explicit String(int _size) : size(0), capacity(_size + 1)
    {
        ptr = (char*)calloc(capacity, sizeof(char));
    }
    // copy constructor
    String(const String& other) : size(other.size), capacity(other.capacity)
    {
        ptr = (char*)calloc(capacity, sizeof(char));
        MemCpy(ptr, other.ptr, size);
    }
    // move constructor 
    String(String&& other) : size(other.size), capacity(other.capacity), ptr(other.ptr)
    {
        other.ptr = nullptr;
    }

    String() : size(0), capacity(32)
    {
        ptr = (char*)calloc(capacity, 1);
    }

    String(const char* _ptr) : size(StrLen(_ptr))
    {
        capacity = size + 32;
        ptr = (char*)calloc(capacity, 1);
        MemCpy(ptr, _ptr, size);
    }

    String(char* _ptr) : size(StrLen(_ptr))
    {
        capacity = size + 32;
        ptr = (char*)calloc(capacity, 1);
        MemCpy(ptr, _ptr, size);
    }
    // asign operator
    String& operator = (const String& right) 
    {
        if (&right != this) 
        {
            CapacityCheck(right.size);
            MemCpy(ptr, right.ptr, right.size);
            size = right.size;
        }
        return *this;
    }

    bool  operator == (const String& b) const { return  Compare(*this, b); }
    bool  operator != (const String& b) const { return !Compare(*this, b); }
    bool  operator == (const char* b)   const { return !strcmp(ptr, b); }
    bool  operator != (const char* b)   const { return  strcmp(ptr, b); }
    char  operator[](int index)         const { ASSERT(index < size && index > 0); return ptr[index]; }
    char& operator[](int index)               { ASSERT(index < size && index > 0); return ptr[index]; }

    char* begin()              { return ptr;        }
    char* end()                { return ptr + size; }
    const char* cbegin() const { return ptr;        }
    const char* cend()   const { return ptr + size; }

    const char* CStr()  const  { return ptr; }
    const char* c_str() const  { return ptr; }
    const char* data()  const  { return ptr; }
          char* data()         { return ptr; }

    static FINLINE bool Compare(const String& a, const String& b)
    {
        return a.size == b.size && !strcmp(a.ptr, b.ptr);
    }

    void Reset()
    {
        MemSet(ptr, 0, size);
        size = 0;
    }

    void Set(const char* str)
    {
        Clear();
        CapacityCheck(StrLen(str));
        MemCpy(ptr, str, StrLen(str));
    }

    void Set(const String& str)
    {
        Clear();
        CapacityCheck(str.size);
        size = str.Length();
        MemCpy(ptr, str.CStr(), str.size);
    }

    void Reserve(int size)
    {
        if (size > capacity) {
            capacity = size;
            ptr = (char*)realloc(ptr, capacity);
        }
    }

    void Clear() { MemSet(ptr, 0, capacity); size = 0; }

    // Append
    void Insert(int index, char value)
    {
        ASSERT(index < size || index > 0);
        const int oldLen = Length();
        Reserve(oldLen + 2);
        char* slow = ptr + oldLen;
        char* fast = slow - 1;
        // shift all right characters to 1 char right
        while (fast >= ptr + index)
        {
            *slow-- = *fast--;
        }
        ptr[index] = value;
        ptr[++size] = '\0';
    }

    void Insert(int index, const char* other)
    {
        ASSERT(index < size || index > 0);
        const int otherLen = (int)StrLen(other);
        const int oldLen = Length();
        Reserve(oldLen + otherLen + 2);
        char* curr = ptr + oldLen - 1;
        while (curr >= ptr + index)
        {
            curr[otherLen + 1] = *curr--;
        }
        MemCpy(ptr + index, other, otherLen);
        size += otherLen;
        ptr[++size] = '\0';
    }

    void Insert(int index, const String& value) { Insert(index, value.CStr()); }

    void Append(int value)   { char buff[16]{}; IntToString(buff, value); Append(buff); }
    void Append(float value) { char buff[16]{}; FloatToString(buff, value); Append(buff); }

    static String From(int value)   { char buff[16]{}; IntToString(buff, value); return buff; }
    static String From(float value) { char buff[16]{}; FloatToString(buff, value); return buff; }

    void operator += (char _char) { AppendChar(_char); }
    void operator += (const String& string) { Append(string); }

    void operator += (int  value)  { Append(value); }
    void operator += (float value) { Append(value); }

    String operator + (char _char)           const { String res(*this); res.AppendChar(_char);  return res; }
    String operator + (const char* string)   const { String res(*this); res.Append(string); return res; }
    String operator + (const String& string) const { String res(*this); res.Append(string); return res; }

    String operator + (int value)   const { String res(*this); res.Append(value); return *this; }
    String operator + (float value) const { String res(*this); res.Append(value); return *this; }

    void AppendChar(char _char)
    {
        CapacityCheck(1);
        ptr[size++] = _char;
        ptr[size] = '\0';
    }

    void Append(const char* _char)
    {
        const int len = StrLen(_char);
        CapacityCheck(len);
#pragma warning(suppress : 4996)
        strncat(ptr + size, _char, len + 1ull);
        size += len;
        ptr[size] = '\0';
    }

    void Append(const char* _char, const int len)
    {
        CapacityCheck(len);
#pragma warning(suppress : 4996)
        strncat(ptr + size, _char, size_t(len) + 1ull);
        size += len;
        ptr[size] = '\0'; // for safety
    }

    String AppendCopy(const char* _char)
    {
        const int len = StrLen(_char);
        char* buffer = (char*)calloc(len + size, sizeof(char));
        MemCpy(buffer, ptr, size);
#pragma warning(suppress : 4996)
        strncat(buffer + size, _char, len + 1ull);
        return String(buffer);
    }

    void Append(const String& string) { Append(string.CStr()); }

    // Find
    int FindIndex(char _char) const
    {
        for (int i = 0; i < size; ++i)
            if (ptr[i] == _char) return i;
        return -1;
    }

    int FindIndex(const char* str, const int len) const
    {
        for (int i = 0; i < size - (len - 1); ++i)
            if (!strncmp(str, ptr + i, len)) return i;
        return -1;
    }

    int FindIndex(const char* str) const
    {
        const int len = (int)StrLen(str);
        for (int i = 0; i < size - (len - 1); ++i)
            if (!strncmp(str, ptr + i, len)) return i;
        return -1;
    }

    // returns true if removed
    void Remove(const char* _char)
    {
        const int otherLen = StrLen(_char);
        const int findIndex = FindIndex(_char, otherLen);

        if (findIndex != -1)
        {
            char* findPos = ptr + findIndex;
            char* otherEnd = findPos + otherLen;
            // shift all right characters to the find pos
            while (*otherEnd) {
                *findPos++ = *otherEnd++;
            }
            // set all removed characters to null terminator
            while (*findPos) {
                *findPos++ = '\0';
            }
            return true;
        }
        return false;
    }

    bool StartsWith(const char* other, int len) const
    {
        ASSERT(len > size);
        while (len--)
            if (other[len] != ptr[len]) 
                return false;
        return true;
    }

    void StartsWith(const String& other) const
    {
        int len = other.Length();
        ASSERT(len > size);
        while (len--)
            if (other[len] != ptr[len]) 
                return false;
        return true;
    }

    bool EndsWith(const char* other, int len) const
    {
        ASSERT(len > size);
        const int diff = Abs(len - size);
        while (len--)
            if (other[len] != ptr[len + diff]) 
                return false;
        return true;
    }

    bool EndsWith(const String& other) const
    {
        int len = other.Length();
        ASSERT(len > size);
        const int diff = Abs(len - size);
        while (len--)
            if (other[len] != ptr[len + diff]) 
                return false;
        return true;
    }

    bool Remove(const String& str)         { return Remove(str.CStr()); }

    bool Contains(char _char)        const { return FindIndex(_char) > 0; }
    bool Contains(const char* _char) const { return FindIndex(_char) > 0; }
    bool Contains(const String& str) const { return FindIndex(str.CStr()) > 0; }

    bool Replace(int start, int end, const char* cstr)
    {
        if (end > size) false;
        const int len = end - start;
        MemCpy(ptr + start, cstr, len);
        return true;
    }

    /// <summary> not suitable for big strings (1k-2k char) you can generate your own algorithm for that </summary>
    bool Replace(const char* from, const char* _new)
    {
        const int fromLen = (int)StrLen(from);
        const int toLen = (int)StrLen(_new);

        int fromIndex = FindIndex(from, fromLen);

        if (fromIndex != -1) // this means finded
        {
            const int newLen = Length() + (fromLen - toLen);

            Reserve(newLen);
            size = newLen;
            const int tailLen = StrLen(ptr + fromIndex + fromLen);
            memmove(ptr + fromIndex + toLen, ptr + fromIndex + fromLen, tailLen + 1 * sizeof(char));
            MemCpy(ptr + fromIndex, _new, toLen * sizeof(char));
            return true;
        }
        return false;
    }

    /// <summary> not suitable for big strings (1k-2k char) you can generate your own algorithm for that </summary>
    /// <returns> number of instance removed </returns>
    int ReplaceAll(const char* old, const char* _new)
    {
        const int _newLen = StrLen(_new);
        bool strResult = true;
        int removedCount = 0;
        while (strResult == true) {
            strResult &= Replace(old, _new);
            ++removedCount;
        }
        return removedCount;
    }

    bool Replace(const String& old, const String& _new)
    {
        return Replace(old.CStr(), _new.CStr());
    }

    String SubString(int begin, int end) const
    {
        if (begin > size) return String();
        const int buffSize = end - begin;
        char* buffer = (char*)calloc(buffSize + 1ull, 1ull);
        MemCpy(buffer, ptr, buffSize);
        return String(buffer);
    }

    int Capacity() const { return capacity; }
    int Length()   const { return size; }
    
private:
    FINLINE bool CapacityCheck(int len)
    {
        if (size + len + 1 >= capacity)
        {
            int newLen = size + len;
            capacity = newLen + (newLen / 2);
            ptr = (char*)realloc(ptr, capacity + 1);
            return true;
        }
        return false;
    }
};

template<> struct Hasher<String<>>
{
    static FINLINE uint64 Hash(const String<>& x)
    {
        return WYHash::Hash(x.CStr(), x.Length());
    }
};
