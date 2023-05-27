#pragma once
#include "Common.hpp"
// #include <iosfwd> // for overriding << operator for std::cout
// #include <sstream>
#include <stdlib.h>
#include "Algorithms.hpp"

inline wchar_t* ToWCharArray(const char* string)
{
    const size_t len = StrLen(string);
    wchar_t* buffer = (wchar_t*)malloc(len + 1ull * sizeof(wchar_t));
#pragma warning(suppress : 4996)
    mbstowcs(buffer, string, len + 1ull);
    return buffer;
}

inline char* ToCharArray(const wchar_t* string)
{
    const size_t len = wcslen(string);
    char* buffer = (char*)malloc(len + 1ull);
#pragma warning(suppress : 4996)
    wcstombs(buffer, string, len + 1ull);
    return buffer;
}

enum class StrResult : int
{
    None = 0, NotFinded = 0, False = 0, // negatives
    Finded = 1, Success = 1, True = 1,  // positives
    IndexOutOfArray = 2                 // errors
};

class String
{
public:
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

    String(int _size) : size(0), capacity(_size + 1)
    {
        ptr = (char*)calloc(capacity, sizeof(char));
    }
    // copy constructor
    String(const String& other) : size(other.size), capacity(other.capacity)
    {
        ptr = (char*)calloc(capacity, sizeof(char));
        memcpy(ptr, other.ptr, size);
    }
    // move constructor 
    String(String&& other) noexcept : size(other.size), capacity(other.capacity), ptr(other.ptr)
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
        memcpy(ptr, _ptr, size);
    }

    String(char* _ptr) : size(StrLen(_ptr))
    {
        capacity = size + 32;
        ptr = (char*)calloc(capacity, 1);
        memcpy(ptr, _ptr, size);
    }
    // asign operator
    String& operator = (const String& right)
    {
        CapacityCheck(right.size);
        memcpy(ptr, right.ptr, right.size);
        size = right.size;
        return *this;
    }

    bool operator == (const String& b) const { return  Compare(*this, b); }
    bool operator != (const String& b) const { return !Compare(*this, b); }
    bool operator == (const char* b)   const { return !strcmp(ptr, b); }
    bool operator != (const char* b)   const { return  strcmp(ptr, b); }
    char operator[](int index)         const { ax_assert(index < size && index > 0); return ptr[index]; }
    char& operator[](int index)              { ax_assert(index < size && index > 0); return ptr[index]; }

    char* begin() { return ptr; }
    char* end() { return ptr + size; }
    const char* cbegin() const { return ptr; }
    const char* cend()   const { return ptr + size; }

    const char* CStr() const { return ptr; };
    const char* c_str() const { return ptr; };
    const char* data() const { return ptr; }
    char* data() { return ptr; }

    static FINLINE bool Compare(const String& a, const String& b)
    {
        return a.size == b.size && !strcmp(a.ptr, b.ptr);
    }

    inline void Reset()
    {
        memset(ptr, 0, size);
        size = 0;
    }

    inline void Set(const char* str)
    {
        Clear();
        CapacityCheck(StrLen(str));
        memcpy(ptr, str, StrLen(str));
    }

    inline void Set(const String& str)
    {
        Clear();
        CapacityCheck(str.size);
        size = str.Length();
        memcpy(ptr, str.CStr(), str.size);
    }

    inline void Reserve(int size)
    {
        if (size > capacity) {
            capacity = size;
            ptr = (char*)realloc(ptr, capacity);
        }
    }

    inline void Clear() { memset(ptr, 0, capacity); size = 0; }

    // Append
    StrResult Insert(int index, char value)
    {
        if (index > size || index < 0) return StrResult::IndexOutOfArray;
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
        return StrResult::Success;
    }

    StrResult Insert(int index, const char* other)
    {
        if (index > size || index < 0) return StrResult::IndexOutOfArray;
        const int otherLen = (int)StrLen(other);
        const int oldLen = Length();
        Reserve(oldLen + otherLen + 2);
        char* curr = ptr + oldLen - 1;
        while (curr >= ptr + index)
        {
            curr[otherLen + 1] = *curr--;
        }
        memcpy(ptr + index, other, otherLen);
        size += otherLen;
        ptr[++size] = '\0';
        return StrResult::Success;
    }

    StrResult Insert(int index, const String& value)
    {
        return Insert(index, value.CStr());
    }

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

    String operator + (int value) const { String res(*this); res.Append(value); return *this; }
    String operator + (float value) const { String res(*this); res.Append(value); return *this; }

    void AppendChar(char _char)
    {
        CapacityCheck(1);
        ptr[size++] = _char;
        ptr[size] = '\0';
    }

    inline void Append(const char* _char)
    {
        const int len = StrLen(_char);
        CapacityCheck(len);
#pragma warning(suppress : 4996)
        strncat(ptr + size, _char, len + 1ull);
        size += len;
        ptr[size] = '\0';
    }

    inline void Append(const char* _char, const int len)
    {
        CapacityCheck(len);
#pragma warning(suppress : 4996)
        strncat(ptr + size, _char, size_t(len) + 1ull);
        size += len;
        ptr[size] = '\0'; // for safety
    }

    inline String AppendCopy(const char* _char)
    {
        const int len = StrLen(_char);
        char* buffer = (char*)calloc(len + size, sizeof(char));
        memcpy(buffer, ptr, size);
#pragma warning(suppress : 4996)
        strncat(buffer + size, _char, len + 1ull);
        return String(buffer);
    }

    inline void Append(const String& string) { Append(string.CStr()); }

    // Find
    int FindIndex(char _char) const
    {
        for (int i = 0; i < size; ++i)
            if (ptr[i] == _char) return i;
        return -1;
    }

    inline int FindIndex(const char* str, const int len) const
    {
        for (int i = 0; i < size - (len - 1); ++i)
            if (!strncmp(str, ptr + i, len)) return i;
        return -1;
    }

    inline int FindIndex(const char* str) const
    {
        const int len = (int)StrLen(str);
        for (int i = 0; i < size - (len - 1); ++i)
            if (!strncmp(str, ptr + i, len)) return i;
        return -1;
    }

    // Remove
    inline StrResult Remove(char _char)
    {
        const int index = FindIndex(_char);
        if (index == -1) return StrResult::NotFinded;
        memmove(ptr + index, ptr + index + 1ull, 1ull);
        ptr[--size] = '\0';
        return StrResult::Success;
    }

    inline StrResult Remove(const char* _char)
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

            return StrResult::Success;
        }
        return StrResult::NotFinded;
    }

    inline StrResult StartsWith(const char* other, int len) const
    {
        if (size < len) StrResult::IndexOutOfArray;

        while (len--)
        {
            if (other[len] != ptr[len]) return StrResult::False;
        }
        return StrResult::True;
    }

    StrResult StartsWith(const String& other) const
    {
        int len = other.Length();
        if (size < len) StrResult::IndexOutOfArray;
        while (len--)
        {
            if (other[len] != ptr[len]) return StrResult::False;
        }
        return StrResult::True;
    }

    inline StrResult EndsWith(const char* other, int len) const
    {
        if (size < len) StrResult::IndexOutOfArray;
        const int diff = Abs(len - size);
        while (len--)
        {
            if (other[len] != ptr[len + diff]) return StrResult::False;
        }
        return StrResult::True;
    }

    StrResult EndsWith(const String& other) const
    {
        int len = other.Length();
        if (size < len) StrResult::IndexOutOfArray;
        const int diff = Abs(len - size);
        while (len--)
        {
            if (other[len] != ptr[len + diff]) return StrResult::False;
        }
        return StrResult::True;
    }

    StrResult Remove(const String& str)
    {
        return Remove(str.CStr());
    }

    StrResult Find(char _char) const { return FindIndex(_char) ? StrResult::Success : StrResult::NotFinded; }

    StrResult Find(const char* _char) const { return FindIndex(_char) ? StrResult::Success : StrResult::NotFinded; }

    StrResult Find(const String& str) const { return FindIndex(str.CStr()) ? StrResult::Success : StrResult::NotFinded; }

    StrResult Replace(int start, int end, const char* cstr)
    {
        if (end > capacity) return StrResult::IndexOutOfArray;
        const int len = end - start;
        memcpy(ptr + start, cstr, len);
        return StrResult::Success;
    }

    /// <summary> not suitable for big strings (1k-2k char) you can generate your own algorithm for that </summary>
    StrResult Replace(const char* from, const char* _new)
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
            memcpy(ptr + fromIndex, _new, toLen * sizeof(char));
            return StrResult::Success;
        }
        return StrResult::NotFinded;
    }

    /// <summary> not suitable for big strings (1k-2k char) you can generate your own algorithm for that </summary>
    /// <returns> number of instance removed </returns>
    int ReplaceAll(const char* old, const char* _new)
    {
        StrResult strResult = StrResult::Success;
        const int _newLen = StrLen(_new);
        int removedCount = 0;
        while (strResult == StrResult::Success) {
            strResult = Replace(old, _new);
            ++removedCount;
        }
        return removedCount;
    }

    StrResult Replace(const String& old, const String& _new)
    {
        return Replace(old.CStr(), _new.CStr());
    }

    String SubString(int begin, int end) const
    {
        if (begin > size) return String();
        const int buffSize = end - begin;
        char* buffer = (char*)calloc(buffSize + 1ull, 1ull);
        memcpy(buffer, ptr, buffSize);
        return String(buffer);
    }
    /*
    friend std::ostream& operator<<(std::ostream& cout, String& wstr) {
        return cout << wstr.CStr();
    }
    */
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
private:
    int capacity = 0;
    int size = 0;
public:
    char* ptr = nullptr;
};

class WString
{
public:
    ~WString()
    {
        if (ptr) { free(ptr); ptr = nullptr;  size = 0; capacity = 0; }
    }
    // copy constructor
    WString(const WString& other) : size(other.size), capacity(other.capacity)
    {
        ptr = (wchar_t*)calloc(capacity, sizeof(wchar_t));
        wmemcpy(ptr, other.ptr, size);
    }
    // move constructor 
    WString(WString&& other) noexcept : size(other.size), capacity(other.capacity), ptr(other.ptr)
    {
        other.ptr = nullptr;
    }

    WString() : size(0), capacity(32)
    {
        ptr = (wchar_t*)calloc(capacity, sizeof(wchar_t));
    }

    WString(int _size) : size(0), capacity(_size)
    {
        ptr = (wchar_t*)calloc(capacity, sizeof(wchar_t));
    }

    WString(const wchar_t* _ptr) : size(WStrLen(_ptr))
    {
        capacity = size + 32;
        ptr = (wchar_t*)calloc(capacity, sizeof(wchar_t));
        wmemcpy(ptr, _ptr, size);
    }

    WString(wchar_t* _ptr) : size(WStrLen(_ptr))
    {
        capacity = size + 32;
        ptr = (wchar_t*)calloc(capacity, sizeof(wchar_t));
        memcpy(ptr, _ptr, size * sizeof(wchar_t));
    }

    WString& operator = (const WString& right)
    {
        CapacityCheck(right.size);
        wmemcpy(ptr, right.ptr, right.size);
        size = right.size;
        return *this;
    }

    bool operator == (const WString& b) { return  Compare(*this, b); }
    bool operator != (const WString& b) { return !Compare(*this, b); }
    wchar_t operator[](int index) const { ax_assert(index < size && index > 0); return ptr[index]; }
    wchar_t& operator[](int index)      { ax_assert(index < size && index > 0); return ptr[index]; }

    wchar_t* begin() { return ptr; }
    wchar_t* end() { return ptr + size; }
    const wchar_t* cbegin() const { return ptr; }
    const wchar_t* cend()   const { return ptr + size; }

    const wchar_t* CStr() const { return ptr; };

    static FINLINE bool Compare(WString a, WString b)
    {
        return a.size == b.size && !wcscmp(a.ptr, b.ptr);
    }

    inline void Reset()
    {
        memset(ptr, 0, size);
        size = 0;
    }

    inline void Set(const wchar_t* str)
    {
        Clear();
        CapacityCheck(WStrLen(str));
        memcpy(ptr, str, WStrLen(str) * sizeof(wchar_t));
    }

    inline void Set(const WString& str)
    {
        Clear();
        CapacityCheck(str.Length());
        size = str.Length();
        memcpy(ptr, str.CStr(), str.Length() * sizeof(wchar_t));
    }

    inline void Reserve(int size)
    {
        if (size > capacity) {
            capacity = size;
            ptr = (wchar_t*)realloc(ptr, capacity * sizeof(wchar_t));
        }
    }

    inline void Clear() { memset(ptr, 0, size * sizeof(wchar_t)); size = 0; }

    // Append
    inline StrResult Insert(int index, wchar_t value)
    {
        if (index > size || index < 0) return StrResult::IndexOutOfArray;
        const int oldLen = Length();
        Reserve(oldLen + 2);
        wchar_t* slow = ptr + oldLen;
        wchar_t* fast = slow - 1;
        // shift all right characters to 1 char right
        while (fast >= ptr + index)
        {
            *slow-- = *fast--;
        }
        ptr[index] = value;
        ptr[++size] = L'\0';
        return StrResult::Success;
    }

    inline StrResult Insert(int index, const wchar_t* other)
    {
        if (index > size || index < 0) return StrResult::IndexOutOfArray;
        const int otherLen = (int)wcslen(other);
        const int oldLen = Length();
        Reserve(oldLen + otherLen + 2);
        wchar_t* curr = ptr + oldLen - 1;
        while (curr >= ptr + index)
        {
            curr[otherLen + 1] = *curr--;
        }
        memcpy(ptr + index, other, otherLen * sizeof(wchar_t));
        size += otherLen;
        ptr[++size] = L'\0';
        return StrResult::Success;
    }

    inline StrResult Insert(int index, const WString& value)
    {
        return Insert(index, value.CStr());
    }

    void Append(int  value) { wchar_t buff[16]; swprintf(buff, L"%d", value);   Append(buff); }
    void Append(float value) { wchar_t buff[16]; swprintf(buff, L"%f", value);   Append(buff); }

    static WString From(int  value) { wchar_t buff[16]; swprintf(buff, L"%d", value); return buff; }
    static WString From(float value) { wchar_t buff[16]; swprintf(buff, L"%f", value); return buff; }

    void operator += (char _char) { AppendChar(_char); }
    void operator += (const WString& string) { Append(string); }

    void operator += (int  value) { Append(value); }
    void operator += (float value) { Append(value); }

    WString operator + (wchar_t _char)         const { WString res = *this; res.Append(_char);  return res; }
    WString operator + (const WString& string) const { WString res = *this; res.Append(string); return res; }

    WString operator + (int  value) const { WString res = *this; res.Append(value); return *this; }
    WString operator + (float value) const { WString res = *this; res.Append(value); return *this; }

    void AppendChar(wchar_t _char)
    {
        CapacityCheck(1);
        ptr[size++] = _char;
        ptr[size] = '\0';
    }

    WString Append(const wchar_t* _char)
    {
        const int len = WStrLen(_char);
        CapacityCheck(len);
#pragma warning(suppress : 4996)
        wcsncat(ptr, _char, size_t(len) + 1ull);
        size += len;
        ptr[size] = '\0';
        return *this;
    }

    WString Append(const wchar_t* _char, const int len)
    {
        CapacityCheck((int)len);
#pragma warning(suppress : 4996)
        wcsncat(ptr, _char, len + 1ull);
        size += len;
        ptr[size] = '\0';
        return *this;
    }

    WString AppendCopy(const wchar_t* _char) const
    {
        const int len = WStrLen(_char);
        wchar_t* buffer = (wchar_t*)calloc(len + size, sizeof(wchar_t));
        wmemcpy(buffer, ptr, (size_t)size);
#pragma warning(suppress : 4996)
        wcsncat(buffer, _char, size_t(len) + 1ull);
        return WString(buffer);
    }

    void Append(const WString& string) { Append(string.CStr()); }

    // Find
    inline int FindIndex(wchar_t _char) const
    {
        for (int i = 0; i < size; ++i)
            if (ptr[i] == _char) return i;
        return false;
    }

    inline int FindIndex(const wchar_t* str) const
    {
        const int len = WStrLen(str);
        for (int i = 0; i < size - (len - 1); ++i)
            if (!wcsncmp(str, ptr + i, len)) return i;
        return -1;
    }

    inline int FindIndex(const wchar_t* str, const size_t wlen) const
    {
        const int len = (int)wlen;
        for (int i = 0; i < size - (len - 1); ++i)
            if (!wcsncmp(str, ptr + i, len)) return i;
        return -1;
    }
    // Remove
    inline StrResult Remove(wchar_t _char)
    {
        const int index = FindIndex(_char);
        if (index == -1) return StrResult::NotFinded;
        memmove(ptr + index, ptr + index + 1, sizeof(wchar_t));
        ptr[--size] = L'\0';
        return StrResult::Success;
    }

    inline StrResult Remove(const wchar_t* _char)
    {
        const int otherLen = WStrLen(_char);
        const int findIndex = FindIndex(_char, otherLen);

        if (findIndex != -1)
        {
            wchar_t* findPos = ptr + findIndex;
            wchar_t* otherEnd = findPos + otherLen;
            // shift all right characters to the find pos
            while (*otherEnd) {
                *findPos++ = *otherEnd++;
            }
            // set all removed characters to null terminator
            while (*findPos) {
                *findPos++ = L'\0';
            }

            return StrResult::Success;
        }
        return StrResult::NotFinded;
    }

    inline StrResult StartsWith(const wchar_t* other, int len) const
    {
        if (size < len) StrResult::IndexOutOfArray;

        while (len--)
        {
            if (other[len] != ptr[len]) return StrResult::False;
        }
        return StrResult::True;
    }

    inline StrResult StartsWith(const WString& other) const
    {
        int len = other.Length();
        if (size < len) StrResult::IndexOutOfArray;

        while (len--)
        {
            if (other[len] != ptr[len]) return StrResult::False;
        }
        return StrResult::True;
    }

    StrResult Remove(const WString& str)
    {
        return Remove(str.CStr());
    }

    StrResult Find(wchar_t _char) const { return FindIndex(_char) ? StrResult::Success : StrResult::NotFinded; }

    StrResult Find(const wchar_t* _char) const { return FindIndex(_char) ? StrResult::Success : StrResult::NotFinded; }

    StrResult Find(const WString& str) const { return FindIndex(str.CStr()) ? StrResult::Success : StrResult::NotFinded; }

    StrResult Replace(int start, int end, const wchar_t* cstr)
    {
        if (end > capacity) return StrResult::IndexOutOfArray;
        const int len = end - start;
        wmemcpy(ptr + start, cstr, len);
        return StrResult::Success;
    }

    StrResult Replace(const wchar_t* from, const wchar_t* _new)
    {
        const int fromLen = (int)wcslen(from);
        const int toLen = (int)wcslen(_new);

        int fromIndex = FindIndex(from, fromLen);

        if (fromIndex != -1) // this means finded
        {
            const int newLen = Length() + (fromLen - toLen);
            Reserve(newLen);
            size = newLen;

            const size_t tailLen = wcslen(ptr + fromIndex + fromLen);
            memmove(ptr + fromIndex + toLen, ptr + fromIndex + fromLen, tailLen + 1ull);
            memcpy(ptr + fromIndex, _new, toLen);
            return StrResult::Success;
        }
        return StrResult::NotFinded;
    }
    /// <summary> not suitable for big strings (1k-2k char) you can generate your own algorithm for that </summary>
    /// <returns> removed instance count </returns>
    int ReplaceAll(const wchar_t* old, const wchar_t* _new)
    {
        const int _newSize = WStrLen(_new);
        int removedCount = 0;
        StrResult strResult = StrResult::Success;

        while (strResult == StrResult::Success) {
            strResult = Replace(old, _new);
            ++removedCount;
        }
        return removedCount;
    }

    StrResult Replace(const WString& old, const WString& _new)
    {
        return Replace(old.CStr(), _new.CStr());
    }

    String ToString() const
    {
        char* characters = (char*)malloc(size + 1ull);
#pragma warning(suppress : 4996)
        wcstombs(characters, ptr, size + 1ull);
        return String(characters);
    }

    WString SubString(int begin, int end) const
    {
        if (begin > size) return WString();
        const int buffLen = end - begin;
        wchar_t* buffer = (wchar_t*)calloc(buffLen + 1ull, sizeof(wchar_t));
        wmemcpy(buffer, ptr, (size_t)buffLen);
        return WString(buffer);
    }
    /*
    friend std::ostream& operator << (std::ostream& cout, const WString& wstr) {
        return cout << ToCharArray(wstr.CStr());
    }
    */
    int Capacity() const { return capacity; }
    int Length()   const { return size; }

private:
    FINLINE void CapacityCheck(int len)
    {
        if (size + len + 1 >= capacity)
        {
            int newLen = size + len;
            capacity = newLen + (newLen / 2);
            ptr = (wchar_t*)realloc(ptr, (size_t)capacity * sizeof(wchar_t));
        }
    }

private:
    int capacity;
    int size;
public:
    wchar_t* ptr = nullptr;
};

inline WString ToWString(const String& string)
{
    wchar_t* buffer = (wchar_t*)malloc(string.Length() + 1ull * sizeof(wchar_t));
#pragma warning(suppress : 4996)
    mbstowcs(buffer, string.CStr(), string.Length() + 1ull);
    return WString(buffer);
}

inline String ToString(const WString& string)
{
    char* buffer = (char*)malloc(string.Length() + 1ull * sizeof(char));
#pragma warning(suppress : 4996)
    wcstombs(buffer, string.CStr(), string.Length() + 1ull);
    return String(buffer);
}

template<> struct Hasher<String>
{
    static FINLINE uint64 Hash(const String& x)
    {
        return WYHash::Hash(x.CStr(), x.Length());
    }
};

template<> struct Hasher<WString>
{
    static FINLINE uint64 Hash(const WString& x)
    {
        return WYHash::Hash(x.CStr(), x.Length());
    }
};