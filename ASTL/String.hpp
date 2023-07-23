
#pragma once

#include "Random.hpp"
#include "Algorithms.hpp"
#include "Memory.hpp"

// https://github.com/WojciechMula/simd-string/blob/master/strcmp.cpp

#ifdef AX_SUPPORT_SSE
inline int StringLength(const char* s)
{
    __m128i* mem = reinterpret_cast<__m128i*>(const_cast<char*>(s));
    const __m128i zeros = _mm_setzero_si128();
    
    for (int result = 0; /**/; mem++, result += 16) 
    {
        const __m128i data = _mm_loadu_si128(mem);
        const uint8_t mode = _SIDD_UBYTE_OPS | _SIDD_CMP_EQUAL_EACH | _SIDD_LEAST_SIGNIFICANT;

        if (_mm_cmpistrc(data, zeros, mode)) 
        {
            int idx = _mm_cmpistri(data, zeros, mode);
            return result + idx;
        }
    }
}

inline int StringCompare(const char* s1, const char* s2)
{
    if (s1 == s2) 
        return 0;
    
    const __m128i* ptr1 = (const __m128i*)s1;
    const __m128i* ptr2 = (const __m128i*)s2;
    
    for (/**/; /**/; ptr1++, ptr2++) {
    
        const __m128i a = _mm_loadu_si128(ptr1);
        const __m128i b = _mm_loadu_si128(ptr2);
    
        const uint8_t mode = _SIDD_UBYTE_OPS | _SIDD_CMP_EQUAL_EACH | _SIDD_NEGATIVE_POLARITY | _SIDD_LEAST_SIGNIFICANT;
      
        if (_mm_cmpistrc(a, b, mode)) {
            // a & b are different (not counting past-zero bytes)
            int idx = _mm_cmpistri(a, b, mode);
    
            const uint8_t b1 = ((char*)ptr1)[idx];
            const uint8_t b2 = ((char*)ptr2)[idx];
    
            if      (b1 < b2) return -1;
            else if (b1 > b2) return +1;
            else  return 0;

        } else if (_mm_cmpistrz(a, b, mode)) {
            // a & b are same, but b contains a zero byte
            break;
        }
    }

    return false;
}
#else
inline int StrLen(const char* s)
{
    const char* begin = s;
    while (*s) s++;
    return s - begin;
}

inline int StrCmp(const char* a, const char* b)
{
    for (; *a && *b; a++, b++)
    {
        if (*a != *b)
        {
            if (*a < *b) return -1;
            else         return +1; // greater
        }
    }
    return *a == *b;// strings are equal
}
#endif

// small string optimization
struct StringAllocator
{
    char stack[16]{};

    char* Allocate(int size)
    {
        if (size < 16)
        {
            return stack;
        }
        return new char[size]{};
    }
    
    void Deallocate(char* ptr, int size)
    {
        if (ptr == stack)
            MemSet<1, 16>(stack, 0, 16);
        else
            delete[] ptr;
    }

    char* Reallocate(char* ptr, int oldCount, int count)
    {
        if (ptr == stack && count < 16)
            return stack;
        
        char* newPtr = new char[count]{};
        MemCpy<1>(newPtr, ptr, oldCount);
        Deallocate(ptr, oldCount);
        return newPtr;
    }
};

class String
{
public:
    char* ptr    = nullptr;
    int capacity = 0;
    int size     = 0;
    StringAllocator allocator{};

public:
    String() : ptr(nullptr), capacity(0), size(0)
    { }

    ~String() { Clear(); }
    
    explicit String(int _capacity) : capacity(_capacity), size(0)
    {
        ptr = allocator.Allocate(_capacity + 1);
    }
     
    String(char* cstr)
    {
        int clen = StringLength(cstr);
        ptr = allocator.stack;
        if (clen >= 16) ptr = allocator.Allocate(clen + 1);
        Copy(ptr, cstr, clen);
        capacity = size = clen;
    }

    String(const char* cstr)
    {
        int clen = StringLength(cstr);
        ptr = allocator.stack;
        if (clen >= 16)  ptr = allocator.Allocate(clen + 1);
        Copy(ptr, cstr, clen);
        capacity = size = clen;
    }

    String(const char* begin, int count)
    {
        ptr = allocator.stack;
        if (count >= 16) ptr = allocator.Allocate(count + 1);
        Copy(ptr, begin, count);
        capacity = size = count;
    }

    // copy constructor
    String(const String& other)
    {
        Asign(other);
    }

    // move constructor
    String(String&& other) noexcept
    {
        capacity = other.capacity;
        size = other.size;
        ptr  = other.ptr == other.allocator.stack ? allocator.stack : other.ptr;
        SmallMemCpy(allocator.stack, other.allocator.stack, 16);
        MemSet<1, 16>(other.allocator.stack, 0, 16);
        other.ptr = nullptr; 
        other.size = other.capacity = 0;
    }

    const char* begin() const { return ptr; }
    const char* end()   const { return ptr + size; }
    char*       begin() { return ptr; }
    char*       end()   { return ptr + size; }

    int  Length() const { return size; }
    bool Empty()  const { return size == 0; }
    
          char* CStr()        { return ptr ? ptr : allocator.stack; }
    const char* CStr()  const { return ptr ? ptr : allocator.stack; }

    const char& operator[] (int index) const { ASSERT(index < size && index > 0); return ptr[index]; }
          char& operator[] (int index)       { ASSERT(index < size && index > 0); return ptr[index]; }

    String& operator = (const String& right) 
    {
        Asign(right);
        return *this;
    }

    bool operator == (const String& other) { return StringCompare(ptr, other.ptr) == 0; }
    bool operator != (const String& other) { return StringCompare(ptr, other.ptr) != 0; }
    bool operator >  (const String& other) { return StringCompare(ptr, other.ptr) == 1; }
    bool operator <  (const String& other) { return StringCompare(ptr, other.ptr) == -1; }
    bool operator >= (const String& other) { return StringCompare(ptr, other.ptr) >= 0; }
    bool operator <= (const String& other) { return StringCompare(ptr, other.ptr) <= 0; }

    String operator + (const String& other) { String cpy(other.size + size + 1);          cpy.Append(ptr); cpy.Append(other.ptr); return (String&&)cpy; }
    String operator + (const char*        other) { String cpy(StringLength(other) + size + 1); cpy.Append(ptr); cpy.Append(other); return (String&&)cpy; }
    
    String& operator += (const String& other) { Append(other.ptr); return *this; }
    String& operator += (const char*   other) { Append(ptr);       return *this; }

    String operator += (float f) { char arr[16]{}; FloatToString(arr, f); Append(arr); return *this; }
    String operator += (int   i) { char arr[16]{}; IntToString(arr, i);   Append(arr); return *this; }
    String operator +  (float f) { char arr[16]{}; FloatToString(arr, f); String cpy = *this; cpy.Append(arr); return (String&&)cpy; }
    String operator +  (int   i) { char arr[16]{}; IntToString(arr, i);   String cpy = *this; cpy.Append(arr); return (String&&)cpy; }

    void Asign(const String& other)
    {
        if (&other == this) return;
        if (other.size >= size) GrowIfNecessarry(other.size - size);

        if (other.size >= 16)
            Copy(ptr, other.ptr, other.size);
        else
        {
            allocator.Deallocate(ptr, size);
            MemCpy<1, 16>(allocator.stack, other.allocator.stack, 16);
            ptr = allocator.stack;
        }

        size = other.size;
    }

    void Append(float f) { char arr[16]{}; FloatToString(arr, f); Append(arr); }
    void Append(int i)   { char arr[16]{}; IntToString(arr, i);   Append(arr); }
    
    static bool CompareN(const char* a, const char* b, int n) 
    {
        for (int i = 0; i < n; i++)
            if (a[i] != b[i]) return false;
        return true;
    }

    int IndexOf(const char* other) const
    {
        int otherLen = StringLength(other);
        for (int i = size - otherLen; i >= 0; i--) {
            if (CompareN(ptr + i, other, otherLen))
                return i;
        }
        return -1;
    }

    int IndexOf(const String& other) const 
    {
        return IndexOf(other.CStr()) != -1;
    }

    bool Contains(const char* other) const
    {
        return IndexOf(other) != -1;
    }
    
    bool Contains(const String& other) const
    {
        return IndexOf(other.CStr()) != -1;
    }

    void Reserve(int size)
    {
        if (!ptr)
        {
            ptr = allocator.Allocate(size);
            capacity = size;
        }
        else if (capacity < size)
        {
            ptr = allocator.Reallocate(ptr, this->size, size);
            capacity = size;
        }
    }

    void Clear()
    {
        if (ptr)
        {
            allocator.Deallocate(ptr, size);
            ptr = nullptr;
            size = capacity = 0;
        }
    }

    void Insert(int index, char c)
    {
        OpenSpace(index, 1);
        ptr[size++] = c;
    }

    void Insert(int index, const char* ptr)
    {
        int size = StringLength(ptr);
        OpenSpace(index, size);
        Copy(this->ptr + index, ptr, size);
    }

    void Insert(int index, const String& other) 
    {
        Insert(index, other.ptr); 
    }
    
    String SubString(int index)
    {
        int count = size - index;
        ASSERT(index + count <= size);
        String ret(count + 2);
        ret.Append(ptr + index, size - count);
        return ret;
    }

    String SubString(int index, int count)
    {
        ASSERT(index + count <= size);
        String ret(count + 2);
        ret.Append(ptr + index, count);
        return ret;
    }

    void Remove(int index, int count = 1)
    {
        RemoveSpace(index, count);
        while (count--)
            ptr[--size] = 0;
    }

    void Remove(const char* find)
    {
        int x = IndexOf(find);
        if (x == -1) return;
        Remove(x, StringLength(find));
    }

    void Replace(const char* find, const char* replace)
    {
        int x = IndexOf(find);
        if (x == -1) return;
        Insert(x, replace);
        int repLen = StringLength(replace);
        Remove(x + repLen, StringLength(find));
    }

    void Replace(const String& find, const String& replace)
    {
        int x = IndexOf(find.c_str());
        if (x == -1) return;
        Insert(x, replace);
        Remove(x + replace.Length(), find.Length());
    }

    void Append(char c)
    {
        GrowIfNecessarry(1);
        ptr[size++] = c;
    }

    void Append(const char* other, int count)
    {
        GrowIfNecessarry(count);
        Copy(ptr + size, other, count);
        size += count;
    }

    void Append(const char* other)
    {
        Append(other, StringLength(other));
    }

    void Append(const String& other) { Append(other.ptr); }

    void RemoveSpace(int _index, int _count)
    {
        ASSERT((_index + _count) <= this->size);
        int newSize = Max(this->size - _count, 0);

        int i = _index;
        int j = _index + _count;

        // *******i****j***   < opens space with incrementing both of the pointers (two pointer algorithm)
        while (j < this->size)
        {
            ptr[i++] = ptr[j++];
        }
    }

    void OpenSpace(int _index, int _count)
    {
        GrowIfNecessarry(_count);

        long i = Min(size + _count, capacity);
        int  j = size;
        ASSERT(i <= INT32_MAX);

        while (j >= _index)
        {
            ptr[--i] = ptr[--j];
        }
        size += _count;
    }

    void GrowIfNecessarry(int _size)
    {
        int newSize = size + _size + 1;
        
        if (newSize < 16)
            return;

        if (newSize >= capacity)
        {
            constexpr int InitialSize = 64;
            newSize = Max(CalculateArrayGrowth(size + _size), InitialSize);
            if (ptr) 
                ptr = allocator.Reallocate(ptr, capacity, newSize);
            else
                ptr = allocator.Allocate(newSize);
            
            for (int i = size; i < newSize; i++)
                ptr[i] = 0;

            capacity    = newSize; 
        }   
    }

    char* c_str()             { return ptr ? ptr : allocator.stack; }
    const char* c_str() const { return ptr ? ptr : allocator.stack; }
#ifdef ASTL_STL_COMPATIBLE
    bool empty()       const  { return size == 0; }
#endif
};

template<> struct Hasher<String>
{
    static FINLINE uint64 Hash(const String& x)
    {
        return WYHash::Hash(x.CStr(), x.Length());
        // return MurmurHash64(x.CStr(), x.Length(), 0xa0761d6478bd642full);
    }
};