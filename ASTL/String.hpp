
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

#include "Profiler.hpp"

// small string optimization
class String
{
public:
    
    union
    {
        char stack[24]{};
        char* ptr[3];
        uint64_t lowHigh[3]; // when using heap high part is =  capacity | (size << 31);
    };

    static constexpr uint32_t LastBit32 = 1u << 31u;

    uint32_t GetCapacity() const
    {
        return uint32_t(lowHigh[1] >> 32ull) & ~LastBit32;
    }
    
    void SetCapacity(int cap) 
    {
        if (!UsingHeap()) return;
        lowHigh[1] &= 0xFFFFFFFFull | (1ull << 63ull);
        lowHigh[1] |= uint64_t(cap) << 32ull;
    }

    uint32_t GetSize() const 
    {
        return !UsingHeap() ? StringLength(stack) : (uint32_t)(lowHigh[1] & 0xFFFFFFFFull);
    }
    
    void SetSize(uint32_t size) 
    {
        if (!UsingHeap()) return;
        lowHigh[1] &= 0xFFFFFFFFull << 32ull; 
        lowHigh[1] |= size;
    }

    bool UsingHeap() const { return lowHigh[1] &  (1ull << 63ull); }
    void SetUsingHeap()    {        lowHigh[1] |= (1ull << 63ull); }

    const char* GetPtr() const { return UsingHeap() ? ptr[0] : stack; }
          char* GetPtr()       { return UsingHeap() ? ptr[0] : stack; }

    void Allocate(int size)
    {
        if (size < 23)
            return;

        int oldSize = GetSize();
        ptr[0] = new char[size] {};
        ptr[2] = nullptr;
        lowHigh[1] = oldSize;

        SetUsingHeap();
    }
    
    void Deallocate()
    {
        if (UsingHeap())
            delete[] GetPtr();

        ptr[0] = ptr[1] = ptr[2] = nullptr;
    }

    void Reallocate(int oldCount, int count)
    {
        if (UsingHeap())
        {
            char* newPtr = new char[count] {};
            char* ptr = this->ptr[0];
            MemCpy<1>(newPtr, ptr, oldCount);
            delete[] ptr;
            this->ptr[0] = newPtr;
        }
        else if (count > 23) // was stack allocating but bigger memory requested
        {
            char* newPtr = new char[count] {};
            MemCpy<1>(newPtr, stack, oldCount);
            ptr[0] = newPtr;
            ptr[2] = nullptr;
            lowHigh[1] = oldCount;
            SetUsingHeap();
        }
    }

public:
    String() 
    { }

    ~String() { Clear(); }
    
    explicit String(int _capacity) 
    {
        Allocate(_capacity + 1);
    }
     
    explicit String(char* cstr)
    {
        int clen = StringLength(cstr);
        Allocate(clen + 1);
        SetCapacity(clen); 
        SetSize(clen);
        Copy(GetPtr(), cstr, clen);
    }

    String(const char* cstr)
    {
        int clen = StringLength(cstr);
        Allocate(clen + 1);
        SetCapacity(clen);
        SetSize(clen);
        SmallMemCpy(GetPtr(), cstr, clen);
    }

    String(const char* begin, int count)
    {
        Allocate(count + 1);
        SetCapacity(count); 
        SetSize(count);
        Copy(GetPtr(), begin, count);
    }

    // copy constructor
    String(const String& other)
    {
        Asign(other);
    }

    // move constructor
    String(String&& other) noexcept
    {
        SmallMemCpy(stack, other.stack, 24);
        other.ptr[0] = other.ptr[1] = other.ptr[2] = nullptr;
    }

    const char* begin() const { return GetPtr(); }
    const char* end()   const { return UsingHeap() ? ptr[0] + (lowHigh[1] & 0xFFFFFFFFull) : stack; }
    char*       begin()       { return GetPtr(); }
    char*       end()         { return UsingHeap() ? ptr[0] + (lowHigh[1] & 0xFFFFFFFFull) : stack; }

    int  Length() const { return GetSize(); }
    bool Empty()  const { return GetSize() == 0; }
    
          char* CStr()        { return GetPtr(); }
    const char* CStr()  const { return GetPtr(); }

    const char& operator[] (int index) const { ASSERT(index < GetSize() && index > 0); return GetPtr()[index]; }
          char& operator[] (int index)       { ASSERT(index < GetSize() && index > 0); return GetPtr()[index]; }

    String& operator = (const String& right) 
    {
        Asign(right);
        return *this;
    }

    bool operator == (const String& other) { return StringCompare(GetPtr(), other.GetPtr()) == 0; }
    bool operator != (const String& other) { return StringCompare(GetPtr(), other.GetPtr()) != 0; }
    bool operator >  (const String& other) { return StringCompare(GetPtr(), other.GetPtr()) == 1; }
    bool operator <  (const String& other) { return StringCompare(GetPtr(), other.GetPtr()) == -1; }
    bool operator >= (const String& other) { return StringCompare(GetPtr(), other.GetPtr()) >= 0; }
    bool operator <= (const String& other) { return StringCompare(GetPtr(), other.GetPtr()) <= 0; }

    String operator + (const String& other) 
    {
        String cpy(other.GetSize() + GetSize() + 1);
        cpy.Append(GetPtr()); 
        cpy.Append(other.GetPtr()); 
        return (String&&)cpy; 
    }

    String operator + (const char* other)
    {
        String cpy(StringLength(other) + GetSize() + 1);
        cpy.Append(GetPtr()); 
        cpy.Append(other); 
        return (String&&)cpy; 
    }
    
    String& operator += (const String& other) { Append(other.GetPtr()); return *this; }
    String& operator += (const char*   other) { Append(other);       return *this; }

    String operator += (float f) { char arr[16]{}; FloatToString(arr, f); Append(arr); return *this; }
    String operator += (int   i) { char arr[16]{}; IntToString(arr, i);   Append(arr); return *this; }
    String operator +  (float f) { char arr[16]{}; FloatToString(arr, f); String cpy = *this; cpy.Append(arr); return (String&&)cpy; }
    String operator +  (int   i) { char arr[16]{}; IntToString(arr, i);   String cpy = *this; cpy.Append(arr); return (String&&)cpy; }

    void Asign(const String& other)
    {
        // TimeFunction
        if (&other == this) 
            return;
        
        uint32_t otherSize = other.GetSize();
        uint32_t size      = GetSize();

        if (otherSize >= size)
            GrowIfNecessarry(otherSize - size);
        else // size > other.size
            MemSet(GetPtr() + otherSize, 0, size - otherSize); // set rest of the characters to 0

        if (UsingHeap())
        {
            Copy(ptr[0], other.GetPtr(), otherSize);
        }
        else 
        {
            MemCpy<1, 24>(GetPtr(), other.stack, 24);
        }
        SetSize(otherSize);
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
        const char* ptr = GetPtr();
        int size = (int)GetSize();

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
        if (GetSize() == 0)
        {
            Allocate(size);
        }
        else if (GetCapacity() < size)
        {
            Reallocate(GetCapacity(), size);
            if (UsingHeap()) SetCapacity(size);
        }
    }

    void Clear()
    {
        Deallocate();
        this->ptr[0] = this->ptr[1] = this->ptr[2] = nullptr;
    }

    void Insert(int index, char c)
    {
        OpenSpace(index, 1);
        int size = GetSize();
        GetPtr()[size++] = c;
        SetSize(size);
    }

    void Insert(int index, const char* ptr)
    {
        int size = StringLength(ptr);
        OpenSpace(index, size);
        Copy(GetPtr() + index, ptr, size);
    }

    void Insert(int index, const String& other) 
    {
        Insert(index, other.GetPtr()); 
    }
    
    String SubString(int index)
    {
        int count = GetSize() - index;
        ASSERT(index + count <= GetSize());
        String ret(count + 2);
        ret.Append(GetPtr() + index, GetSize() - count);
        return ret;
    }

    String SubString(int index, int count)
    {
        ASSERT(index + count <= GetSize());
        String ret(count + 2);
        ret.Append(GetPtr() + index, count);
        return ret;
    }

    void Remove(int index, int count = 1)
    {
        RemoveSpace(index, count);
        char* ptr = GetPtr();
        int size = (int)GetSize();
        while (count--)
            ptr[--size] = 0;
        SetSize(size);
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
        uint32_t size = GetSize();
        GetPtr()[size++] = c;
        SetSize(size);
    }

    void Append(const char* other, int count)
    {
        GrowIfNecessarry(count);
        uint32_t size = GetSize();
        Copy(GetPtr() + size, other, count);
        SetSize(size + count);
    }

    void Append(const char* other)
    {
        Append(other, StringLength(other));
    }

    void Append(const String& other) { Append(other.GetPtr()); }

    void RemoveSpace(int _index, int _count)
    {
        uint32_t size = this->GetSize();
        ASSERT((_index + _count) <= size);
        int newSize = Max(this->GetSize() - _count, 0u);

        int i = _index;
        int j = _index + _count;

        char* ptr = GetPtr();
        // *******i****j***   < opens space with incrementing both of the pointers (two pointer algorithm)
        while (j < size)
        {
            ptr[i++] = ptr[j++];
        }
    }

    void OpenSpace(int _index, int _count)
    {
        GrowIfNecessarry(_count);
        uint32_t size = this->GetSize();
        long i = Min((long)size + _count, (long)GetCapacity());
        int  j = size;
        ASSERT(i <= INT32_MAX);
        char* ptr = GetPtr();

        while (j >= _index)
        {
            ptr[--i] = ptr[--j];
        }
        SetSize(size + _count);
    }

    void GrowIfNecessarry(int _size)
    {
        uint32_t size = GetSize();
        int newSize = size + _size + 1;
        
        if (newSize <= 23)
            return;

        uint32_t capacity = GetCapacity();
        if (newSize >= capacity)
        {
            constexpr int InitialSize = 48;
            newSize = Max(CalculateArrayGrowth(size + _size), InitialSize);
            
            if (size != 0) 
                Reallocate(size, newSize);
            else
                Allocate(newSize);

            SetCapacity(newSize); 
        }   
    }

    char* c_str()             { return GetPtr(); }
    const char* c_str() const { return GetPtr(); }
#ifdef ASTL_STL_COMPATIBLE
    bool empty()       const  { return GetSize()  == 0; }
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