/********************************************************************************
*    Purpose:                                                                   *
*        Used for manipulating strings.                                         *
*    Author:                                                                    *
*        Anilcan Gulkaya 2023 anilcangulkaya7@gmail.com github @benanil         *
********************************************************************************/
#pragma once

#include "Random.hpp"
#include "Algorithms.hpp"
#include "Memory.hpp"

AX_NAMESPACE 

inline bool StringEqual(const char *a, const char *b, int n) 
{
    for (int i = 0; i < n; i++)
        if (a[i] != b[i])
            return false;
    return true;
}

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

// returns 0 if equal
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

inline char* FindCharInString(const char* s, int c) 
{
    ASSERT(c >= 0 && c < 256);

    __m128i* mem = reinterpret_cast<__m128i*>(const_cast<char*>(s));
    const __m128i set = _mm_setr_epi8(c, 0, 0, 0, 0, 0, 0, 0, 
                                      0, 0, 0, 0, 0, 0, 0, 0);

    const uint8_t mode = _SIDD_UBYTE_OPS | _SIDD_CMP_EQUAL_ANY | _SIDD_LEAST_SIGNIFICANT;
        
    for (/**/; /**/; mem++) {

        const __m128i chunk = _mm_loadu_si128(mem);

        if (_mm_cmpistrc(set, chunk, mode)) {
            // there is character c in a chunk
            const int idx = _mm_cmpistri(set, chunk, mode);

            return reinterpret_cast<char*>(mem) + idx;;
        } else if (_mm_cmpistrz(set, chunk, mode)) {
            // there is zero byte in a chunk
            break;
        }
    }

    return nullptr;
}

#else
inline int StringLength(const char* s)
{
    const char* begin = s;
    while (*s) s++;
    return s - begin;
}

inline int StringCompare(const char* a, const char* b)
{
    for (; *a && *b; a++, b++)
    {
        if (*a != *b)
        {
            if (*a < *b) return -1;
            else         return +1; // greater
        }
    }
    return !(*a == *b);// strings are equal
}

inline const char* FindCharInString(const char *s, int c)
{
    int idx = 0;
    while (s[idx])
        if (s[idx++] == c)
            return s + (--idx);
    return nullptr;
}

#endif

// https://github.com/lemire/fastvalidate-utf-8/blob/master/include/simdasciicheck.h
inline bool ValidateAscii(const char *src, uint64_t len) 
{
    uint64_t i = 0;
    int error_mask = 0;
#ifdef AX_SUPPORT_SSE
    __m128i has_error = _mm_setzero_si128();
    if (len >= 16) 
    {
        for (; i <= len - 16; i += 16) 
        {
            has_error = _mm_or_si128(has_error, _mm_loadu_si128((const __m128i *)(src + i)));
        }
    }
    error_mask = _mm_movemask_epi8(has_error);
#endif

    char tail_has_error = 0;
    for (; i < len; i++) 
        tail_has_error |= src[i];
    
    error_mask |= (tail_has_error & 0x80);
    return !error_mask;
}

inline bool IsUTF8ASCII(const char* string, uint64_t size)
{
    const unsigned char * bytes = (const unsigned char *)string;
    for (int i = 0; i < size; i++)
    {
        // use bytes[0] <= 0x7F to allow ASCII control characters
        if (!(bytes[i] == 0x09 || bytes[i] == 0x0A || bytes[i] == 0x0D ||
              (0x20 <= bytes[i] && bytes[i] <= 0x7E)))
            return false;
    }
    return true;
}

inline bool IsUTF8(char c) 
{
    return (c & 0xC0) != 0x80; 
}

inline unsigned UTF8NextChar(char *s, int *i)
{
    unsigned ch = 0;
    int sz = 0;
    static const unsigned offsetsFromUTF8[6] = {
        0x00000000U, 0x00003080U, 0x000E2080U,
        0x03C82080U, 0xFA082080U, 0x82082080U
    };

    do {
        ch <<= 6;
        ch += (unsigned char)s[(*i)++];
        sz++;
    } while (s[*i] && !IsUTF8(s[*i]));
    ch -= offsetsFromUTF8[sz-1];
    return ch;
}

// small string optimization
class String
{
public:
    
    uint64_t  lowHigh[3]{}; // when using heap high part is =  capacity | (size << 31);

    static const uint32_t LastBit32 = 1u << 31u;
    static const uint64_t  LastBit64 = 1ull << 63ull;

    uint32_t GetCapacity() const
    {
        return UsingHeap() ? uint32_t(lowHigh[1] >> 32ull) & ~LastBit32 : 23;
    }
    
    void SetCapacity(int cap) 
    {
        if (!UsingHeap()) return;
        lowHigh[1] &= 0xFFFFFFFFull | LastBit64;
        lowHigh[1] |= uint64_t (cap) << 32ull;
    }

    uint32_t GetSize() const 
    {
        return !UsingHeap() ? StringLength((char*)lowHigh) : (uint32_t)(lowHigh[1] & 0xFFFFFFFFull);
    }
    
    void SetSize(uint32_t size) 
    {
        if (!UsingHeap()) return;
        lowHigh[1] &= 0xFFFFFFFFull << 32ull; 
        lowHigh[1] |= size;
    }

    bool UsingHeap() const { return lowHigh[1] &  LastBit64; }
    void SetUsingHeap()    {        lowHigh[1] |= LastBit64; }

    const char* GetPtr() const { return UsingHeap() ? (char*)lowHigh[0] : (char*)lowHigh; }
          char* GetPtr()       { return UsingHeap() ? (char*)lowHigh[0] : (char*)lowHigh; }

    void Allocate(int size)
    {
        if (size < 23)
            return;

        int oldSize = GetSize();
        char* ptr = new char[size]{};
        lowHigh[0] = (uint64_t )ptr;
        lowHigh[1] = oldSize;
        lowHigh[2] = 0;

        SetUsingHeap();
    }
    
    void Deallocate()
    {
        if (UsingHeap())
            delete[] GetPtr();

        lowHigh[0] = lowHigh[1] = lowHigh[2] = 0;
    }

    void Reallocate(int oldCount, int count)
    {
        if (UsingHeap())
        {
            char* newPtr = new char[count] {};
            char* ptr = (char*)this->lowHigh[0];
            MemCpy<1>(newPtr, ptr, oldCount);
            delete[] ptr;
            this->lowHigh[0] = (uint64_t )newPtr;
        }
        else if (count > 23) // was stack allocating but bigger memory requested
        {
            char* newPtr = new char[count] {};
            SmallMemCpy(newPtr, lowHigh, oldCount);
            lowHigh[0] = (uint64_t )newPtr;
            lowHigh[2] = 0;
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
        SmallMemCpy(lowHigh, other.lowHigh, 24);
        other.lowHigh[0] = other.lowHigh[1] = other.lowHigh[2] = 0;
    }

    const char* begin() const { return GetPtr(); }
    const char* end()   const { return UsingHeap() ? (const char*)(lowHigh[0] + (lowHigh[1] & 0xFFFFFFFFull)) : (char*)lowHigh; }
    char*       begin()       { return GetPtr(); }
    char*       end()         { return UsingHeap() ? (char*)(lowHigh[0] + (lowHigh[1] & 0xFFFFFFFFull)) : (char*)lowHigh; }

    int  Length() const { return GetSize(); }
    bool Empty()  const { return GetSize() == 0; }
    
          char* CStr()        { return GetPtr(); }
    const char* CStr()  const { return GetPtr(); }

    const char& operator[] (uint32_t index) const { ASSERT(index < GetSize() && index > 0); return GetPtr()[index]; }
          char& operator[] (uint32_t index)       { ASSERT(index < GetSize() && index > 0); return GetPtr()[index]; }

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
            Copy((char*)lowHigh[0], other.GetPtr(), otherSize);
        }
        else 
        {
            SmallMemCpy(GetPtr(), other.lowHigh, 24);
        }
        SetSize(otherSize);
    }

    void Append(float f) { char arr[16]{}; FloatToString(arr, f); Append(arr); }
    void Append(int i)   { char arr[16]{}; IntToString(arr, i);   Append(arr); }
    

    int IndexOf(const char* other) const
    {
        int otherLen = StringLength(other);
        const char* ptr = GetPtr();
        int size = (int)GetSize();

        for (int i = size - otherLen; i >= 0; i--) {
            if (StringEqual(ptr + i, other, otherLen))
                return i;
        }
        return -1;
    }

    int IndexOf(const String& other)   const { return IndexOf(other.CStr()) != -1; }
    bool Contains(const char* other)   const { return IndexOf(other) != -1; }
    bool Contains(const String& other) const { return IndexOf(other.CStr()) != -1; }

    void Reserve(uint32_t size)
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
        lowHigh[0] = lowHigh[1] = lowHigh[2] = 0;
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
    
    String SubString(uint32_t index)
    {
        uint32_t count = GetSize() - index;
        ASSERT(index + count <= GetSize());
        String ret(count + 2);
        ret.Append(GetPtr() + index, GetSize() - count);
        return ret;
    }

    String SubString(uint32_t index, int count)
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

    void RemoveSpace(uint32_t _index, int _count)
    {
        uint32_t size = this->GetSize();
        ASSERT((_index + _count) <= size);
        uint32_t newSize = MAX(this->GetSize() - _count, 0u);

        uint32_t i = _index;
        uint32_t j = _index + _count;

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
        long i = MIN((long)size + _count, (long)GetCapacity());
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
        uint32_t newSize = size + _size + 1;
        
        if (newSize <= 23)
            return;

        uint32_t capacity = GetCapacity();
        if (newSize >= capacity)
        {
            const int InitialSize = 48;
            newSize = MAX(CalculateArrayGrowth(size + _size), InitialSize);
            
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
    static __forceinline uint64_t Hash(const String& x)
    {
        return WYHash::Hash(x.CStr(), x.Length());
        // return MurmurHash64(x.CStr(), x.Length(), 0xa0761d6478bd642full);
    }
};

struct StringView
{
  const char* ptr;
  int size;

  StringView(const char* _ptr)
  : ptr(_ptr), size(StringLength(ptr))
  { }

  StringView(const String& string)
  : ptr(string.GetPtr()), size(string.GetSize())
  { }

  char operator[](int index)
  {
    return ptr[index];
  }
};


#if 1 // AX_SUPPORT_SSE
/* based on the valid_utf8 routine from the PCRE library by Philip Hazel
length is in bytes, since without knowing whether the string is valid
it's hard to know how many characters there are! */
// returns 1 if ASCII UTF8, returns 2 if non ASCII UTF8, if not utf8 returns 0
inline int UTF8Valid(const char *str, uint64_t length)
{
    const unsigned char *p, *pend = (unsigned char*)str + length;
    unsigned char c;
    int ret = 1; /* ASCII */
    uint64_t ab;

    static const char trailingBytesForUTF8[128] = {
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
        2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, 3,3,3,3,3,3,3,3,4,4,4,4,5,5,5,5
    };

    for (p = (unsigned char*)str; p < pend; p++) {
        c = *p;
        if (c < 128)
            continue;
        ret = 2; /* non-ASCII UTF-8 */
        if ((c & 0xc0) != 0xc0)
            return 0;
        ab = trailingBytesForUTF8[c-128];
        if (length < ab)
            return 0;
        length -= ab;

        p++;
        /* Check top bits in the second byte */
        if ((*p & 0xc0) != 0x80)
            return 0;

        /* Check for overlong sequences for each different length */
        switch (ab) {
            case 1:
                /* Check for xx00 000x */
                if ((c & 0x3e) == 0) return 0;
                continue;   /* We know there aren't any more bytes to check */
            case 2:
                /* Check for 1110 0000, xx0x xxxx */
                if (c == 0xe0 && (*p & 0x20) == 0) return 0;
                break;

            case 3:
                /* Check for 1111 0000, xx00 xxxx */
                if (c == 0xf0 && (*p & 0x30) == 0) return 0;
                break;
            case 4:
                /* Check for 1111 1000, xx00 0xxx */
                if (c == 0xf8 && (*p & 0x38) == 0) return 0;
                break;
            case 5:
                /* Check for leading 0xfe or 0xff, and then for 1111 1100, xx00 00xx */
                if (c == 0xfe || c == 0xff || (c == 0xfc && (*p & 0x3c) == 0)) return 0;
                break;
        }

        /* Check for valid bytes after the 2nd, if any; all must start 10 */
        while (--ab > 0) {
            if ((*(++p) & 0xc0) != 0x80) return 0;
        }
    }
    return ret;
}
#else
// all byte values must be no larger than 0xF4
inline void checkSmallerThan0xF4(__m128i current_bytes, __m128i *has_error) 
{
    // unsigned, saturates to 0 below max
    *has_error = _mm_or_si128(*has_error, _mm_subs_epu8(current_bytes, _mm_set1_epi8(0xF4)));
}

inline __m128i continuationLengths(__m128i high_nibbles) {
    return _mm_shuffle_epi8(
        _mm_setr_epi8(1, 1, 1, 1, 1, 1, 1, 1, // 0xxx (ASCII)
                      0, 0, 0, 0,             // 10xx (continuation)
                      2, 2,                   // 110x
                      3,                      // 1110
                      4), // 1111, next should be 0 (not checked here)
        high_nibbles);
}

inline __m128i carryContinuations(__m128i initial_lengths, __m128i previous_carries) 
{
    __m128i right1 = _mm_subs_epu8(_mm_alignr_epi8(initial_lengths, previous_carries, 16 - 1), _mm_set1_epi8(1));
    __m128i sum = _mm_add_epi8(initial_lengths, right1);
    __m128i right2 = _mm_subs_epu8(_mm_alignr_epi8(sum, previous_carries, 16 - 2), _mm_set1_epi8(2));
    return _mm_add_epi8(sum, right2);
}

inline void checkContinuations(__m128i initial_lengths, __m128i carries, __m128i *has_error) {
    // overlap || underlap
    // carry > length && length > 0 || !(carry > length) && !(length > 0)
    // (carries > length) == (lengths > 0)
    __m128i overunder = _mm_cmpeq_epi8(_mm_cmpgt_epi8(carries, initial_lengths), _mm_cmpgt_epi8(initial_lengths, _mm_setzero_si128()));
    *has_error = _mm_or_si128(*has_error, overunder);
}

// when 0xED is found, next byte must be no larger than 0x9F
// when 0xF4 is found, next byte must be no larger than 0x8F
// next byte must be continuation, ie sign bit is set, so signed < is ok
inline void checkFirstContinuationMax(__m128i current_bytes, __m128i off1_current_bytes, __m128i *has_error) {
    __m128i maskED = _mm_cmpeq_epi8(off1_current_bytes, _mm_set1_epi8(0xED));
    __m128i maskF4 = _mm_cmpeq_epi8(off1_current_bytes, _mm_set1_epi8(0xF4));
    __m128i badfollowED = _mm_and_si128(_mm_cmpgt_epi8(current_bytes, _mm_set1_epi8(0x9F)), maskED);
    __m128i badfollowF4 = _mm_and_si128(_mm_cmpgt_epi8(current_bytes, _mm_set1_epi8(0x8F)), maskF4);
    *has_error = _mm_or_si128(*has_error, _mm_or_si128(badfollowED, badfollowF4));
}

// map off1_hibits => error condition
// hibits     off1    cur
// C       => < C2 && true
// E       => < E1 && < A0
// F       => < F1 && < 90
// else      false && false
inline void checkOverlong(__m128i current_bytes, __m128i off1_current_bytes, __m128i hibits, __m128i previous_hibits, __m128i *has_error) {
    __m128i off1_hibits = _mm_alignr_epi8(hibits, previous_hibits, 16 - 1);
    __m128i initial_mins = _mm_shuffle_epi8(
                           _mm_setr_epi8(-128, -128, -128, -128, -128, -128, -128, -128, -128, -128,
                                         -128, -128, // 10xx => false
                                         0xC2, -128, // 110x
                                         0xE1,       // 1110
                                         0xF1),off1_hibits);

    __m128i initial_under = _mm_cmpgt_epi8(initial_mins, off1_current_bytes);

    __m128i second_mins = _mm_shuffle_epi8(
                          _mm_setr_epi8(-128, -128, -128, -128, -128, -128, -128, -128, -128, -128,
                                        -128, -128, // 10xx => false
                                        127, 127,   // 110x => true
                                        0xA0,       // 1110
                                        0x90), off1_hibits);
    __m128i second_under = _mm_cmpgt_epi8(second_mins, current_bytes);
    *has_error = _mm_or_si128(*has_error, _mm_and_si128(initial_under, second_under));
}

struct processed_utf_bytes 
{
    __m128i rawbytes;
    __m128i high_nibbles;
    __m128i carried_continuations;
};

inline void count_nibbles(__m128i bytes, processed_utf_bytes *answer) 
{
    answer->rawbytes = bytes;
    answer->high_nibbles = _mm_and_si128(_mm_srli_epi16(bytes, 4), _mm_set1_epi8(0x0F));
}

// check whether the current bytes are valid UTF-8
// at the end of the function, previous gets updated
inline processed_utf_bytes checkUTF8Bytes(__m128i current_bytes, processed_utf_bytes* previous, __m128i *has_error) 
{
    struct processed_utf_bytes pb;
    count_nibbles(current_bytes, &pb);

    checkSmallerThan0xF4(current_bytes, has_error);

    __m128i initial_lengths = continuationLengths(pb.high_nibbles);
    pb.carried_continuations = carryContinuations(initial_lengths, previous->carried_continuations);

    checkContinuations(initial_lengths, pb.carried_continuations, has_error);

    __m128i off1_current_bytes = _mm_alignr_epi8(pb.rawbytes, previous->rawbytes, 16 - 1);
    checkFirstContinuationMax(current_bytes, off1_current_bytes, has_error);

    checkOverlong(current_bytes, off1_current_bytes, pb.high_nibbles, previous->high_nibbles, has_error);
    return pb;
}

inline bool UTF8Valid(const char *src, uint64_t len) 
{
    uint64_t i = 0;
    __m128i has_error = _mm_setzero_si128();
    processed_utf_bytes previous;
    previous.rawbytes = _mm_setzero_si128();
    previous.high_nibbles = _mm_setzero_si128();
    previous.carried_continuations = _mm_setzero_si128();
    
    if (len >= 16) {
        for (; i <= len - 16; i += 16) {
            __m128i current_bytes = _mm_loadu_si128((const __m128i *)(src + i));
            previous = checkUTF8Bytes(current_bytes, &previous, &has_error);
        }
    }

    // last part
    if (i < len) {
        char buffer[16]{};
        SmallMemCpy(buffer, src + i, len - i);
        __m128i current_bytes = _mm_loadu_si128((const __m128i *)(buffer));
        previous = checkUTF8Bytes(current_bytes, &previous, &has_error);
    } else {
        has_error =
        _mm_or_si128(_mm_cmpgt_epi8(previous.carried_continuations,
                                    _mm_setr_epi8(9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 1)),
                                    has_error);
    }

    return _mm_testz_si128(has_error, has_error);
}

#endif // has sse

AX_END_NAMESPACE 

