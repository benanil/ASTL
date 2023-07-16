
#pragma once

#include "Common.hpp"
 
template<typename T> struct RemoveRef      { typedef T Type; };
template<typename T> struct RemoveRef<T&>  { typedef T Type; };
template<typename T> struct RemoveRef<T&&> { typedef T Type; };
template<typename T> struct RemovePtr      { typedef T Type; };
template<typename T> struct RemovePtr<T*>  { typedef T Type; };

template<typename T>
FINLINE typename RemoveRef<T>::Type&& Move(T&& obj)
{
  typedef typename RemoveRef<T>::Type CastType;
  return (CastType&&)obj;
}

template<typename T>
FINLINE T&& Forward(typename RemoveRef<T>::Type& obj) { return (T&&)obj; }

template<typename T>
FINLINE T&& Forward(typename RemoveRef<T>::Type&& obj) { return (T&&)obj; }

template<class T, class U = T>
FINLINE constexpr T Exchange(T& obj, U&& new_value)
{
  T old_value = Move(obj);
  obj = Forward<U>(new_value);
  return old_value;
}

// Aligns the size by the machine word.
inline uint64_t Align(uint64_t n, uint64_t alignment = 64) 
{
    return (n + alignment - 1) & ~(alignment - 1);
}

// if you want simd version of memmove it is here: https://hackmd.io/@AndybnA/0410
inline void* MemMove(void *dest, const void *src, size_t n)
{
    char* d       = (char*) dest; 
    const char* s = (const char*)src;

    if (dest >= src && dest < s + n) /* overlap */ 
        for (int i = n - 1; i >= 0; i--) 
            d[i] = s[i];
    else 
        for (size_t i = 0; i < n; i++)
            d[i] = s[i];
    
    return dest;
}

inline void MemSetAligned64(void* dst, unsigned char val, uint64_t sizeInBytes)
{
    // we want an offset because the while loop below iterates over 4 uint64 at a time
    const uint64_t offset = (sizeInBytes >> 5) & 31;  // 5 = divide by 32 because 4 x uint64_t, & 31 is mod by 32
    const uint64_t* end = (uint64_t*)((char*)dst + sizeInBytes - offset);
    uint64_t* dp = (uint64_t*)dst;
    uint64_t  d4 = (uint64_t)val;
#ifdef _MSC_VER
    __stosb((unsigned char*)&d4, val, 8);
#else
    d4 |= d4 << 8ull; d4 |= d4 << 16ull; d4 |= d4 << 32ull;
#endif

    while (dp < end)
    {
        dp[0] = dp[1] = dp[2] = dp[3] = d4;
        dp += 4;
    }

    end = (uint64_t*)((char*)dst + sizeInBytes);
    while (dp < end)
        *dp++ = d4;
}

inline void MemSetAligned32(void* dst, unsigned char val, uint64_t sizeInBytes)
{
    const uint32_t offset = (sizeInBytes >> 4) & 15; // 4 = divide by 16 because 4 x uint32_t, & 15 is mod by 16
    const uint32_t* end = (uint32_t*)((char*)dst + sizeInBytes - offset);
    uint32_t* dp = (uint32_t*)dst;
    uint32_t  d4 = (uint32_t)val;
#ifdef _MSC_VER
    __stosb((unsigned char*)&d4, val, 4);
#else
    d4 |= d4 << 8ull; d4 |= d4 << 16ull;
#endif
    while (dp < end)
    {
        dp[0] = dp[1] = dp[2] = dp[3] = d4;
        dp += 4;
    }

    end = (uint32_t*)((char*)dst + sizeInBytes);
    while (dp < end)
        *dp++ = d4;
}
#include <memory.h>
// use size for struct/class types such as Vector3 and Matrix4, 
// leave as zero size constant for big arrays or unknown size arrays
template<int alignment = 0, int size = 0>
inline void MemSet(void* dst, unsigned char val, uint64_t sizeInBytes)
{
                                       memset(dst, val, sizeInBytes);
                                       return;
    if constexpr (size)
    #ifdef _MSC_VER
        __stosb((unsigned char*)dst, val, size);
    #else
        __builtin_memset(dst, val, sizeInBytes);
    #endif
    else if (alignment == 8)
        MemSetAligned64(dst, val, sizeInBytes);
    else if (alignment == 4)
        MemSetAligned32(dst, val, sizeInBytes);
    else
    {
        if ((((uint64_t)dst) & 7) == 0)       MemSetAligned64(dst, val, sizeInBytes);
        else if ((((uint64_t)dst) & 3) == 0)  MemSetAligned32(dst, val, sizeInBytes);
        else
        {
            char* dp = (char*)dst;
            while (sizeInBytes--)
                *dp++ = val;
        }
    }
}

inline void MemCpyAligned64(void* dst, const void* src, uint64_t sizeInBytes)
{
    // 5 = divide by 32 because 4 x uint64_t, & 31 is mod by 32
    const uint64_t offset = (sizeInBytes >> 5) & 31; 
    uint64_t*       dp  = (uint64_t*)dst;
    const uint64_t* sp  = (const uint64_t*)src;
    const uint64_t* end = (const uint64_t*)((char*)src + sizeInBytes - offset);
        
    while (sp < end)
    {
        dp[0] = sp[0];
        dp[1] = sp[1];
        dp[2] = sp[2];
        dp[3] = sp[3];
        dp += 4, sp += 4;
    }
    end = (const uint64_t*)((char*)src + sizeInBytes);

    while (dp < end)
        *dp++ = *sp++;
}

inline void MemCpyAligned32(void* dst, const void* src, uint64_t sizeInBytes)
{
    const uint32_t offset = (sizeInBytes >> 4) & 15; 
    uint32_t*       dp  = (uint32_t*)dst;
    const uint32_t* sp  = (const uint32_t*)src;
    const uint32_t* end = (const uint32_t*)((char*)src + sizeInBytes - offset);
        
    while (sp < end)
    {
        dp[0] = sp[0];
        dp[1] = sp[1];
        dp[2] = sp[2];
        dp[3] = sp[3];
        dp += 4, sp += 4;
    }
    
    end = (const uint32_t*)((char*)src + sizeInBytes);
    while (dp < end)
        *dp++ = *sp++;
}

// use size for struct/class types such as Vector3 and Matrix4,
// and use MemCpy for big arrays or unknown size arrays
template<int alignment = 0, int size = 0>
inline void MemCpy(void* dst, const void* src, uint64_t sizeInBytes)
{
                                       memcpy(dst, src, sizeInBytes);
                                       return;
    if constexpr (size != 0)
#ifdef _MSC_VER
    __movsb((unsigned char*)dst, (unsigned char const*)src, size);
#else
    __builtin_memcpy(dst, src, size);
#endif
    else if (alignment == 8)
        MemCpyAligned64(dst, src, sizeInBytes);
    else if (alignment == 4)
        MemCpyAligned32(dst, src, sizeInBytes);
    else
    {
        const uint64_t alignXor = uint64_t(dst) ^ uint64_t(src); // to check if both has same alignment
        if (!(alignXor & 7))
            MemCpyAligned64(dst, src, sizeInBytes);
        else if (!(alignXor & 3))
            MemCpyAligned32(dst, src, sizeInBytes);
        else
        {
            const char* cend = (char*)((char*)src + sizeInBytes);
            const char* scp = (const char*)src;
            char* dcp = (char*)dst;
        
            while (scp < cend) *dcp++ = *scp++;
        }
    }
}


// todo remove template from allocators

template<typename T> 
struct Allocator 
{
    static constexpr bool IsPOD = false;

    T* Allocate(int count) const {
        return new T[count]{};
    }
    
    T* AllocateUninitialized(int count) const {
        return new T[count]; 
    }
    
    void Deallocate(T* ptr, int count) const {
        delete[] ptr;
    }
    
    T* Reallocate(T* ptr, int oldCount, int count) const
    {
        T* old  = ptr;
        T* _new = new T[count];
        while (oldCount > count) // delete remaining memory if new count is less than count (reduced)
        {
            old[--oldCount].~T();
        }
        for (int i = 0; i < Min(count, oldCount); ++i)
        {
            _new[i] = Forward<T>(old[i]);
        }
        delete[] old;
        return _new;
    }
};

template<typename T>
struct MallocAllocator
{
    static constexpr bool IsPOD = true;

    T* Allocate(int count) const {
        return new T[count]{};
    }

    T* AllocateUninitialized(int count) const {
        return new T[count];
    }
      
    void Deallocate(T* ptr, int count) const {
        delete[] ptr;
    }
      
    T* Reallocate(T* ptr, int oldCount, int count) const
    {
        T* old = ptr;
        T* nev = AllocateUninitialized(count);
        MemCpy<alignof(T)>(nev, old, oldCount * sizeof(T));
        Deallocate(old, oldCount);
        return nev;
    }
};

template<typename T>
struct FixedSizeGrowableAllocator
{
    static constexpr bool IsPOD = false;

    static constexpr int InitialSize = NextPowerOf2(512 / Min((int)sizeof(T), 128));

    struct Fragment
    {
        Fragment* next;
        T*        ptr ;
        int64_t   size; // used like a index until we fill the fragment
    };

    int currentCapacity   = 0;
    Fragment* base    = nullptr;
    Fragment* current = nullptr;

    FixedSizeGrowableAllocator()
    {
        currentCapacity = InitialSize;
        base = new Fragment;
        current = base;
        base->next  = nullptr;
        base->ptr   = new T[InitialSize];
        base->size  = 0;
    }

    ~FixedSizeGrowableAllocator()
    {
        while (base)
        {
            delete[] base->ptr;
            Fragment* oldBase = base;
            base = base->next;
            delete oldBase;
        }
    }

    // copy constructor. no need move constructor
    FixedSizeGrowableAllocator(const FixedSizeGrowableAllocator& other)
    {
        if (!other.base) return;

        int64_t totalSize = 0l;
        Fragment* start = other.base;

        while (start)
        {
            totalSize += start->size;
            start = start->next;
        }

        currentCapacity = NextPowerOf2(totalSize);
        
        // even though other contains multiple fragments we will fit all data into one fragment
        base = new Fragment;
        base->next = nullptr;
        base->ptr  = new T[currentCapacity]; 
        base->size = totalSize;
        current = base;
        
        // copy other's memory to ours
        T* curr = base->ptr;
        start   = other.base;
        while (start)
        {
            Copy(curr, start->ptr, start->size);
            curr += start->size;
            start = start->next;
        }
    }

    void CheckFixGrow(int count)
    {
        if (current->size + count >= currentCapacity)
        {
            currentCapacity <<= 1;
            current = new Fragment();
            current->next = nullptr;
            current->ptr  = new T[currentCapacity];
            current->size = 0;
        }
    }

    T* Allocate(int count)
    {
        CheckFixGrow(count);
        T* ptr = current->ptr + current->size;
        current->size += count;
        for (int i = 0; i < count; i++)
            ptr[i].T();
        return ptr;
    }

    T* AllocateUninitialized(int count)
    {
        CheckFixGrow(count);
        T* ptr = current->ptr + current->size;
        current->size += count;
        return ptr;
    }

    void Deallocate(T* ptr, int count)
    {
        for (int i = 0; i < count; i++)
            ptr[i].~T();
    }

    T* Reallocate(T* ptr, int oldCount, int count)
    {
        Deallocate(ptr, oldCount);
        return Allocate(count);
    }
};

// todo add ScopedPtr
// todo add ScopedFn
// todo SharedPtr in different hpp file 