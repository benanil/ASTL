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

inline void MemSetAligned64(void* dst, char val, uint64_t sizeInBytes)
{
    const uint64_t* end = (uint64_t*)((char*)dst + sizeInBytes) - 4;
    uint64_t* dp = (uint64_t*)dst;
    uint64_t  d4 = (uint64_t)val;
    d4 |= d4 << 8ull; d4 |= d4 << 16ull; d4 |= d4 << 32ull;

    while (dp < end)
    {
        dp[0] = dp[1] = dp[2] = dp[3] = d4;
        dp += 4;
    }
    end += 4;
    while (dp < end) *dp++ = d4;
}

inline void MemSetAligned32(void* dst, char val, uint64_t sizeInBytes)
{
    const uint32_t* end = (uint32_t*)((char*)dst + sizeInBytes) - 4;
    uint32_t* dp = (uint32_t*)dst;
    uint32_t  d4 = (uint32_t)val;
    d4 |= d4 << 8ull; d4 |= d4 << 16ull;

    while (dp < end)
    {
        dp[0] = dp[1] = dp[2] = dp[3] = d4;
        dp += 4;
    }
    end += 4;
    while (dp < end) *dp++ = d4;
}

template<int alignment = 0, int size = 0>
inline void MemSet(void* dst, char val, uint64_t sizeInBytes)
{
    if constexpr (size)
    // use this for struct/class types such as Vector3 and Matrix4, and use MemSet for big arrays or unknown size arrays
    #ifdef _MSC_VER
        __stosb((unsigned char*)dst, (unsigned char)val, size);
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
            while (sizeInBytes)
            {
                *dp++ = val;
                sizeInBytes--;
            }
        }
    }
}

inline void MemCpyAligned64(void* dst, const void* src, uint64_t sizeInBytes)
{
    uint64_t*       dp  = (uint64_t*)dst;
    const uint64_t* sp  = (const uint64_t*)src;
    const uint64_t* end = (const uint64_t*)(src) + (sizeInBytes / 8) - 4;
        
    while (sp < end)
    {
        dp[0] = sp[0];
        dp[1] = sp[1];
        dp[2] = sp[2];
        dp[3] = sp[3];
        dp += 4, sp += 4;
    }
    end += 4;
    while (sp < end) *dp++ = *sp++;
}

inline void MemCpyAligned32(void* dst, const void* src, uint64_t sizeInBytes)
{
    uint32_t*       dp  = (uint32_t*)dst;
    const uint32_t* sp  = (const uint32_t*)src;
    const uint32_t* end = (const uint32_t*)(src) + (sizeInBytes / 4) - 4;
        
    while (sp < end)
    {
        dp[0] = sp[0];
        dp[1] = sp[1];
        dp[2] = sp[2];
        dp[3] = sp[3];
        dp += 4, sp += 4;
    }
    end += 4;
    while (sp < end) *dp++ = *sp++;
}

template<int alignment = 0, int size = 0>
inline void MemCpy(void* dst, const void* src, uint64_t sizeInBytes)
{
    // use this for struct/class types such as Vector3 and Matrix4,
    // and use MemCpy for big arrays or unknown size arrays
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
        const uint64_t alignXor = uint64_t(dst) ^ uint64_t(src);
        if (!(alignXor & 7))
            MemCpyAligned64(dst, src, sizeInBytes);
        else if (!(alignXor & 3))
            MemCpyAligned32(dst, src, sizeInBytes);
        else
        {
            const char* cend = (char*)(src + sizeInBytes);
            const char* scp = (const char*)src;
            char* dcp = (char*)dst;
        
            while (scp < cend) *dcp++ = *scp++;
        }
    }
}

// todo and memmove

template<typename T> 
struct Allocator 
{
    static constexpr bool IsPOD = false;

    T* Allocate(int count) {
        return new T[count]{};
    }
    
    template<typename ...Args>
    T* Allocate(Args&... args)
    {
        return new T(Forward<Args>(args)...);
    }

    T* AllocateUninitialized(int count) {
        return (T*)::operator new(sizeof(T) * count); // for now I close uninitialized 
    }
    
    void Deallocate(T* ptr, int count) {
        ::operator delete(ptr, count * sizeof(T));
    }
    
    T* Reallocate(T* ptr, int oldCount, int count)
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

    T* Allocate(int count) {
        return new T[count];
    }
      
    template<typename ...Args>
    T* Allocate(Args&... args)
    {
        return new T(Forward<Args>(args)...);
    }

    T* AllocateUninitialized(int count) {
        return new T[count];
    }
      
    void Deallocate(T* ptr, int count) {
        delete[] ptr;
    }
      
    T* Reallocate(T* ptr, int oldCount, int count)
    {
        T* old = ptr;
        T* nev = AllocateUninitialized(count);
        MemCpy<alignof(T)>(nev, old, count);
        Deallocate(old, oldCount);
        return nev;
    }
};

template<typename T>
struct FixedSizeGrowableAllocator
{
    static constexpr int InitialSize = 512 / Min((int)sizeof(T), 128);

    struct Fragment
    {
        Fragment* next;
        T*        ptr;
    };

    int currIndex     = 0;
    int currentSize   = 0;
    Fragment* base    = nullptr;
    Fragment* current = nullptr;

    FixedSizeGrowableAllocator()
    {
        currentSize = InitialSize;
        currIndex = 0;
        base = new Fragment();
        current = base;
        base->next  = nullptr;
        base->ptr   = new T[InitialSize];
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

    void CheckFixGrow(int count)
    {
        if (currIndex + count >= currentSize)
        {
            currIndex = 0;
            currentSize <<= 1;
            current = new Fragment();
            current->next = nullptr;
            current->ptr = new T[currentSize];
        }
    }

    template<typename ...Args>
    T* Allocate(Args&&... args)
    {
        CheckFixGrow(1);
        T* ptr = current->ptr + currIndex;
        new (ptr) T(Forward<Args>(args)...);
        currIndex++;
        return ptr;
    }

    T* Allocate(int count)
    {
        CheckFixGrow(count);
        T* ptr     = current->ptr + currIndex;
        currIndex += count;
        for (int i = 0; i < count; i++)
            ptr[i].T();
        return ptr;
    }

    T* AllocateUninitialized(int count)
    {
        CheckFixGrow(count);
        T* ptr     = current->ptr + currIndex;
        currIndex += count;
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