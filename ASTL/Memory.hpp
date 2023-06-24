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

#ifdef AX_SUPPORT_AVX2
AXGLOBALCONST uint64_t g_256MemMask[8] = { ~0ull, ~0ull, ~0ull, ~0ull, 0ull, 0ull, 0ull, 0ull };
#elif AX_SUPPORT_SSE2
AXGLOBALCONST uint32_t g_128MemMask[8] = { ~0u, ~0u, ~0u, ~0u, 0u, 0u, 0u, 0u };
#endif

void MemSet(void* dst, char val, uint64_t sizeInBytes)
{
    if (sizeInBytes <= 512)
    {
#ifdef AX_SUPPORT_AVX2
        __stosb((unsigned char*)dst, BitCast<unsigned char>(val), sizeInBytes);
#else
        __builtin_memset(dst, val, sizeInBytes);
#endif
        return;
    }

    uint64_t* dp = (uint64_t*)dst;
    uint64_t  d4 = (uint64_t)val;
    d4 |= d4 << 8ull; d4 |= d4 << 16ull; d4 |= d4 << 32ull;

    while (sizeInBytes >= (sizeof(uint64_t) * 4))
    {
        dp[0] = dp[1] = dp[2] = dp[3] = d4;
        dp += 4;
        sizeInBytes -= sizeof(uint64_t) * 4;
    }
    
    char* dcp = (char*)dp;
    while (sizeInBytes)
    {
        *dcp++ = val;
        sizeInBytes--;
    }
}

inline void MemCpy(void* dst, const void* src, uint64_t sizeInBytes)
{
    if (sizeInBytes <= 512) // if data is small we copy it without loop
    {
#ifdef _MSC_VER
        __movsb((unsigned char*)dst, (unsigned char const*)src, sizeInBytes);
#else
        __builtin_memcpy(dst, src, size);
#endif
        return;
    }
    uint64_t* dp = (uint64_t*)dst;
    const uint64_t* sp = (const uint64_t*)src;

    while (sizeInBytes >= (sizeof(uint64_t) * 4))
    {
        dp[0] = sp[0];
        dp[1] = sp[1];
        dp[2] = sp[2];
        dp[3] = sp[3];
        dp += 4; sp += 4;
        sizeInBytes -= sizeof(uint64_t) * 4;
    }

    char* dcp = (char*)dp;
    const char* scp = (const char*)sp;

    while (sizeInBytes)
    {
        *dcp++ = *scp++;
        sizeInBytes--;
    }
}

template<typename T> 
struct Allocator 
{
    T* Allocate(int count) {
        return new T[count]{};
    }
    
    T* AllocateUninitialized(int count) {
        // return (T*)::operator new(sizeof(T) * count); // for now I close uninitialized 
        return new T[count]{};
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
    T* Allocate(int count) {
        return new T[count];
    }
      
    T* AllocateUninitialized(int count) {
        return (T*)::operator new(sizeof(T) * count);
    }
      
    void Deallocate(T* ptr, int count) {
        ::operator delete(ptr, count * sizeof(T));
    }
      
    T* Reallocate(T* ptr, int oldCount, int count)
    {
        T* old = ptr;
        T* nev = AllocateUninitialized(count);
        MemCpy(nev, old, count);
        Deallocate(old, oldCount);
        return nev;
    }
};

// todo add ScopedPtr
// todo add ScopedFn
// todo SharedPtr in different hpp file 