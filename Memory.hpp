
/****************************************************************
*    Purpose: General purpose memory functions                  *
*    Author : Anilcan Gulkaya 2023 anilcangulkaya7@gmail.com    *
****************************************************************/

// T&& Move(T&& obj);
// T&& Forward(T& obj);
// T Exchange(T& obj, U&& new_value);
// uint64_t AlignAddress(uint64_t addr, uint64_t align)
// T* AlignPointer(T* ptr, uint64_t align)
// void* AllocAligned(uint64_t bytes, uint64_t align)
// void FreeAligned(void* pMem)
// void MemSet(void* dst, unsigned char val, uint64_t  sizeInBytes)
// void MemCpy(void* dst, const void* src, uint64_t  sizeInBytes)
// struct Allocator;
// struct MallocAllocator;
// struct FixedSizeGrowableAllocator;
// struct StackAllocator;

#pragma once

#include "Common.hpp"

AX_NAMESPACE 

template<typename T> struct RemoveRef      { typedef T Type; };
template<typename T> struct RemoveRef<T&>  { typedef T Type; };
template<typename T> struct RemoveRef<T&&> { typedef T Type; };
template<typename T> struct RemovePtr      { typedef T Type; };
template<typename T> struct RemovePtr<T*>  { typedef T Type; };

template<typename T>
__forceinline typename RemoveRef<T>::Type&& Move(T&& obj)
{
  typedef typename RemoveRef<T>::Type CastType;
  return (CastType&&)obj;
}

template<typename T>
__forceinline T&& Forward(typename RemoveRef<T>::Type& obj) { return (T&&)obj; }

template<typename T>
__forceinline T&& Forward(typename RemoveRef<T>::Type&& obj) { return (T&&)obj; }

template<class T, class U = T>
__forceinline __constexpr T Exchange(T& obj, U&& new_value)
{
  T old_value = Move(obj);
  obj = Forward<U>(new_value);
  return old_value;
}

// Aligned memory codes taken from Game Engine Architecture book by Jason Gregory

// Shift the given address upwards if/as necessary to// ensure it is aligned to the given number of bytes.
inline uint64_t AlignAddress(uint64_t addr, uint64_t align)
{
    const uint64_t  mask = align - 1;
    ASSERT((align & mask) == 0); // pwr of 2
    return (addr + mask) & ~mask;
}

// Shift the given pointer upwards if/as necessary to// ensure it is aligned to the given number of bytes.
template<typename T>
inline T* AlignPointer(T* ptr, uint64_t align)
{
    const uint64_t  addr = (uint64_t )ptr;
    const uint64_t  addrAligned = AlignAddress(addr, align);
    return (T*)addrAligned;
}

// Aligned allocation function. IMPORTANT: 'align'// must be a power of 2 (typically 4, 8 or 16).
inline void* AllocAligned(uint64_t bytes, uint64_t align)
{
    // Allocate 'align' more bytes than we need.
    uint64_t  actualBytes = bytes + align;
    // Allocate unaligned block.
    uint8* pRawMem = new uint8[actualBytes];
    // Align the block. If no alignment occurred,// shift it up the full 'align' bytes so we// always have room to store the shift.
    uint8* pAlignedMem = AlignPointer(pRawMem, align);
    
    if (pAlignedMem == pRawMem)
        pAlignedMem += align;
    // Determine the shift, and store it.// (This works for up to 256-byte alignment.)
    uint8 shift = (uint8)(pAlignedMem - pRawMem);
    ASSERT(shift > 0 && shift <= 256);
    pAlignedMem[-1] = (uint8)(shift & 0xFF);
    return pAlignedMem;
}

inline void FreeAligned(void* pMem)
{
    ASSERT(pMem);
    // Convert to U8 pointer.
    uint8* pAlignedMem = (uint8*)pMem;
    // Extract the shift.
    uint64_t  shift = pAlignedMem[-1];
    if (shift == 0)
        shift = 256;
    // Back up to the actual allocated address,
    uint8* pRawMem = pAlignedMem - shift;
    delete[] pRawMem;
}

// if you want memmove it is here with simd version: https://hackmd.io/@AndybnA/0410

inline void MemSetAligned64(void* RESTRICT dst, unsigned char val, uint64_t sizeInBytes)
{
    // we want an offset because the while loop below iterates over 4 uint64_t at a time
    const uint64_t * end = (uint64_t *)((char*)dst + (sizeInBytes >> 3));
    uint64_t * dp = (uint64_t *)dst;
#ifdef _MSC_VER
    uint64_t   d4;
    __stosb((unsigned char*)&d4, val, 8);
#else
    uint64_t   d4 = val; d4 |= d4 << 8ull; d4 |= d4 << 16ull; d4 |= d4 << 32ull;
#endif

    while (dp < end)
    {
        dp[0] = dp[1] = dp[2] = dp[3] = d4;
        dp += 4;
    }
   
    switch (sizeInBytes & 7)
    {
        case 7: *dp++ = d4;
        case 6: *dp++ = d4;
        case 5: *dp++ = d4;
        case 4: *dp++ = d4;
        case 3: *dp++ = d4;
        case 2: *dp++ = d4;
        case 1: *dp   = d4;
    };
}

inline void MemSetAligned32(uint32_t* RESTRICT dst, unsigned char val, uint64_t sizeInBytes)
{
    const uint32_t* end = (uint32_t*)((char*)dst + (sizeInBytes >> 2));
    uint32_t* dp = (uint32_t*)dst;
#ifdef _MSC_VER
    uint32_t  d4;
    __stosb((unsigned char*)&d4, val, 4);
#else
    uint32_t  d4 = val; d4 |= d4 << 8ull; d4 |= d4 << 16ull;
#endif
    while (dp < end)
    {
        dp[0] = dp[1] = dp[2] = dp[3] = d4;
        dp += 4;
    }
    
    switch (sizeInBytes & 3)
    {
        case 3: *dp++ = d4;
        case 2: *dp++ = d4;
        case 1: *dp   = d4;
    };
}

// use size for struct/class types such as Vector3 and Matrix4, 
// leave as zero size constant for big arrays or unknown size arrays
template<int alignment = 0>
inline void MemSet(void* dst, unsigned char val, uint64_t sizeInBytes)
{
    if_constexpr (alignment == 8)
        MemSetAligned64(dst, val, sizeInBytes);
    else if (alignment == 4)
        MemSetAligned32((uint32_t*)dst, val, sizeInBytes);
    else
    {
        uint64_t  uptr = (uint64_t )dst;
        if (!(uptr & 7) && uptr >= 8) MemSetAligned64(dst, val, sizeInBytes);
        else if (!(uptr & 3) && uptr >= 4)  MemSetAligned32((uint32_t*)dst, val, sizeInBytes);
        else
        {
            unsigned char* dp = (unsigned char*)dst;
            while (sizeInBytes--)
                *dp++ = val;
        }
    }
}

inline void MemCpyAligned64(void* dst, const void* RESTRICT src, uint64_t sizeInBytes)
{
    uint64_t *       dp  = (uint64_t *)dst;
    const uint64_t * sp  = (const uint64_t *)src;
    const uint64_t * end = (const uint64_t *)((char*)src) + (sizeInBytes >> 3);
        
    while (sp < end)
    {
        dp[0] = sp[0];
        dp[1] = sp[1];
        dp[2] = sp[2];
        dp[3] = sp[3];
        dp += 4, sp += 4;
    }

    SmallMemCpy(dp, sp, sizeInBytes & 7);
}

inline void MemCpyAligned32(uint32_t* dst, const uint32_t* RESTRICT src, uint64_t sizeInBytes)
{
    uint32_t*       dp  = (uint32_t*)dst;
    const uint32_t* sp  = (const uint32_t*)src;
    const uint32_t* end = (const uint32_t*)((char*)src) + (sizeInBytes >> 2);
        
    while (sp < end)
    {
        dp[0] = sp[0];
        dp[1] = sp[1];
        dp[2] = sp[2];
        dp[3] = sp[3];
        dp += 4, sp += 4;
    }
    
    switch (sizeInBytes & 3)
    {
        case 3: *dp++ = *sp++;
        case 2: *dp++ = *sp++;
        case 1: *dp++ = *sp++;
    };
}

// use size for structs and classes such as Vector3 and Matrix4,
// and use MemCpy for big arrays or unknown size arrays
template<int alignment = 0>
inline void MemCpy(void* dst, const void* RESTRICT src, uint64_t sizeInBytes)
{
    if_constexpr (alignment == 8)
        MemCpyAligned64(dst, src, sizeInBytes);
    else if (alignment == 4)
        MemCpyAligned32((uint32_t*)dst, (const uint32_t*)src, sizeInBytes);
    else
    {
        const char* cend = (char*)((char*)src + sizeInBytes);
        const char* scp = (const char*)src;
        char* dcp = (char*)dst;
        
        while (scp < cend) *dcp++ = *scp++;
    }
}

// some kind of object allocator
template<typename T> 
struct Allocator 
{
    static const bool IsPod = false;
    // we don't want to use same initial size for all data types because we want 
    // more initial size for small data types such as byte, short, int but for bigger value types we want less initial size
    static const int InitialSize = 512 / MIN((int)sizeof(T), 128);

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
        for (int i = 0; i < MIN(count, oldCount); ++i)
        {
            _new[i] = (T&&)old[i];
        }
        delete[] old;
        return _new;
    }
};

template<typename T>
struct MallocAllocator
{
    static const bool IsPod = true;
    static const int InitialSize = 512 / MIN((int)sizeof(T), 128);

    T* Allocate(int count) const 
    {
        return new T[count]{};
    }

    T* AllocateUninitialized(int count) const 
    {
        return new T[count];
    }
      
    void Deallocate(T* ptr, int count) const 
    {
        delete[] ptr;
    }
      
    T* Reallocate(T* ptr, int oldCount, int count) const
    {
        T* old = ptr;
        T* nev = AllocateUninitialized(count);
        MemCpy<alignof(T)>(nev, old, MIN(oldCount, count) * sizeof(T));
        Deallocate(old, oldCount);
        return nev;
    }
};

template<typename T>
struct FixedSizeGrowableAllocator
{
    static const bool IsPOD = false;
    
    static __constexpr int InitialSize()
    {
        return NextPowerOf2(512 / MIN((int)sizeof(T), 128));
    } 

    struct Fragment
    {
        Fragment* next;
        T* ptr;
        int64_t   size; // used like a index until we fill the fragment
    };

    int currentCapacity = 0;
    Fragment* base = nullptr;
    Fragment* current = nullptr;

    FixedSizeGrowableAllocator(int initialSize)
    {
        // WARNING initial size must be power of two
        ASSERT((initialSize & (initialSize - 1)) == 0);
        currentCapacity = initialSize;
        base = new Fragment;
        current = base;
        base->next = nullptr;
        base->ptr = new T[initialSize];
        base->size = 0;
    }

    FixedSizeGrowableAllocator() : FixedSizeGrowableAllocator(InitialSize())
    { }

    ~FixedSizeGrowableAllocator()
    {
        if (!base) return;
        while (base)
        {
            delete[] base->ptr;
            Fragment* oldBase = base;
            base = base->next;
            delete oldBase;
        }
    }
    
    // move constructor
    FixedSizeGrowableAllocator(FixedSizeGrowableAllocator&& other)
    {
        currentCapacity = other.currentCapacity;
        base = other.base;
        current = other.current;
        other.currentCapacity = 0;
        other.base = other.current = nullptr;
    }

    // copy constructor.
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
        base->ptr = new T[currentCapacity];
        base->size = totalSize;
        current = base;

        // copy other's memory to ours
        T* curr = base->ptr;
        start = other.base;
        while (start)
        {
            Copy(curr, start->ptr, start->size);
            curr += start->size;
            start = start->next;
        }
    }

    void* TakeOwnership()
    {
        void* result = base;
        base = nullptr;
        return result;
    }

    void CheckFixGrow(int count)
    {
        if (current->size + count >= currentCapacity)
        {
            while (currentCapacity < current->size + count)
                currentCapacity <<= 1;

            current->next = new Fragment{};
            current = current->next;
            current->next = nullptr;
            current->ptr = new T[currentCapacity];
            current->size = 0;
        }
    }

    T* Allocate(int count)
    {
        CheckFixGrow(count);
        T* ptr = current->ptr + current->size;
        current->size += count;
        T def{};
        for (int i = 0; i < count; i++)
            ptr[i] = def;
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

// stack allocation for data structures
// this is a lot faster than heap allocations such as malloc
// fixed size, recommended to use for small allocations
// for safety if we grow too much this will use heap memory
template<typename T, int capacity, bool isPod = false>
struct StackAllocator
{
    static const bool IsPod = isPod;
    static const int InitialSize = capacity;

    T arr[capacity]{};

    T* Allocate(int count)
    {
        if (count <= capacity) return arr;
        return new T[count]{};
    }

    T* AllocateUninitialized(int count) 
    {
        if (count <= capacity) return arr;
        return new T[count];
    }
      
    void Deallocate(T* ptr, int count) const 
    {
        // is stack allocated?
        if (ptr >= arr && ptr <= arr + capacity)
        {
            for (int i = 0; i < count; i++)
                ptr[i].~T();
        }
        else
        {
            delete[] ptr;
        }
    }
      
    T* Reallocate(T* ptr, int oldCount, int count) const
    {
        bool stackAllocated = ptr >= arr && ptr <= arr + capacity;
        T* old = ptr, *_new;

        if (stackAllocated && count > capacity) 
        {
            // warning stack size exceeded consider using greater capacity or use heap memory
            _new = new T[count];
            for (int i = 0; i < oldCount; i++)
                _new[i] = (T&&)old[i];
            return _new;
        }
        else if (!stackAllocated)
        {
            _new = new T[count];
            if_constexpr (IsPod)
            {
                MemCpy<alignof(T)>(_new, old, MIN(oldCount, count) * sizeof(T));
            }
            else
            {
                for (int i = 0; i < MIN(count, oldCount); ++i)
                {
                    _new[i] = (T&&)old[i];
                }
            }
            delete[] old;
            return _new;
        }
        // if stack allocated and size is not greater than capacity code does nothing
        return old;
    }
};

template<typename T>
struct ScopedPtr
{
    T* ptr;
    
    ScopedPtr(T* x) : ptr(x) {}
    ~ScopedPtr() { delete[] ptr; ptr = nullptr; }
    
    T* operator->() { return ptr; }
    operator T*() const { return ptr; }

          T& operator[](int n)       { return ptr[n]; }
    const T& operator[](int n) const { return ptr[n]; }

};

// todo add ScopedFn
// todo SharedPtr in different hpp file 

AX_END_NAMESPACE 