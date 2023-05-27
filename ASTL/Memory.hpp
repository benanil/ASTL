#pragma once

#include "Common.hpp"
#include <memory.h> // < this is here for memcpy and memset

template<typename T>
struct RemoveRef { typedef T Type; };

template<typename T>
struct RemoveRef<T&> { typedef T Type; };

template<typename T>
struct RemoveRef<T&&> { typedef T Type; };

template<typename T>
struct RemovePtr { typedef T Type; };

template<typename T>
struct RemovePtr<T*> { typedef T Type; };

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

template<typename T>
FINLINE void SwapPrimitive(T& left, T& right)
{
  T temp = Move(left);
  left = Move(right);
  right = Move(temp);
}

template<class T, class U = T>
FINLINE constexpr T Exchange(T& obj, U&& new_value)
{
  T old_value = Move(obj);
  obj = Forward<U>(new_value);
  return old_value;
}

FINLINE void MemSet(void* dst, char val, size_t sizeInBytes)
{
    memset(dst, val, sizeInBytes);
}

FINLINE void MemCpy(void* dst, const void* src, size_t sizeInBytes)
{
    memcpy(dst, src, sizeInBytes);
}

template<typename T> 
struct Allocator 
{
    T* Allocate(int count)
    {
        return new T[count];
    }
    
    T* AllocateUninitialized(int count)
    {
        return new T[count];
    }
    
    void Deallocate(T* ptr, size_t count)
    {
        delete [] ptr;
    }
    
    T* Reallocate(T* ptr, size_t oldCount, size_t count)
    {
        T* old = ptr;
        T* _new = new T[count];
        for (int i = 0; i < Min(oldCount, count); ++i)
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
  T* Allocate(int count)
  {
    return (T*)malloc(sizeof(T) * count);
  }
    
  T* AllocateUninitialized(int count)
  {
    return (T*)calloc(sizeof(T) * count);
  }
    
  void Deallocate(T* ptr, size_t count)
  {
    free(ptr);
  }
    
  T* Reallocate(T* ptr, size_t oldCount, size_t count)
  {
    return (T*)realloc(ptr, sizeof(T) * count);
  }
};
