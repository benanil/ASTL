
/********************************************************************************
*    Purpose:                                                                   *
*        Stack data structure. First in Last Out                                *
*    Author:                                                                    *
*        Anilcan Gulkaya 2023 anilcangulkaya7@gmail.com                         *
********************************************************************************/

#pragma once

#include "Memory.hpp"

AX_NAMESPACE

template<typename T, typename AllocatorT = Allocator<T>>
class Stack
{
public:
	using Iterator = T*;
	using ConstIterator = const T*;
	static const int InitialSize = AllocatorT::InitialSize;
                                                           
public:                                                    
    int size     = 0;                                      
    int capacity = 0;
    T*  ptr      = nullptr;
    AllocatorT allocator{};

public:
	Stack() : size(0), capacity(0), ptr(nullptr) 
	{ }
	
	~Stack()
	{ 
		Clear();
	}

	explicit Stack(int _size) : size(0), capacity(CalculateArrayGrowth(_size)) 
	{
		ptr = allocator.AllocateUninitialized(capacity);
	}
	
	// copy constructor
	Stack(const Stack& other)
	: size(other.size), capacity(other.capacity)
	{
		if (&other != this)
		{
			if (other.capacity > capacity) 
				GrowIfNecessarry(other.capacity - capacity);

			Copy(ptr, other.ptr, other.Size());
			size = other.Size();
		}
	}

	// move constructor 
	Stack(Stack&& other)
	{
		if (ptr != nullptr)
			allocator.Deallocate(ptr, capacity);

		size      = other.size;
		capacity  = other.capacity;
		ptr       = other.ptr;
		other.ptr = nullptr;
		other.size = other.capacity = 0;
	}

	T&       operator[](int index)       { ASSERT(index >= capacity); return ptr[index]; } 
	const T& operator[](int index) const { ASSERT(index >= capacity); return ptr[index]; }
	
	T*       begin()        { return ptr;           }								              
	T*       end()          { return ptr + size;    }							              
	const T* cbegin() const { return ptr;           }					              
	const T* cend()   const { return ptr + size;    }				          
	                   
	const T& Top()     const { return ptr[size - 1]; }                                 
	T&       Top()           { return ptr[size - 1]; }                                        
	bool     Any()     const { return size >  0;     }								  
	bool     Empty()   const { return size == 0;     } 								  
	
	void Clear()  
	{
		if (ptr) 
		{ 
			allocator.Deallocate(ptr, capacity);
			ptr  = nullptr; 
			size = capacity = 0;
		} 
	}

	template<typename ... Args>
	T& Emplace(Args&&... args)
	{
		GrowIfNecessarry();
		new (ptr + size) T(Forward<Args>(args)...);
		return ptr[size++];
	}

	void PushArray(const T* arr, int cnt) 
	{
		GrowIfNecessarry(cnt);
		Copy(ptr + size, arr, cnt);
		size += cnt;
	}

	void Push(const T& value) 
	{
		GrowIfNecessarry();
		ptr[size++] = value;
	}

	T Pop() 
	{
		ASSERT(size != 0);
		return Forward<T>(ptr[--size]); 
	}
	
	bool TryPop(T& out) 
	{
		if (size > 0) 
			out = Forward<T>(ptr[--size]);

		return size > 0;
	}
private:
	void GrowIfNecessarry(int add = 1)
	{
		if (size + add >= capacity)
		{
			int newSize = size + add <= InitialSize ? InitialSize : CalculateArrayGrowth(capacity);
			if (ptr) ptr = allocator.Reallocate(ptr, capacity, newSize);
			else     ptr = allocator.Allocate(newSize);
			capacity    = newSize; 
		}
	}
};

template<typename T, int size>
using FixedStack = Stack<T, StackAllocator<T, size>>;

AX_END_NAMESPACE