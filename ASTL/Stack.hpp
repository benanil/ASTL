#pragma once

#include "Common.hpp"

AX_NAMESPACE

template<typename T, typename AllocatorT = Allocator<T>>
class Stack
{
public:
	using Iterator = T*;
	using ConstIterator = const T*;

    // we don't want to use same initial size for all data types because we want 
	// more initial size for small data types such as byte, short, int but for bigger value types we want less initial size
	// if this initial size is big for your needs, use static array please.[ byte    |  InitialSize ]
	static const int InitialSize = 384 / MIN((int)sizeof(T), 128);    //  1         384
                                                                          //  2         192
public:                                                                   //  4         96
    int size     = 0;                                                     //  8         48
    int capacity = 0;
    T*  ptr      = nullptr;
    AllocatorT allocator{};

public:
	Stack() : size(0), capacity(0), ptr(nullptr) 
	{ }
	
	~Stack()
	{ 
		if (ptr) 
		{ 
			allocator.Deallocate(ptr, capacity);
			ptr  = nullptr; 
			size = capacity = 0;
		} 
	}

	explicit Stack(int _size) : size(0), capacity(CalculateArrayGrowth(_size)) 
	{
		ptr = allocator.AllocateUninitialized(capacity);
	}
	
	Stack(int _size, int count) : size(count), capacity(CalculateArrayGrowth(_size)) 
	{
		ptr = allocator.Allocate(capacity);
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
	bool     Empty() const { return size == 0;     } 								  
	
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
	void GrowIfNecessarry()
	{
		if (size + 1 >= capacity)
		{
			int newSize = MAX(CalculateArrayGrowth(capacity), InitialSize);
			if (ptr) ptr = allocator.Reallocate(ptr, capacity, newSize);
			else     ptr = allocator.Allocate(newSize);
			capacity    = newSize; 
		}
	}
};

// recommended using static array instead of this.
#ifdef false
template<typename T, int MAXSize>
class FixedStack
{
private:
	int size = 0;
public:
	T arr[MAXSize];
	__constexpr FixedStack() { }

	__constexpr void Push(const T& value) { arr[size++] = value; }
	__constexpr T&   Pop()                { return arr[--size];  }
	__constexpr bool Any()   const        { return size > 0;     }
	__constexpr bool Empty() const        { return size == 0;    }
	
	__constexpr T&   operator[](int index) { return arr[index];  }
	__constexpr T    GetLast()             { return arr[size - 1]; }
	__constexpr T    GetFirst()            { return arr[0]; }
	__constexpr void Clear()               { memset(arr, 0, MAXSize * sizeof(T)); size = 0; }

	__constexpr T* begin()              { return arr;        }
	__constexpr T* end()                { return arr + size; }
	__constexpr const T* cbegin() const { return arr;        }
	__constexpr const T* cend()   const { return arr + size; }
};
#endif

AX_END_NAMESPACE