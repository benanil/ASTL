#pragma once

#include "Memory.hpp"
#include "Algorithms.hpp"

// I use power of two increasing to not use modulo operator
// because and operator is faster

template<typename ValueT,
         typename AllocatorT = Allocator<ValueT>>
class Queue
{
public:
	struct Iterator 
	{
		Iterator(ValueT* ptr, uint capacity, uint rear) 
		: m_ptr(ptr), m_mod(capacity-1), m_rear(rear) {}
		~Iterator() {}

		ValueT&   operator * ()    { return m_ptr[m_rear]; }
		ValueT*   operator ->()    { return m_ptr + m_rear; }
		Iterator& operator ++()    { m_rear = (m_rear + 1) & m_mod; return *this; }
		Iterator  operator ++(int) { Iterator tmp = *this; ++(*this); return tmp; }
		
		friend bool operator == (const Iterator& a, const Iterator& b)
        { return (a.m_rear == b.m_rear) & a.m_ptr == b.m_ptr; };
		
		friend bool operator != (const Iterator& a, const Iterator& b)
        { return (a.m_rear != b.m_rear) | (a.m_ptr == b.m_ptr); };  

	private:
		ValueT* m_ptr;
		uint m_mod, m_rear;
	};

	struct ConstIterator 
	{
		ConstIterator (const ValueT* ptr, uint capacity, uint rear) 
		: m_ptr(ptr), m_mod(capacity-1), m_rear(rear) {}
		~ConstIterator() {}

		const ValueT&  operator * ()  const { return m_ptr[m_rear]; }
		const ValueT*  operator ->()  const { return m_ptr + m_rear; }
		ConstIterator& operator ++()        { m_rear = (m_rear + 1) & m_mod; return *this; }  
		ConstIterator  operator ++(int)     { ConstIterator tmp = *this; ++(*this); return tmp; }
		
		friend bool operator == (const ConstIterator& a, const ConstIterator & b)
        { return (a.m_rear == b.m_rear) & (a.m_ptr == b.m_ptr); };
		
		friend bool operator != (const ConstIterator& a, const ConstIterator & b)
        { return (a.m_rear != b.m_rear) | (a.m_ptr != b.m_ptr); };

	private:
		const ValueT* m_ptr;
		uint m_mod, m_rear;
	};

	using iterator = Iterator; // stl compatible
	using const_iterator = ConstIterator;

private:
	// we don't want to use same initial size for all data types because we want 
	// more initial size for small data types such as byte, short, int but for bigger value types we want less initial size
	static constexpr int InitialSize = NextPowerOf2(512 / Min((int)sizeof(ValueT), 128));
							                                               
	ValueT* ptr  = nullptr;                                                
	uint capacity = InitialSize;                                            
	uint front    = 0;
	uint rear     = 0;
	AllocatorT allocator{};
public:
	Queue() : ptr(nullptr), capacity(0), front(0), rear(0)
	{ }
	
	~Queue()
	{
		if (ptr != nullptr)
		{
			allocator.Deallocate(ptr, capacity);
			ptr      = nullptr;
			capacity = front = rear = 0;
		}
	}

	explicit Queue(uint _capacity) : capacity(NextPowerOf2(_capacity + 1)), front(0), rear(0)
	{
		ptr = allocator.AllocateUninitialized((int)capacity);
	}

	Queue(uint _capacity, uint count) : capacity(NextPowerOf2(_capacity + 1)), front(count), rear(0)
	{
		ptr = allocator.Allocate((int)capacity);
	}

	Queue(uint _size, const ValueT& val) 
	: capacity(NextPowerOf2(_size + (_size / 2))), front(0), rear(0)
	{
		ptr = allocator.AllocateUninitialized(capacity);
		
		for (; front < _size; front++)
		{
			ptr[front] = val;
		}
	}

	Queue(ValueT* p, uint size)
	: capacity(NextPowerOf2(size + (size / 2))), front(size), rear(0)
	{
		ptr = allocator.AllocateUninitialized(capacity);
		Copy(ptr, p, size);
	}

	// copy constructor
	Queue(const Queue& other) 
	: capacity(other.capacity), front(other.front), rear(other.rear)
	{
		Copy(ptr, other.ptr, other.capacity);
	}

	Queue& operator = (Queue const& other)
	{
		if (&other != this)
		{
			uint newSize = other.Size();
			GrowIfNecessary(newSize);
			const ValueT* otherPtr  = other.ptr;
			const uint otherCapacity = other.capacity;
			uint otherRear = other.rear;
			uint e = other.capacity-1;

			for (int i = 0; i < newSize; ++i)
			{
				ptr[i].~ValueT();
				ptr[i] = otherPtr[otherRear];
				otherRear = (otherRear + 1) & e;
			}
		}
		return *this;
	}
	
	Queue& operator = (Queue&& other) noexcept
	{
		if (&other != this)
		{
		    if (ptr != nullptr)
		   		allocator.Deallocate(ptr, capacity);

			ptr            = other.ptr;
			capacity       = other.capacity;
			front          = other.front;
			rear           = other.rear;
			other.capacity = other.front = other.rear = 0;
			other.ptr      = nullptr;
		}
		return *this;
	}

	bool operator == (const Queue& other) const
	{
		if (&other == this) return true;
		if (other.Size() != Size()) return false;
		
		for (ConstIterator ob = other.cbegin(), tb = cbegin();
           ob != other.cend(); ob++, tb++)
		{
			if (*ob != *tb)
				return false;
		}
		return true;
	}

	bool operator != (const Queue& other) const
	{
		return !(*this == other);
	}

	void Clear()
	{
		if (capacity != InitialSize)
		{
			ptr      = allocator.Reallocate(ptr, InitialSize, capacity);
			capacity =  front = rear = 0;
		}
	}

	Iterator begin()              { return { ptr, capacity, rear  }; }
	Iterator end()                { return { ptr, capacity, front }; }
	ConstIterator cbegin() const  { return { ptr, capacity, rear  }; }
	ConstIterator cend()   const  { return { ptr, capacity, front }; }

	bool Any()     const { return Size() > 0;  }
	bool IsEmpty() const { return Size() == 0; }

	template<typename... Args>
	ValueT& Emplace(Args&&... args)
	{
		GrowIfNecessary(1);
		new (ptr + front) ValueT(Forward<Args>(args)...);
		uint idx = front;
		front   = IncrementIndex(front);
		return ptr[idx];
	}

	void Enqueue(const ValueT& value)
	{
		GrowIfNecessary(1);
		ptr[front] = value;
		front   = IncrementIndex(front);
	}

	void Enqueue(const ValueT* begin, uint count)
	{
		GrowIfNecessary(count);
		uint f = front, e = capacity-1;
		for (uint i = 0; i < count; i++)
		{
			ptr[f++] = begin[i];
			f &= e;
		}
		front += count;
	}

	void Enqueue(const ValueT* begin, const ValueT* end)
	{
		Enqueue(begin, PointerDistance(begin, end));
	}

	// returns true if size is enough
	bool TryDequeue(ValueT* result, uint count)
	{
		if (Size() + count > capacity)
		{
			return false;
		}
		uint r = rear, e = capacity-1;

		for (uint i = 0; i < count; i++)
		{
			result[i] = Forward<ValueT>(ptr[r++]);
			r &= e;
		}
		rear = r;
		return true;
	}

	void Dequeue(ValueT* result, int count)
	{
		ASSERT(Size() + count <= capacity);
		uint r = rear, e = capacity-1;
		for (uint i = 0; i < count; i++)
		{
			result[i] = Forward<ValueT>(ptr[r++]);
			r &= e; // better than modulo 
		}
		rear = r;
	}

	ValueT Dequeue()
	{	
		ASSERT(rear != front);
		ValueT val = Forward<ValueT>(ptr[rear]);
		rear = IncrementIndex(rear);
		return val; 
	}

	bool TryDequeue(ValueT& out)
	{
		if (AX_UNLIKELY(rear == front)) 
			return false;
		out = Forward<ValueT>(ptr[rear]);
		rear = IncrementIndex(rear);
		return true;
	}

	uint Size() const 
	{
		if (front < rear)
		{
			return front + (capacity - rear);
		}
		return front - rear; 
	}

	// todo: shrink to fit

private:

	uint IncrementIndex(uint x) const
	{
		return (x + 1) & (capacity-1);
	}

	uint CalculateQueueGrowth() const
	{
		ASSERT(capacity != (1 << 31)); // max size
		return capacity << 1;
	}

	void GrowIfNecessary(uint _size)
	{
		uint size = Size();
		if (AX_LIKELY(size + _size < capacity))
		{
			return; // no need to grow
		}

		const uint newCapacity = Max(CalculateQueueGrowth(), InitialSize);
		if (ptr)  ptr = allocator.Reallocate(ptr, capacity, newCapacity);
		else      ptr = allocator.Allocate(newCapacity);

		// unify front and rear, if they are seperate.
		if (front < rear)
		{
			for (uint i = 0; i < front; i++)
			{
				ptr[capacity++] = Forward<ValueT>(ptr[i]);
			}
		}
		
		// move everything to the right
		for (uint i = 0; i < size; ++i)
		{
			ptr[newCapacity - 1 - i] = Forward<ValueT>(ptr[capacity - 2 - i]);
		}

		rear     = newCapacity - size; 
		front    = 0;
		capacity = newCapacity;
	}
};


/////   PriorityQueue /////


enum Compare
{
	Compare_Less    = 0,
	Compare_Equal   = 1,
	Compare_Greater = 2
};

template<typename T, Compare op>
struct Comparer
{
	static bool Compare(const T& a, const T& b)
	{
		if constexpr (op == Compare_Less)    return a < b;
		if constexpr (op == Compare_Equal)   return a == b;
		if constexpr (op == Compare_Greater) return a > b;
		else static_assert(true && "undefined operator");
	}
};

template<typename T, 
         typename ComparerT comp = Comparer<T, Compare_Less>,
         typename AllocatorT     = Allocator<T>>
class PriorityQueue
{
    // we don't want to use same initial size for all data types because we want 
	// more initial size for small data types such as byte, short, int but for bigger value types we want less initial size
	static constexpr int InitialSize = 512 / Min((int)sizeof(T), 128);
public:
	T*      array    = nullptr;
	int     size     = 0;
	int     capacity = 0;
	AllocatorT allocator;

	PriorityQueue() {}
	
	~PriorityQueue()
	{
		if (array)
		{
			allocator.Deallocate(array, capacity);
			array = nullptr;
			size  = capacity = 0;
		}
	}

	explicit Stack(int _size) : size(0), capacity(CalculateArrayGrowth(_size)) 
	{
		array = allocator.AllocateUninitialized(capacity);
	}
	
	PriorityQueue& operator=(const PriorityQueue& other) 	
	{
	    if (this != &other) 
		{
		    if (other.size > capacity) 
				GrowIfNecessarry(other.size - capacity);
			Copy(array, other.array, other.size);
	        size = other.size;
	    }
	    return *this;
	}

	// move constructor 
	PriorityQueue(PriorityQueue&& other)
	{
		if (ptr != nullptr)
			allocator.Deallocate(ptr, capacity);

		size      = other.size;
		capacity  = other.capacity;
		ptr       = other.ptr;
		other.ptr = nullptr;
		other.size = other.capacity = 0;
	}

	void GrowIfNecessarry(int adition)
	{
		if (size + adition >= capacity)
		{
			int newCapacity = Max(CalculateArrayGrowth(size + adition), InitialSize);
			if (array)
				array = allocator.Reallocate(array, capacity, newCapacity);
			else
				array = allocator.Allocate(newCapacity);
			capacity = newCapacity;
		}
	}

	int InsertFixup(int i, const T& myval)
	{
		while (i > 0) 
		{
			int p  = (i - 1) >> 1;
			if (!ComparerT::Compare(myval, array[p]))
				break;
            
			array[i] = (T&&)array[p];
			i = p;
		}
		return i;
	}

	void Push(const T& myval) 
	{
		GrowIfNecessarry(1);
		int i = size++;
		i = InsertFixup(i, myval);
		array[i] = (T&&)myval;
	}

	template<typename ... Args>
	void Emplace(Args&&... args) 
	{
		GrowIfNecessarry(1);
		T myval(Forward<Args>(args)...);

		int i = size++;
		i = InsertFixup(i, myval);
		
		array[i] = (T&&)myval;
	}

	bool Remove(const T& myval) 
	{
		for (int i = 0; i < this->size; i++) 
		{
			if (this->array[i] == myval) 
			{
				RemoveAt(i);
				return true;
			}
		}
		return false;
	}

	T RemoveAt(int index) 
	{
		ASSERT(!(index > this->size - 1 || index < 0));
		PercolateUp(index, true);
		return this->Pop();
	}

	void PercolateUp(int i, bool force = false)
	{
		T myval = (T&&)array[i];
		while (i > 0) 
		{
			int p  = (i - 1) >> 1;
			
			if (!force && !ComparerT::Compare(myval, array[p]))
 				break;
			
			array[i] = (T&&)array[p];
			i = p;
    	}
		array[i] = (T&&)myval;
	}

	void PercolateDown(int i) 
	{
		int size  = this->size;
		int hsize = this->size >> 1;
		T ai = (T&&)array[i];

		while (i < hsize)
		{
			int l = (i << 1) + 1;
			int r = l + 1;
			
			if (r < size)
			{
				if (ComparerT::Compare(array[r], array[l]))
				{
					l = r;
				}
			}
			
			if (!ComparerT::Compare(array[l], ai))
				break;
			
			array[i] = (T&&)array[l];
			i = l;
		}
		array[i] = (T&&)ai;
	}

	const T* begin() const { return this->array; }
	const T* end()   const { return this->array + size; }
	T* begin() { return this->array; }
	T* end() { return this->array + size; }

	const T& Peek() const 
	{
		return array[0];
	}

	T Pop() 
	{
		ASSERT(size != 0);
		T ans = (T&&)array[0];
		if (size-- > 1) 
		{
			array[0] = array[size];
			PercolateDown(0);
		} 
		return ans;
	}

	bool IsEmpty() const { return size == 0; }
	bool     Any() const { return size >  0; }					
	
	void Clear()  
	{
		if (capacity != InitialSize)
		{
			if (array)
				allocator.Deallocate(array, capacity);
			size = capacity = 0; 
		}
	}
};