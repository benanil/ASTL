#pragma once

#include "Memory.hpp"
#include "Algorithms.hpp"

// todo: add priority queue here (min heap) https://github.com/lemire/FastPriorityQueue.js/blob/master/FastPriorityQueue.js
//       use power of two increasing to not use mod operator
//       make queue power of two and use & operator instead of modulo

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
        { return a.m_rear == b.m_rear && a.m_ptr == b.m_ptr; };
		
		friend bool operator != (const Iterator& a, const Iterator& b)
        { return a.m_rear != b.m_rear || a.m_ptr == b.m_ptr; };  

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
        { return a.m_rear == b.m_rear && a.m_ptr == b.m_ptr; };
		
		friend bool operator != (const ConstIterator& a, const ConstIterator & b)
        { return a.m_rear != b.m_rear || a.m_ptr != b.m_ptr; };

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
		ptr[front++] = value;
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
			r &= e;
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
		if (UINT32_MAX-capacity < capacity)
			return UINT32_MAX;
		return capacity + capacity;
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