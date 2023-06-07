#pragma once

#include "Memory.hpp"
#include "Algorithms.hpp"

// todo: add priority queue here (min heap)

template<typename ValueT,
         typename AllocatorT = Allocator<ValueT>>
class Queue
{
public:
	struct Iterator 
	{
		Iterator(ValueT* ptr, int capacity, int rear) 
		: m_ptr(ptr), m_capacity(capacity), m_rear(rear) {}
		~Iterator() {}

		ValueT&   operator * ()    { return m_ptr[m_rear]; }
		ValueT*   operator ->()    { return m_ptr + m_rear; }
		Iterator& operator ++()    { m_rear = (m_rear + 1) % m_capacity; return *this; }
		Iterator  operator ++(int) { Iterator tmp = *this; ++(*this); return tmp; }
		
		friend bool operator == (const Iterator& a, const Iterator& b)
        { return a.m_rear == b.m_rear && a.m_ptr == b.m_ptr; };
		
		friend bool operator != (const Iterator& a, const Iterator& b)
        { return a.m_rear != b.m_rear || a.m_ptr == b.m_ptr; };  

	private:
		ValueT* m_ptr;
		int m_capacity, m_rear;
	};

	struct ConstIterator 
	{
		ConstIterator (const ValueT* ptr, int capacity, int rear) 
		: m_ptr(ptr), m_capacity(capacity), m_rear(rear) {}
		~ConstIterator() {}

		const ValueT&  operator * ()  const { return m_ptr[m_rear]; }
		const ValueT*  operator ->()  const { return m_ptr + m_rear; }
		ConstIterator& operator ++()        { m_rear = (m_rear + 1) % m_capacity; return *this; }  
		ConstIterator  operator ++(int)     { ConstIterator tmp = *this; ++(*this); return tmp; }
		
		friend bool operator == (const ConstIterator& a, const ConstIterator & b)
        { return a.m_rear == b.m_rear && a.m_ptr == b.m_ptr; };
		
		friend bool operator != (const ConstIterator& a, const ConstIterator & b)
        { return a.m_rear != b.m_rear || a.m_ptr != b.m_ptr; };

	private:
		const ValueT* m_ptr;
		int m_capacity, m_rear;
	};

	using iterator = Iterator; // stl compatible
	using const_iterator = ConstIterator;

private:
	// we don't want to use same initial size for all data types because we want 
	// more initial size for small data types such as byte, short, int but for bigger value types we want less initial size
	// if this initial size is big for your needs, use static array please.[ byte    |  InitialSize ]
	static constexpr int InitialSize = 384 / Min((int)sizeof(ValueT), 128);//  1         384
							                                               //  2         192
	ValueT* ptr  = nullptr;                                                //  4         96
	int capacity = InitialSize;                                            //  8         48         : minimum initial size is 8
	int front    = 0;
	int rear     = 0;
	AllocatorT allocator{};
public:
	Queue() : front(0), rear(0), capacity(InitialSize)
	{
		ptr = allocator.AllocateUninitialized(capacity);
	}
	
	~Queue()
	{
		if (ptr != nullptr)
		{
			allocator.Deallocate(ptr, capacity);
			ptr      = nullptr;
			capacity = front = rear = 0;
		}
	}

	explicit Queue(int _capacity) : capacity(_capacity + 1), front(0), rear(0)
	{
		ptr = allocator.AllocateUninitialized(capacity);
	}

	Queue(int _capacity, int count) : capacity(_capacity + 1), front(count), rear(0)
	{
		ptr = allocator.Allocate(capacity);
	}

	Queue(int _size, const ValueT& val) 
	: capacity(_size + (_size / 2)), front(0), rear(0)
	{
		ptr = allocator.AllocateUninitialized(capacity);
		
		for (; front < _size; front++)
		{
			ptr[front] = val;
		}
	}

	Queue(ValueT* p, int size)
	: capacity(size + (size / 2)), front(size), rear(0)
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
		    int newSize = other.Size();
			GrowIfNecessary(newSize);
			const ValueT* otherPtr  = other.ptr;
			const int otherCapacity = other.capacity;
			int otherRear = other.rear;

			for (int i = 0; i < newSize; ++i)
			{
				ptr[i].~ValueT();
				ptr[i] = otherPtr[otherRear];
				otherRear = (otherRear + 1) % otherCapacity;
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
		ptr[front].~ValueT();
		new (ptr + front) ValueT(Forward<Args>(args)...);
		int idx = front;
		front   = IncrementIndex(front);
		return ptr[idx];
	}

	void Enqueue(const ValueT& value)
	{
		GrowIfNecessary(1);
		ptr[front++] = value;
	}

	void Enqueue(const ValueT* begin, int count)
	{
		GrowIfNecessary(count);
		int f = front;
		for (int i = 0; i < count; i++)
		{
			ptr[f] = begin[i];
			f      = IncrementIndex(f);
		}
		front += count;
	}

	void Enqueue(const ValueT* begin, const ValueT* end)
	{
		Enqueue(begin, PointerDistance(begin, end));
	}

	// returns true if size is enough
	bool TryDequeue(ValueT* result, int count)
	{
		if (Size() + count > capacity)
		{
			return false;
		}
		int r = rear;
		for (int i = 0; i < count; i++)
		{
			result[i] = Forward<ValueT>(ptr[r]);
			r         = IncrementIndex(r);
		}
		rear = r;
		return true;
	}

	void Dequeue(ValueT* result, int count)
	{
		ASSERT(Size() + count <= capacity);
		int r = rear;
		for (int i = 0; i < count; i++)
		{
			result[i] = Forward<ValueT>(ptr[r]);
			r         = IncrementIndex(r);
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

	int Size() const 
	{
		if (front < rear)
		{
			return front + (capacity - rear);
		}
		return front - rear; 
	}

	// todo: shrink to fit

private:

	int IncrementIndex(int x) const
	{
		return (x + 1) % capacity;
	}

	int CalculateArrayGrowth() const
	{
		if (INT32_MAX-capacity < capacity)
			return INT32_MAX;
		return capacity + capacity;
	}

	void GrowIfNecessary(int _size)
	{
		int size = Size();
		if (AX_LIKELY(size + _size < capacity))
		{
			return; // no need to grow
		}

		const int newCapacity = Max(CalculateArrayGrowth(), 8);
		ptr                   = allocator.Reallocate(ptr, capacity, newCapacity);
		
		// unify front and rear, if they are seperate.
		if (front < rear)
		{
			for (int i = 0; i < front; i++)
			{
				ptr[capacity++] = Forward<ValueT>(ptr[i]);
			}
		}
		
		// move everything to the right
		for (int i = 0; i < size; ++i)
		{
			ptr[newCapacity - 1 - i] = Forward<ValueT>(ptr[capacity - 2 - i]);
		}

		rear     = newCapacity - size; 
		front    = 0;
		capacity = newCapacity;
	}
};