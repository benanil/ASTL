#pragma once

#include "Memory.hpp"
#include "Algorithms.hpp"

AX_NAMESPACE

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
		
		friend bool operator == (const ConstIterator& a, const ConstIterator& b)
        { return (a.m_rear == b.m_rear) & (a.m_ptr == b.m_ptr); };
		
		friend bool operator != (const ConstIterator& a, const ConstIterator& b)
        { return (a.m_rear != b.m_rear) | (a.m_ptr != b.m_ptr); };

	private:
		const ValueT* m_ptr;
		uint m_mod, m_rear;
	};

	using iterator = Iterator; // stl compatible
	using const_iterator = ConstIterator;
	
	static __constexpr int InitialSize()
	{
		return IsPowerOfTwo(AllocatorT::InitialSize) == false ?
           NextPowerOf2(AllocatorT::InitialSize) : AllocatorT::InitialSize;
	} 

private:
							                                               
	ValueT* ptr   = nullptr;                                                
	uint capacity = 0;                                            
	uint front    = 0;
	uint rear     = 0;
	uint size     = 0;	
	AllocatorT allocator{};
public:
	Queue() : ptr(nullptr), capacity(0), front(0), rear(0), size(0)
	{ }
	
	~Queue()
	{
		if (ptr != nullptr)
		{
			allocator.Deallocate(ptr, capacity);
			ptr      = nullptr;
			capacity = front = rear = size = 0;
		}
	}

	explicit Queue(int _capacity) : capacity(NextPowerOf2(_capacity + 1)), front(0), rear(0), size(0)
	{
		ptr = allocator.AllocateUninitialized((int)capacity);
	}

	Queue(int _capacity, uint count) : capacity((uint)NextPowerOf2(_capacity + 1)), front(count), rear(0), size(0)
	{
		ptr = allocator.Allocate((int)capacity);
	}

	Queue(int _size, const ValueT& val) 
	: capacity((uint)NextPowerOf2(_size)), front(0), rear(0), size(_size)
	{
		ptr = allocator.AllocateUninitialized(capacity);
		
		for (; front < (uint)_size; front++)
		{
			ptr[front] = val;
		}
	}

	Queue(ValueT* p, int _size)
	: capacity((uint)NextPowerOf2(_size)), front((uint)size), rear(0), size(_size)
	{
		ptr = allocator.AllocateUninitialized(capacity);
		Copy(ptr, p, size);
	}

	// copy constructor
	Queue(const Queue& other) 
	: capacity(other.capacity), front(other.front), rear(other.rear), size(other.size)
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
			front = newSize;
			rear  = 0;
			size  = newSize;
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
			size           = other.size;
			other.capacity = other.front = other.rear = other.size = 0;
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
		if (capacity != InitialSize())
		{
			ptr      = allocator.Reallocate(ptr, InitialSize(), capacity);
			capacity = front = rear = size = 0;
		}
	}

	Iterator begin()              { return { ptr, capacity, rear  }; }
	Iterator end()                { return { ptr, capacity, front }; }
	ConstIterator cbegin() const  { return { ptr, capacity, rear  }; }
	ConstIterator cend()   const  { return { ptr, capacity, front }; }

	bool Any()   const { return Size() > 0;  }
	bool Empty() const { return Size() == 0; }

	template<typename... Args>
	ValueT& Emplace(Args&&... args)
	{
		GrowIfNecessary(1);
		new (ptr + front) ValueT(Forward<Args>(args)...);
		uint idx = front;
		front = IncrementIndex(front);
		size++;
		return ptr[idx];
	}

	void Enqueue(const ValueT& value)
	{
		GrowIfNecessary(1);
		ptr[front] = value;
		front = IncrementIndex(front);
		size++;
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
		size  += count;
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
		size -= count;
		return true;
	}

	void Dequeue(ValueT* result, uint count)
	{
		ASSERT(Size() + count <= capacity);
		uint r = rear, e = capacity-1;
		for (uint i = 0; i < count; i++)
		{
			result[i] = Forward<ValueT>(ptr[r++]);
			r &= e; // better than modulo 
		}
		rear = r;
		size -= count;
	}

	ValueT Dequeue()
	{	
		ASSERT(rear != front);
		ValueT val = Forward<ValueT>(ptr[rear]);
		rear = IncrementIndex(rear);
		size--;
		return val; 
	}

	bool TryDequeue(ValueT& out)
	{
		if (AX_UNLIKELY(rear == front)) 
			return false;
		out = Forward<ValueT>(ptr[rear]);
		rear = IncrementIndex(rear);
		size--;
		return true;
	}

	uint Size() const { return size; }
	
	// todo: shrink to fit

private:

	uint IncrementIndex(uint x) const
	{
		return (x + 1) & (capacity-1);
	}

	void GrowIfNecessary(uint _size)
	{
		ASSERT(capacity != (1 << 31)); // max size
		uint newSize = NextPowerOf2((int)(size + _size));
		
		if (AX_LIKELY(newSize <= capacity))
		{
			return; // no need to grow
		}

		const uint newCapacity = newSize <= InitialSize() ? InitialSize() : newSize;
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
			ptr[newCapacity - 1 - i] = Forward<ValueT>(ptr[capacity - 1 - i]);
		}
		capacity = newCapacity;

		rear     = newCapacity - size; 
		front    = 0;
		capacity = newCapacity;
	}
};


/////   PriorityQueue /////


enum PQCompare
{
	PQCompare_Less    = 0,
	PQCompare_Greater = 1
};

// https://github.com/lemire/FastPriorityQueue.js/blob/master/FastPriorityQueue.js

template<typename T, 
	PQCompare compare = PQCompare_Less,
	typename AllocatorT = Allocator<T>>
class PriorityQueue 
{
public:
	// we don't want to use same initial size for all data types because we want 
	// more initial size for small data types such as byte, short, int but for bigger value types we want less initial size
	static const int InitialSize = AllocatorT::InitialSize;
	
	T*         heap     = nullptr;
	int        size     = 0;
	int        capacity = 0;
	AllocatorT allocator;
	
	PriorityQueue() : heap(nullptr), size(0), capacity(0) {}
	
	~PriorityQueue()
	{
		Clear();
	}

	explicit PriorityQueue(int _size) : size(0), capacity(CalculateArrayGrowth(_size)) 
	{
		heap = allocator.AllocateUninitialized(capacity);
	}

	PriorityQueue(const T* begin, const T* end) : heap(nullptr), size(0), capacity(0) 
	{
		Push(begin, end); 
	}

	PriorityQueue& operator=(const PriorityQueue& other) 	
	{
	    if (this != &other) 
		{
		    if (other.size > capacity) 
				GrowIfNecessarry(other.size - capacity);
			Copy(heap, other.heap, other.size);
	        size = other.size;
	    }
	    return *this;
	}
private:

	static bool Compare(const T& a, const T& b)
	{
		if_constexpr (compare == PQCompare_Less)    return a < b;
		if_constexpr (compare == PQCompare_Greater) return a > b;
		else { ASSERT(0); return 0; }
	}

	void GrowIfNecessarry(int adition)
	{
		if (size + adition >= capacity)
		{
			int newCapacity = size + adition <= InitialSize ? InitialSize : CalculateArrayGrowth(size + adition);
			if (heap)
				heap = allocator.Reallocate(heap, capacity, newCapacity);
			else
				heap = allocator.Allocate(newCapacity);
			capacity = newCapacity;
		}
	}

	void HeapifyUp(int index, bool force = false) 
	{
		while (index != 0)
		{
			int parent = (index - 1) >> 1;
			
			if (Compare(heap[parent], heap[index]))
			{
				Swap(heap[parent], heap[index]);
				index = parent;
			}
			else break;
		}
	}

	void HeapifyDown(int index) 
	{
		while (true)
		{
			int leftChild  = (index << 1) + 1;
			int rightChild = leftChild + 1;
			int largest    = index;

			if (leftChild < size && !Compare(heap[leftChild], heap[largest]))
				largest = leftChild;

			if (rightChild < size && !Compare(heap[rightChild], heap[largest]))
				largest = rightChild;

			if (largest != index) 
			{
				Swap(heap[index], heap[largest]);
				index = largest;
			}
			else break;
		}
	}

public:

	bool Empty() const { return size == 0; }

	void Clear()
	{
		if (heap)
		{
			allocator.Deallocate(heap, capacity);
			heap = nullptr;
			size  = capacity = 0;
		}
	}

	void Push(const T& value) 
	{
		GrowIfNecessarry(1);
		heap[size++] = value; 
		HeapifyUp(size - 1); 
	}
	
	T RemoveAt(int index)
	{
		ASSERT(!(index > this->size - 1 || index < 0));
		T res = (T&&)heap[index];
		HeapifyUp(index, true);
		size--;
		return res;
	}

	template<typename... Args>
	void Emplace(Args&&... args) 
	{
		GrowIfNecessarry(1);
		T t(Forward<Args>(args)...);
		heap[size] = (T&&)t;
		HeapifyUp(size++);
	}

	void Pop()
	{
		ASSERT(!Empty());
		Swap(heap[0], heap[--size]);
		heap[size].~T();
		HeapifyDown(0);
	}

	const T& Top() const {
		ASSERT(!Empty());
		return heap[0];
	}

#ifdef ASTL_STL_COMPATIBLE
	bool empty() const { return size == 0; }
	void push(const T& value) { Push(value); }
	template<typename... Args>
	void emplace(Args && ... args) { Emplace(Forward<Args>(args)...); }
	void pop() { Pop(); }
	void clear() { Clear(); }
	const T& top() const { return Top(); }
#endif
};

template<typename T, int capacity>
using StackQueue = Queue<T, StackAllocator<T, capacity>>;

template<typename T, PQCompare compare, int capacity>
using StackPriorityQueue = PriorityQueue<T, compare, StackAllocator<T, capacity>>;

AX_END_NAMESPACE