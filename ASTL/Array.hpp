#pragma once

#include "Memory.hpp"
#include "Algorithms.hpp"

template<typename ValueT, 
	     typename AllocatorT = Allocator<ValueT>>
class Array
{
 public:
    using Iterator = ValueT*;
    using ConstIterator = const ValueT*;

    // we don't want to use same initial size for all data types because we want 
    // more initial size for small data types such as byte, short, int but for bigger value types we want less initial size
    // if this initial size is big for your needs, use static array please.[ byte    |  InitialSize ]
    static constexpr int InitialSize = 384 / Min((int)sizeof(ValueT), 128);//  1         384
                                                                           //  2         192
private:                                                                   //  4         96
    ValueT* arr    = nullptr;                                              //  8         48
    int m_count    = 0;
    int m_capacity = 0;
    AllocatorT allocator{};

public:
	Array() : arr(nullptr), m_count(0), m_capacity(0), allocator()
	{ }
	
	~Array()
	{
		if (arr != nullptr)
		{
			allocator.Deallocate(arr, m_capacity);
			m_count    = 0;
			m_capacity = 0;
			arr        = nullptr;
		}
	}
	
	explicit Array(const AllocatorT& _allocator) : arr(nullptr), m_count(0), m_capacity(0), allocator(_allocator)
	{ }
	
	explicit Array(Array const& other) : Array(other.Data(), other.Size()) {}
	
	explicit Array(uint _capacity) : m_capacity(_capacity)
	{
		arr = allocator.AllocateUninitialized(m_capacity);
	}
	
	Array(uint _capacity, uint _count) 
	: m_capacity(_capacity), m_count(_count)
	{
		arr = allocator.Allocate(m_capacity);
	}
	
	Array(uint _count, const ValueT& val) 
	: m_capacity(_count + (_count / 2)), m_count(_count)
	{
		arr = allocator.AllocateUninitialized(m_capacity);
		FillN(arr, val, m_count);
	}
	
	Array(ConstIterator _begin, ConstIterator _end)
	{
		m_count    = PointerDistance(_begin, _end);
		m_capacity = m_count + (m_count / 2);
		arr        = allocator.AllocateUninitialized(m_capacity);
		Copy(arr, _begin, m_count);
	}
	
	Array(ConstIterator _begin, int _size) 
	  : m_count(_size), m_capacity(_size + (_size / 2))
	{
		arr = allocator.AllocateUninitialized(m_capacity);
		Copy(arr, _begin, _size);
	}
	
	Array& operator = (Array const& other)
	{
		if (&other != this)
		{
			GrowIfNecessary(other.Size());
			Copy(arr, other.arr, other.Size());
			m_count = other.Size();
		}
		return *this;
	}
	
	Array& operator = (Array&& other) noexcept
	{
		if (&other != this)
		{
		    if (arr != nullptr)
		        allocator.Deallocate(arr, m_capacity);
			arr              = other.arr;
			m_count          = other.Size();
			m_capacity       = other.m_capacity;
			other.m_capacity = 0;
			other.m_count    = 0;
			other.arr        = nullptr;
		}
		return *this;
	}

	bool operator == (const Array& other) const
	{
		if (&other == this) return true;
		if (other.Size() != Size()) return false;
		
		for (int i = 0; i < other.Size(); ++i)
		{
			if (arr[i] != other.arr[i]) 
				return false;
		}
		return true;
	}

	bool operator != (const Array& other) const
	{
		return !(other == *this);
	}

	ValueT& operator[](int index) 
	{
		ASSERT(index >= 0 && index < m_count);
		return arr[index]; 
	}
	
	const ValueT& operator[](int index) const 
	{
		ASSERT(index >= 0 && index < m_count);
		return arr[index]; 
	}

	Array  operator + (const ValueT& other) const { Array _new = *this; _new.Add(other); return _new; }

	Array& operator += (const ValueT& type)  { Add(type);                           return *this; } 

	Array& operator += (const Array& other)  { Add(other.cbegin(), other.Size());   return *this; }

	void Add(const ValueT& value)
	{
		GrowIfNecessary(m_count + 1);
		arr[m_count++] = value;
	}

	void Add(ConstIterator _begin, int _size)
	{
		ASSERT(m_count > INT32_MAX - _size); // array max size reached
		GrowIfNecessary(m_count + _size);
		for (int i = 0; i < _size; i++)
		{
			arr[m_count++] = *_begin++;
		}
	}

	void Add(ConstIterator _begin, ConstIterator _end)
	{
		Add(_begin, PointerDistance(_begin, _end));
	}

	void Insert(int index, const ValueT& value)
	{
		OpenSpace(index, 1);
		arr[index] = value;
	}

	void Insert(ConstIterator _it, const ValueT* _begin, const ValueT* _end)
	{
		Insert(PointerDistance(begin(), _it), _begin, _end);
	}

	void Insert(int index, const ValueT* _begin, const ValueT* _end)
	{
		uint size = PointerDistance(_begin, _end);
		OpenSpace(index, size);
		Copy(arr, _begin, size);
	}

	void InsertUnordered(int index, const ValueT& value)
	{
		GrowIfNecessary(m_count + 1);
		arr[m_count++] = Forward<ValueT>(arr[index]);
		arr[index] = value;
	}

	// <summary> returns number of elements removed </summary>
	template<typename Func_t>
	int RemoveAll(const Func_t& match)
	{
		int freeIndex = 0;   // the first free slot in items array

		// Find the first item which needs to be removed.
		while (freeIndex < m_count && match(arr[freeIndex]) == false) freeIndex++;
		if (freeIndex >= m_count) return 0;

		int current = freeIndex + 1;
		while (current < m_count)
		{
			// Find the first item which needs to be kept.
			while (current < m_count && match(arr[current]) == true) current++;

			if (current < m_count)
			{
				// copy item to the free slot.
				arr[freeIndex++] = Forward<ValueT>(arr[current++]);
			}
		}

		int numCleared = m_count - freeIndex;
		ClearN(arr + freeIndex, numCleared);
		m_count = freeIndex;
		return numCleared; // removed item count
	}

	template<typename... Args>
	ValueT& EmplaceAt(int index, Args&&... args)
	{
		OpenSpace(index, 1);
		arr[index].~ValueT();
		new(&arr[index]) ValueT(Forward<Args>(args)...);
		return arr[index];
	}

	template<typename... Args>
	ValueT& EmplaceBack(Args&&... args)
	{
		GrowIfNecessary(m_count + 1);
		new(&arr[m_count]) ValueT(Forward<Args>(args)...);
		return arr[m_count++];
	}

	void Clear()
	{
		Resize(0);
	}

	void RemoveAt(Iterator _it, int _count = 1)
	{
		RemoveSpace(PointerDistance(begin(), _it), _count);
	}

	void RemoveAt(int _index, int _count = 1)
	{
		RemoveSpace(_index, _count);
	}

	void RemoveBack() { RemoveSpace(m_count - 1, 1); }

	void Remove(const ValueT& val)
	{
		int index = IndexOf(arr, val, m_count);
		if (index != -1)
		{
			RemoveSpace(index, 1);
		}
	}

	// faster than remove
	void RemoveUnordered(int index)
	{
		arr[index].~ValueT();
		Swap(arr[index], arr[--m_count]);
	}

	void Resize(int _size) 
	{
		if (_size > m_capacity)
		{
			GrowIfNecessary(_size);
		}
		else if (_size < m_capacity)
		{
			if constexpr (!AllocatorT::IsPOD)
			{
				for (int i = m_count - 1; i >= _size; --i)
				{
					arr[i].~ValueT();
				}
			}
			ReduceIfNecessary(_size);
		}
		m_count = _size;
	}

	// you can use this function after removing elements in the array, this will reduct the amount of memory used
	void ShrinkToFit() 
	{
		int newCapacity = CalculateArrayGrowth(m_count);
		if (m_capacity > newCapacity)
		{
			arr = allocator.Reallocate(arr, m_capacity, newCapacity);
		}
	}

	const ValueT& Back()     const { return arr[m_count - 1];  }
	ValueT&       Back()           { return arr[m_count - 1];  }
	void          PopBack()        { arr[--m_count].~ValueT(); }
	
	const ValueT& Front()    const { return arr[0];          }
	ValueT&       Front()          { return arr[0];          }
	void          PopFront()       { RemoveAt(0);            }
	
	int           Size()     const { return m_count;         }
	int           Capacity() const { return m_capacity;      }
	ValueT*       Data()           { return arr;             }
	const ValueT* Data()     const { return arr;             }
	bool          Empty()    const { return m_count == 0;    }
	                         
	ConstIterator cbegin()   const { return arr;            }
	ConstIterator cend()     const { return arr + m_count;  }
	ConstIterator end()      const { return arr + m_count;  }
	ConstIterator begin()    const { return arr;            }
	Iterator      begin()          { return arr;            }
	Iterator      end()            { return arr + m_count;  }
	
	#ifdef ASTL_STL_COMPATIBLE
	using iterator = ValueT*; 
	using const_iterator = const ValueT*;
	void erase(Iterator _it, int _count = 1) { Remove(_it, _count); }
	void insert(int index, const ValueT& value) { Insert(index, value); }
	void push_back(const ValueT& value) { Add(value); }
	template<typename... Args>
	ValueT& emplace_back(Args&&... args) 
	{ return EmplaceBack(Forward<Args>(args)...); }
	int size()     const { return m_count; }
	int capacity() const { return m_capacity; }
	ValueT* data() { return arr; }
	const ValueT* data() const { return arr; }
	bool empty() const { return size == 0; }
	const ValueT& back() const { return arr[m_count - 1]; }
	ValueT& back() { return arr[m_count - 1]; }
	ValueT pop_back() { return arr[--m_count]; }
	void resize(int _size) { Resize(_size); }
	void shrink_to_fit() { ShrinkToFit(); }
	#endif

private:

	void GrowIfNecessary(int _size)
	{
		ASSERT(_size <= INT32_MAX);
		if (AX_UNLIKELY(_size > m_capacity))
		{
			size_t oldCapacity = m_capacity;
			m_capacity = Max(CalculateArrayGrowth(_size), InitialSize);
			
			if (arr) // array can be nullptr (first initialization)
				arr = allocator.Reallocate(arr, oldCapacity, m_capacity);
			else
				arr = allocator.Allocate(m_capacity);
		}
	}

	void OpenSpace(int _index, int _count)
	{
		long i = m_count + _count;
		int  j = m_count;
		ASSERT(i <= INT32_MAX);

		GrowIfNecessary((int)i);
		while (j >= _index)
		{
			arr[--i] = Forward<ValueT>(arr[--j]);
		}
		m_count += _count;
	}

	bool ReduceIfNecessary(int _newSize)
	{
		constexpr uint32 minBytesToAllocate    = 32u * 1024u; // 32 kilobyte of memory quarter of your cpu L1 chace (2023)
		constexpr uint32 minElementsToReduce   = (uint32)Max(minBytesToAllocate / sizeof(ValueT), (size_t)InitialSize);
		const int halfCap                      = m_capacity / 2;
		
		// you can use _newSize as new capacity or you can use m_capacity / 2
		// but I choose half way between _newSize and halfCap because _newSize might be too low
		_newSize = Max(minElementsToReduce, uint32(CalculateArrayGrowth(_newSize) + halfCap) / 2);
		
		// if size is too small no need to reallocate
		// resize if requested size is smaller than (halfCap - CalculateArrayGrowth(_newSize)) 
		if ((_newSize > minElementsToReduce && _newSize < halfCap) ||
			// initial capacity requested and current capacity is not equal to initial size (if size is not already min)
			(m_capacity != minElementsToReduce && _newSize == minElementsToReduce)) 
		{
			arr = allocator.Reallocate(arr, _newSize, m_capacity);
			m_capacity = _newSize;
			return true;
		}
		return false;
	}

	void RemoveSpace(int _index, int _count)
	{
		ASSERT((_index + _count) <= m_count);
		int newSize = Max(m_count - _count, 0);

		int i = _index;
		int j = _index + _count;

		// *******i****j***   < opens space with incrementing both of the pointers (two pointer algorithm)
		while (j < m_count)
		{
			arr[i++] = Forward<ValueT>(arr[j++]);
		}
		
		bool reduced = ReduceIfNecessary(newSize);

		if constexpr (!AllocatorT::IsPOD)
		{
			if (!reduced)
			{
				for (int i = newSize; i < m_count; ++i)
				{
					arr[i].~ValueT();
				}
			}
		}
		m_count = newSize;
	}
};