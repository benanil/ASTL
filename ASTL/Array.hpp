#pragma once
#include "Memory.hpp"
#include "Algorithms.hpp"

template<typename ValueT,
         typename AllocatorT = Allocator<ValueT>>
class Array
{
	using Iterator = ValueT*;
	using ConstIterator = const ValueT*;
	using iterator = ValueT*; // stl compatible
	using const_iterator = const ValueT*;

	static constexpr int InitialSize = 8;

private:
	ValueT* arr = nullptr;
	int count = 0;
	int m_capacity = InitialSize;
	AllocatorT allocator{};
public:
	Array()
	{
	    arr = allocator.Allocate(m_capacity);
	}
	
	~Array()
	{
		if (arr != nullptr)
		{
			allocator.Deallocate(arr, m_capacity);
			count = 0;
			m_capacity = InitialSize;
			arr = nullptr;
		}
	}
	
	explicit Array(const AllocatorT& _allocator) : allocator(_allocator)
	{
		arr = allocator.Allocate(m_capacity);
	}
	
	explicit Array(Array const& other) : Array(other.Data(), other.Size()) {}
	
	explicit Array(uint _capacity) : m_capacity(_capacity)
	{
		arr = allocator.Allocate(m_capacity);
	}
	
	Array(uint _capacity, uint _count) 
	: m_capacity(_capacity), count(_count)
	{
		arr = allocator.Allocate(m_capacity);
	}
	
	Array(uint _count, const ValueT& val) 
	: m_capacity(_count + (_count / 2)), count(_count)
	{
		arr = allocator.AllocateUninitialized(m_capacity);
		FillN(arr, val, count);
	}
	
	Array(ConstIterator _begin, ConstIterator _end)
	{
		count = PointerDistance(_begin, _end);
		m_capacity  = count + (count / 2);
		arr = allocator.Allocate(m_capacity);
		Copy(arr, _begin, count);
	}
	
	Array(ConstIterator _begin, int _size) 
	  : count(_size), m_capacity(_size + (_size / 2))
	{
		arr = allocator.Allocate(m_capacity);
		Copy(arr, _begin, _size);
	}
	
	Array& operator = (Array const& other)
	{
		if (&other != this)
		{
			GrowIfNecessary(other.Size());
			Copy(arr, other.arr, other.Size());
			count = other.Size();
	  }
	  return *this;
	}
	
	Array& operator = (Array&& other) noexcept
	{
		if (&other != this)
		{
			arr = other.arr;
			count = other.Size();
			m_capacity = other.m_capacity;
			other.m_capacity = InitialSize;
			other.count = 0;
			other.arr = nullptr;
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
		ax_assert(index >= 0 && index < count);
		return arr[index]; 
	}
	
	const ValueT& operator[](int index) const 
	{
		ax_assert(index >= 0 && index < count);
		return arr[index]; 
	}

	Array operator+(const ValueT& other)
	{
		Array _new = *this;
		_new.Add(other);
		return _new;
	}

	Array& operator+=(const ValueT& type)
	{
		Add(type);
		return *this;
	}

	Array& operator+=(const Array& other)
	{
		Add(other.cbegin(), other.Size());
		return *this;
	}

	void Add(const ValueT& value)
	{
		GrowIfNecessary(count + 1);
		arr[count++] = value;
	}

	void Add(ConstIterator _begin, int _size)
	{
		GrowIfNecessary(count + _size);
		for (int i = 0; i < _size; i++)
		{
			arr[count++] = *_begin++;
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

	void InsertUninitialized(int index)
	{
		OpenSpace(index, 1);
		arr[index].~ValueT();
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
		GrowIfNecessary(count + 1);
		int end = count + 1;
		arr[end] = Forward<ValueT>(arr[index]);
		arr[index] = value;
	}

	// <summary> returns number of elements removed </summary>
	template<typename Func_t>
	int RemoveAll(const Func_t& match)
	{
		int freeIndex = 0;   // the first free slot in items array

		// Find the first item which needs to be removed.
		while (freeIndex < count && match(arr[freeIndex]) == false) freeIndex++;
		if (freeIndex >= count) return 0;

		int current = freeIndex + 1;
		while (current < count)
		{
			// Find the first item which needs to be kept.
			while (current < count && match(arr[current]) == true) current++;

			if (current < count)
			{
				// copy item to the free slot.
				arr[freeIndex++] = Forward<ValueT>(arr[current++]);
			}
		}

		int numCleared = count - freeIndex;
		ClearN(arr + freeIndex, numCleared);
		count = freeIndex;
		return numCleared; // removed item count
	}

	template<typename... Args>
	ValueT& EmplaceAt(int index, Args&&... args)
	{
		OpenSpace(index, 1);
		arr[count].~ValueT();
		new(&arr[index]) ValueT(Forward<Args>(args)...);
		return arr[index];
	}

	template<typename... Args>
	ValueT& EmplaceBack(Args&&... args)
	{
		GrowIfNecessary(count + 1);
		arr[count].~ValueT();
		new(&arr[count]) ValueT(Forward<Args>(args)...);
		return arr[count++];
	}

	template<typename... Args>
	ValueT& emplace_back(Args&&... args)
	{
		return EmplaceBack(Forward<Args>(args)...);
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

	void Remove(const ValueT& val)
	{
		int index = IndexOf(arr, val, count);
		if (index != -1)
		{
			RemoveSpace(index, 1);
		}
	}

	// faster than remove
	void RemoveUnordered(int index)
	{
		arr[index].~ValueT();
		SwapPrimitive(arr[index], arr[--count]);
	}

	void Resize(int _size) 
	{
		if (_size > m_capacity)
		{
			GrowIfNecessary(_size);
		}
		else if (_size < m_capacity)
		{
			for (int i = count - 1; i >= _size; --i)
			{
				arr[i].~ValueT();
			}
			ReduceIfNecessary(_size);
		}
		count = _size;
	}

	// you can use this function after removing elements in the array, this will reduct the amount of memory used
	void ShrinkToFit() 
	{
		ReduceIfNecessary(count);
	}

	const ValueT& Back() const { return arr[count - 1]; }
	ValueT& Back() { return arr[count - 1]; }
	void PopBack() { arr[--count].~ValueT(); }
	
	const ValueT& Front() const { return arr[0]; }
	ValueT& Fack()    { return arr[0]; }
	void PopFront() { RemoveAt(0); }
	
	int Size()     const { return count;     }
	int Capacity() const { return m_capacity; }
	ValueT* Data() { return arr; }
	const ValueT* Data() const { return arr; }
	bool Empty() const { return count == 0; }
	
	ConstIterator  cbegin() const { return arr; }
	ConstIterator  cend()   const { return arr + count; }
	ConstIterator  end()    const { return arr + count; }
	ConstIterator  begin()  const { return arr; }
	Iterator  begin() { return arr; }
	Iterator  end() { return arr + count; }
	
	#ifdef ASTL_STL_COMPATIBLE
	void erase(Iterator _it, int _count = 1) { Remove(_it, _count); }
	void insert(int index, const ValueT& value) { Insert(index, value); }
	void push_back(const ValueT& value) { Add(value); }
	int size()     const { return count; }
	int capacity() const { return m_capacity; }
	ValueT* data() { return arr; }
	const ValueT* data() const { return arr; }
	bool empty() const { return size == 0; }
	const ValueT& back() const { return arr[count - 1]; }
	ValueT& back() { return arr[count - 1]; }
	ValueT pop_back() { return arr[--count]; }
	void resize(int _size) { Resize(_size); }
	void shrink_to_fit() { ShrinkToFit(); }
	#endif

private:

	void GrowIfNecessary(int _size)
	{
		if (AX_UNLIKELY(_size > m_capacity))
		{
			size_t oldCapacity = m_capacity;
			m_capacity = CalculateGrowth(_size);
			arr = allocator.Reallocate(arr, oldCapacity, m_capacity);
		}
	}

	void OpenSpace(int _index, int _count)
	{
		int i = count + _count;
		int j = Max(count, 0);

		GrowIfNecessary(i);
		while (j >= _index)
		{
			arr[--j] = Forward<ValueT>(arr[--i]);
		}
		count += _count;
	}

	bool ReduceIfNecessary(int _newSize)
	{
		// leave space 
		_newSize = CalculateGrowth(_newSize);

		int oldCapacity = m_capacity;
		// only reallocate if new size is greater than initial size(8)
		// and requested size is smaller than (half of the capacity + growth addition) 
		if (_newSize > InitialSize && _newSize < (m_capacity / 2))
		{
			// you can use _newSize as new capacity or you can use m_capacity / 2
			// but I choose half way between _newSize and m_capacity / 2 because _newSize might be too low
			m_capacity = (_newSize + (m_capacity / 2)) / 2;  
			arr = allocator.Reallocate(arr, oldCapacity, m_capacity);
			return true;
		}
		else if (_newSize <= InitialSize && m_capacity != InitialSize) // probably 0, check if is it already initial size
		{
			m_capacity = InitialSize;  
			arr = allocator.Reallocate(arr, oldCapacity, m_capacity);
			return true;
		}
		return false;
	}

	void RemoveSpace(int _index, int _count)
	{
		int newSize = Max(count - _count, 0);

		int i = _index;
		int j = _index + _count;

		// *******i****j***   < opens space with incrementing both of the pointers (two pointer algorithm)
		while (j < count)
		{
			arr[i++] = Forward<ValueT>(arr[j++]);
		}
		
		if (!ReduceIfNecessary(newSize))
		{
			for (int i = newSize; i < count; ++i)
			{
				arr[i].~ValueT();
			}
		}
		count = newSize;
	}

	static int CalculateGrowth(int _size)
	{
		const int addition = _size / 2;
		if (AX_UNLIKELY(_size > (INT32_MAX - addition))) {
			return INT32_MAX; // growth would overflow
		}
		return _size + addition; // growth is sufficient
	}
};