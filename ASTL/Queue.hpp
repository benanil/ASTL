#pragma once
#include "Memory.hpp"

// todo: fix this
template<typename T>
class Queue
{
public:
	~Queue()
	{
		if (ptr)
		{
			Clear();
			free(ptr);
			ptr = nullptr;
		}
	}

	Queue() : front(0), rear(0), capacity(32)
	{
		ptr = (T*)calloc(capacity, sizeof(T));
	}

	Queue(int _size) : front(0), rear(0), capacity(_size)
	{
		ptr = (T*)calloc(capacity, sizeof(T));
	}

	// copy constructor
	Queue(const Queue& other)
		: capacity(other.capacity),
		ptr((T*)calloc(capacity, sizeof(T))),
		front(other.front),
		rear(other.rear)
	{
		memcpy(ptr, other.ptr, capacity);
	}
	// move constructor
	Queue(Queue&& other)
		: capacity(other.capacity),
		ptr(other.ptr),
		front(other.front),
		rear(other.rear)
	{
		other.ptr = nullptr;
	}

	void Clear()
	{
		memset(ptr, 0, sizeof(T) * capacity);
		rear = front = 0;
	}

	T* begin() { return ptr + rear; }
	T* end() { return ptr + front; }
	const T* cbegin() const { return ptr + rear; }
	const T* cend()   const { return ptr + front; }

	inline bool Any() const { return GetSize() > 0; }

	void Enqueue(T value)
	{
		if (front + 1 >= capacity)
		{
			capacity += 32;
			const int size = GetSize();
			ptr = (T*)realloc(ptr, capacity * sizeof(T));
			memmove(ptr, ptr + rear, size * sizeof(T));
			rear = 0;
			front = size;
		}
		ptr[front++] = value;
	}

	void Enqueue(T* begin, T* end)
	{
		const int count = PointerDistance(begin, end) * sizeof(T);

		if (front + count >= capacity)
		{
			capacity += (capacity / 2) + count;
			const int size = GetSize();
			ptr = (T*)realloc(ptr, capacity * sizeof(T));
			memmove(ptr, ptr + rear, size * sizeof(T));
			rear = 0;
			front = size;
		}
		memcpy(ptr + rear + front, begin, count * sizeof(T));
	}

	// returns true if size is enough
	bool Dequeue(T** result, int count)
	{
		if (GetSize() < count) return false;
		*result = (T*)malloc(sizeof(T) * count);
		memcpy(result, ptr + rear, sizeof(count));
		rear += count;
		return true;
	}

	T Dequeue() {
		return ptr[rear++];
	}

	bool TryDequeue(T& out)
	{
		if (GetSize() > 0) out = Dequeue();
		return GetSize();
	}

public:
	T* ptr;
	inline int GetSize() const { return front - rear; }
private:
	int capacity;
	int front = 0;
	int rear = 0;
};