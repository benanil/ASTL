#pragma once
#include "Common.hpp"

template<typename T>
class Stack
{
public:
	~Stack() { if (ptr) { Clear();  free(ptr);  ptr = nullptr; } }

	Stack()
		: size(0), capacity(32) {
		ptr = (T*)calloc(capacity, sizeof(T));
	}
	Stack(int _size)
		: size(0), capacity(_size) {
		ptr = (T*)calloc(capacity, sizeof(T));
	}

	// copy constructor
	Stack(const Stack& other) : size(other.size), capacity(other.capacity), ptr((char*)calloc(capacity, sizeof(T)))
	{
		memcpy(ptr, other.ptr, size);
	}
	// move constructor 
	Stack(Stack&& other) : size(other.size), capacity(other.capacity), ptr(other.ptr)
	{
		other.ptr = nullptr;
	}

	T& operator[](int index) { return ptr[index]; }                       
	const T& operator[](int index) const { return ptr[index]; }           
	T* begin() { return ptr; }								              
	T* end() { return ptr + size; }							              
	const T* cbegin() const { return ptr; }					              
	const T* cend() const { return ptr + size; }				          
	T GetLast() { return ptr[size - 1]; }                                 
	T GetFirst(){ return ptr[0]; }                                        
	inline bool Any() { return size > 0; }								  
	void Clear()  { memset(ptr, 0, capacity * sizeof(T)); size = 0; }

	void Push(T value) {
		if (size + 1 == capacity) {
			capacity += 32;
			ptr = (T*)realloc(ptr, capacity * sizeof(T));
		}
		ptr[size++] = value;
	}

	T Pop() { return ptr[--size]; }
	bool TryPop(T& out) {
		if (size > 0) out = ptr[--size];
		return size > 0;
	}

public:
	int size;
	int capacity;
	T* ptr;
};

template<typename T, int MaxSize>
class FixedStack
{
private:
	int size = 0;
public:
	T arr[MaxSize];
	constexpr FixedStack() { }

	constexpr void Push(const T& value) { arr[size++] = value; }
	constexpr T& Pop() { return arr[--size]; }

	constexpr bool Any() const { return size > 0; }
	
	constexpr T& operator[](int index) { return arr[index]; }
	constexpr T GetLast() { return arr[size - 1]; }
	constexpr T GetFirst() { return arr[0]; }
	constexpr void Clear() { memset(arr, 0, MaxSize * sizeof(T)); size = 0; }

	constexpr T* begin() { return arr; }
	constexpr T* end() { return arr + size; }
	constexpr const T* cbegin() const { return arr; }
	constexpr const T* cend()   const { return arr + size; }
};

