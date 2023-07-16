#pragma once

inline bool IsNumber(char a) { return a <= '9' && a >= '0'; };
inline bool IsLower(char a) { return a >= 'a' && a <= 'z'; };
inline bool IsUpper(char a) { return a >= 'A' && a <= 'Z'; };
inline bool ToLower(char a) { return a < 'a' ? a + ('A' - 'a') : a; }
inline bool ToUpper(char a) { return a > 'Z' ? a - 'a' + 'A' : a; }
// is alphabetical character?
inline bool IsChar(char a) { return IsUpper(a) || IsLower(a); };
inline bool IsWhitespace(char c) { return (c == ' ' || c == '\t' || c == '\r'); }

template<typename T>
inline void Swap(T& a, T& b)
{
	T temp = (T&&)a;
	a = (T&&)b;
	b = (T&&)temp;
}

template<typename T>
inline void BubbleSort(T* arr, int len)
{
    for (int i = 0; i < len-1; ++i)
    {
        for (int j = 0; j < len-i-1; ++j)
        {
            if (arr[j] > arr[j + 1])
                Swap(arr[j], arr[j + 1]);
        }
    }
}

// smaller stack size compared to quicksort
template<typename T>
inline void ShellSort(T* arr, int n)
{
	for (int gap = n / 2; gap > 0; gap /= 2)
	{
		for (int i = gap; i < n; ++i)
		{
			T temp = arr[i];
			int j = i;
			for (; j >= gap && arr[j - gap] > temp; j -= gap)
				arr[j] = arr[j - gap];

			arr[j] = temp;
		}
	}
}


// left is beginning of the list(usually 0)
// right is end of the list usually size -1
// for simd: https://github.com/WojciechMula/simd-sort
template<typename T>
inline void QuickSort(T* arr, int left, int right)
{
	int i, j;
	while (right > left)
	{
		j = right;
		i = left - 1;
		T v = arr[right];

		while (true)
		{
			do i++; while (arr[i] < v && i < j);
			do j--; while (arr[j] > v && i < j);

			if (i >= j) break;
			Swap(arr[i], arr[j]);
		}

		Swap(arr[i], arr[right]);

		if ((i - 1 - left) <= (right - i - 1))
		{
			QuickSort(arr, left, i - 1);
			left = i + 1;
		}
		else
		{
			QuickSort(arr, i + 1, right);
			right = i - 1;
		}
	}
}


template<typename T>
inline void Reverse(T* begin, T* end)
{
	while (begin < end)
	{
		Swap(*begin, *end);
		begin++;
		end--;
	}
}

// this is classic binary search. or you can use this:
// optimized for integers. https://www.youtube.com/watch?v=1RIPMQQRBWk
// Optimizing Binary Search - Sergey Slotin - CppCon 2022
template<typename T>
inline T* BinarySearch(T* begin, int len, T value)
{
	int low = 0;
	int high = len;

	while (low < high)
	{
        T mid = (low + high) / 2;
        if (begin[mid] < value) low = mid + 1;
        else if (begin[mid] > value) high = mid - 1;
        else return begin + mid; // (begin[mid] == value)
	}
	return nullptr;
}

// String to number functions
inline int ParseNumber(const char*& curr)
{
	while (*curr && (*curr != '-' && !IsNumber(*curr))) curr++; // skip whitespace
    
    int val = 0l;
	bool negative = false;
	
    if (*curr == '-') curr++, negative = true;
	
    while (*curr > '\n' && IsNumber(*curr))
		val = val * 10 + (*curr++ - '0');
    
    return negative ? -val : val;
}

inline bool IsParsable(const char*& curr)
{   // additional checks
	if (*curr == 0 || *curr == '\n') return false;
	while (IsWhitespace(*curr)) curr++;
	if (!IsNumber(*curr) || *curr != '-') return false;
	return true;
}

inline void ParseFloat(float* f, char*& ptr)
{
    while (IsWhitespace(*ptr)) ptr++;
	
    double sign = 1.0;
    if(*ptr == '-') sign = -1.0, ptr++; 

    double num = 0.0;

    while (IsNumber(*ptr)) 
        num = 10.0 * num + (double)(*ptr++ - '0');

    if (*ptr == '.') ptr++;

    double fra = 0.0, div = 1.0;

    while (IsNumber(*ptr))
        fra = 10.0f * fra + (double)(*ptr++ - '0'), div *= 10.0f;
	
    num += fra / div;
    *f = (float)(sign * num);
}

inline int Pow10(int x) 
{
  int res = x != 0;
  while (x-- > 0) res *= 10;
  return res;
}

inline int Log10(int n) // float version of this is declared in math.h
{
  return n < 10 ? 0 : 1 + Log10(n / 10);
}

// time complexity O(numDigits(x)), space complexity O(1)
// @returns number of characters added
inline int IntToString(char* ptr, int x, int afterPoint = 0)
{
    int len = Log10(x);
    int blen = len;
    int size = 0;

    while (++blen < afterPoint)
    {
        ptr[size++] = '0';
        ptr[size] = '\0';
    }
    
    len = Pow10(len) + 1;
   
    while (len)
    {
        int digit = x / len;
        ptr[size++] = char(digit + '0');
        ptr[size] = '\0';
        x -= len * digit;
        len /= 10;
    }
    return size;
}

// converts floating point to string
// @param afterpoint num digits after friction 
// @returns number of characters added
inline int FloatToString(char* ptr, float f, int afterpoint = 4)
{
    int iPart = (int)f;
    int numChars = IntToString(ptr, iPart);
    float fPart = f - iPart;
    ptr[numChars++] = '.';
    int power = Pow10(afterpoint); 
    return numChars + IntToString(ptr + numChars, int(fPart * power), afterpoint);
}

inline bool StartsWith(const char*& curr, const char* str)
{
    const char* currStart = curr;
    while (IsWhitespace(*curr)) curr++;
    if (*curr != *str) return false;
    while (*str && *curr++ == *str++);
    bool isEqual = *str == 0;
    if (!isEqual) curr = currStart;
    return isEqual;
}

template<typename T> inline void Fill(T* begin, T* end, const T& val)
{
	while (begin < end) *begin++ = val;
}

template<typename T> inline void FillN(T* arr, T val, int n)
{
	for (int i = 0; i < n; ++i) arr[i] = val;
}

template<typename T> inline bool Contains(const T* arr, const T& val, int n)
{
    for (int i = 0; i < n; ++i) 
        if (arr[i] == val) return true;
    return false;
}

// returns the index if found, otherwise returns -1
template<typename T> inline int IndexOf(const T* arr, const T& val, int n)
{
    for (int i = 0; i < n; ++i) 
        if (arr[i] == val) return i;
    return -1;
}

// returns number of val contained in arr
template<typename T> inline int CountIf(const T* arr, const T& val, int n)
{
    int count = 0;
    for (int i = 0; i < n; ++i) 
        count += arr[i] == val;
    return count;
}

template<typename T> inline void Copy(T* dst, const T* src, int n)
{
	for (int i = 0; i < n; ++i)
		dst[i] = src[i];
}

template<typename T> inline void MoveArray(T* dst, T* src, int n)
{
  for (int i = 0; i < n; ++i)
	dst[i] = (T&&)src[i];
}

template<typename T> inline void ClearN(T* src, int n)
{
    for (int i = 0; i < n; ++i)
        src[i].~T();
}
