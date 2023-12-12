
/****************************************************************
*    Purpose: Usefull General algorithms, parsing, sorting etc. *
*    Author : Anilcan Gulkaya 2023 anilcangulkaya7@gmail.com    *
****************************************************************/

#pragma once

#ifdef AX_USE_NAMESPACE
namespace ax {
#endif

inline constexpr bool IsNumber(char a) { return a <= '9' && a >= '0'; };
inline constexpr bool IsLower(char a) { return a >= 'a' && a <= 'z'; };
inline constexpr bool IsUpper(char a) { return a >= 'A' && a <= 'Z'; };
inline constexpr bool ToLower(char a) { return a < 'a' ? a + ('A' - 'a') : a; }
inline constexpr bool ToUpper(char a) { return a > 'Z' ? a - 'a' + 'A' : a; }
// is alphabetical character?
inline constexpr bool IsChar(char a) { return IsUpper(a) || IsLower(a); };
inline constexpr bool IsWhitespace(char c) { return c <= ' '; }

template<typename T>
inline constexpr void Swap(T& a, T& b)
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
        T mid = (low + high) >> 1;
        if (begin[mid] < value) low = mid + 1;
        else if (begin[mid] > value) high = mid - 1;
        else return begin + mid; // (begin[mid] == value)
	}
	return nullptr;
}

// String to number functions
inline int ParseNumber(const char*& ptr)
{
  const char* curr = ptr;
  while (*curr && (*curr != '-' && !IsNumber(*curr))) curr++; // skip whitespace
    
  int val = 0l;
  bool negative = false;
  
  if (*curr == '-') curr++, negative = true;
  
  while (*curr > '\n' && IsNumber(*curr))
    val = val * 10 + (*curr++ - '0');
  
  ptr = curr;
  return negative ? -val : val;
}

inline int ParsePositiveNumber(const char*& ptr)
{
    const char* curr = ptr;
    while (*curr && !IsNumber(*curr)) curr++; // skip whitespace

    int val = 0;
    while (*curr > '\n' && IsNumber(*curr))
        val = val * 10 + (*curr++ - '0');
    ptr = curr;
    return val;
}

inline bool IsParsable(const char* curr)
{   // additional checks
	if (*curr == 0 || *curr == '\n') return false;
	if (!IsNumber(*curr) || *curr != '-') return false;
	return true;
}

inline float ParseFloat(const char*& text)
{
    const char* ptr = text;
    while (!IsNumber(*ptr) && *ptr != '-') ptr++;
	
    double sign = 1.0;
    if(*ptr == '-') sign = -1.0, ptr++; 

    double num = 0.0;

    while (IsNumber(*ptr)) 
        num = 10.0 * num + (double)(*ptr++ - '0');

    if (*ptr == '.') ptr++;

    double fra = 0.0, div = 1.0;

    while (IsNumber(*ptr) && div < 1e8) // 1e8 is 1 and 8 zero 100000000
        fra = 10.0f * fra + (double)(*ptr++ - '0'), div *= 10.0f;
	
    while (IsNumber(*ptr)) ptr++;

    num += fra / div;
    text = ptr;
    return (float)(sign * num);
}

#ifndef AX_NO_UNROLL
    #if defined(__clang__)
    #   define AX_NO_UNROLL _Pragma("clang loop unroll(disable)") _Pragma("clang loop vectorize(disable)")
    #elif defined(__GNUC__) >= 8
    #   define AX_NO_UNROLL _Pragma("GCC unroll 0")
    #elif defined(_MSC_VER)
    #   define AX_NO_UNROLL __pragma(loop(no_vector))
    #else
    #   define AX_NO_UNROLL
    #endif  
#endif

inline unsigned int Log10Algo(unsigned int n)
{                                                                        
	if (n <= 9) return 0;
	if (n <= 99) return 1;
	if (n <= 999) return 2;
	if (n <= 9999) return 3;
	if (n <= 99999) return 4;
	if (n <= 999999) return 5;
	if (n <= 9999999) return 6;
	if (n <= 99999999) return 7;
	if (n <= 999999999) return 8;
	if (n <= 2147483647) return 9; // 2147483647 = int max
	return 0;
}                                                                                        
                                                                     
// time complexity O(numDigits(x)), space complexity O(1)                                
// @returns number of characters added                                                   
inline int IntToString(char* ptr, int x, int afterPoint=0)                                      
{                                                                                        
    int size = 0;
    if (x < 0) ptr[size++] = '-', x = 0-x;
    
	int numDigits = Log10Algo(x);
    int blen = numDigits;
    
    AX_NO_UNROLL while (++blen <= afterPoint)
    {
        ptr[size++] = '0';
    }

	unsigned int const PowersOf10[10] = { 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000 };
	numDigits = PowersOf10[numDigits];

    while (numDigits)
    {
        int digit = x / numDigits;
        ptr[size++] = (char)(digit + '0');
        x -= numDigits * digit;
        numDigits /= 10;
    }

    ptr[size] = '\0';
    return size;
}

inline int Pow10(int x) 
{
	int res = x != 0;
	while (x-- > 0) res *= 10;
	return res;
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
    return numChars + IntToString(ptr + numChars, int(fPart * power), afterpoint-1);
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

template<typename T> inline void ConstructN(T* src, int n)
{
    T def{};
    for (int i = 0; i < n; ++i) 
        src[i] = def;
}

#ifdef AX_USE_NAMESPACE
}
#endif