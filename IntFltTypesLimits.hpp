#pragma once

#if !defined(_STDINT) || !defined(_INTTYPES)
typedef unsigned char      uint8   ;
typedef unsigned char      uint8_t ;
typedef unsigned short     uint16  ;
typedef unsigned short     uint16_t;
typedef unsigned int       uint32  ;
typedef unsigned int       uint32_t;
typedef unsigned long long uint64;
typedef unsigned long long uint64_t;
typedef char      int8   ;
typedef signed char int8_t;
typedef short     int16  ;
typedef short     int16_t;
typedef int       int32  ;
typedef int       int32_t;
typedef long long int64  ;
typedef long long int64_t;

typedef signed long long int intptr_t;
typedef unsigned long long int uintptr_t;

#define INT8_MIN   (-127 - 1)
#define INT16_MIN  (-32767 - 1)
#define INT32_MIN  (-2147483647 - 1)
#define INT64_MIN  (-9223372036854775807 - 1)
#define INT8_MAX   127
#define INT16_MAX  32767
#define INT32_MAX  2147483647
#define INT64_MAX  9223372036854775807
                   
#define UINT8_MAX  0xffui8
#define UINT16_MAX 0xffffui16
#define UINT32_MAX 0xffffffffui32
#define UINT64_MAX 0xffffffffffffffffui64

#define INT8_C(x)    (x)
#define INT16_C(x)   (x)
#define INT32_C(x)   (x)
#define INT64_C(x)   (x ## LL)

#define UINT8_C(x)   (x)
#define UINT16_C(x)  (x)
#define UINT32_C(x)  (x ## U)
#define UINT64_C(x)  (x ## ULL)

#define INTMAX_C(x)  INT64_C(x)
#define UINTMAX_C(x) UINT64_C(x)

#endif // defined stdint.h

typedef uint8_t  uchar;
typedef uint16_t ushort;
typedef uint32_t uint;
typedef uint64_t ulong;


#ifndef _INC_FLOAT // include guard for 3rd party interop

#define DBL_MAX          1.7976931348623158e+308 // max value
#define DBL_MIN          2.2250738585072014e-308 // min positive value

#define FLT_MAX          3.402823466e+38F        // max value
#define FLT_MIN          1.175494351e-38F        // min normalized positive value

#endif