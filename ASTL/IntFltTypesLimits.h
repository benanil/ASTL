#pragma once

#ifndef _STDINT

#if defined(_MSC_VER)
#if _MSC_VER < 1300
typedef unsigned char      uint8   ;
typedef unsigned char      uint8_t ;
typedef unsigned short     uint16  ;
typedef unsigned short     uint16_t;
typedef unsigned int       uint32  ;
typedef unsigned int       uint32_t;
typedef unsigned long long uint64  ;
typedef unsigned long long uint64_t;
typedef char      int8   ;
typedef char      int8_t ;
typedef short     int16  ;
typedef short     int16_t;
typedef int       int32  ;
typedef int       int32_t;
typedef long long int64  ;
typedef long long int64_t;
#else
typedef unsigned __int8  uint8   ;
typedef unsigned __int8  uint8_t ;
typedef unsigned __int16 uint16  ;
typedef unsigned __int16 uint16_t;
typedef unsigned __int32 uint32  ;
typedef unsigned __int32 uint32_t;
typedef unsigned __int64 uint64  ;
typedef unsigned __int64 uint64_t;
typedef signed __int8  int8   ;
typedef signed __int8  int8_t ;
typedef signed __int16 int16  ;
typedef signed __int16 int16_t;
typedef signed __int32 int32  ;
typedef signed __int32 int32_t;
typedef signed __int64 int64  ;
typedef signed __int64 int64_t;
#endif
#else
#include <stdint.h>
typedef uint8_t  uint8 ;
typedef int8_t   int8  ;
typedef uint16_t uint16;
typedef int16_t  int16 ;
typedef uint32_t uint32;
typedef int32_t  int32 ;
typedef uint64_t unt64 ;
typedef int64_t  int64 ;
#endif

#define INT8_MIN   (-127i8 - 1)
#define INT16_MIN  (-32767i16 - 1)
#define INT32_MIN  (-2147483647i32 - 1)
#define INT64_MIN  (-9223372036854775807i64 - 1)
#define INT8_MAX   127i8
#define INT16_MAX  32767i16
#define INT32_MAX  2147483647i32
#define INT64_MAX  9223372036854775807i64
                   
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

typedef uint16_t ushort;
typedef uint32_t uint;
typedef uint64_t ulong;

#ifndef _INC_FLOAT // include guard for 3rd party interop

#define DBL_MAX          1.7976931348623158e+308 // max value
#define DBL_MIN          2.2250738585072014e-308 // min positive value

#define FLT_MAX          3.402823466e+38F        // max value
#define FLT_MIN          1.175494351e-38F        // min normalized positive value

#endif