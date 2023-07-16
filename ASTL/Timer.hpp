#pragma once
#include <stdio.h>
#include "Common.hpp"

#ifndef NDEBUG
#	define CSTIMER(message) Timer timer = Timer(message);
#else
#   define CSTIMER(message) Timer timer = Timer(message);
#endif

FINLINE uint64_t ReadTimeStampCount()
{
#ifdef _MSC_VER
    return __rdtsc();
#else
    uint64_t tsc;
    __asm__ volatile("rdtsc" : "=A" (tsc));
    return tsc;
#endif
}

struct Timer
{
	uint64_t start_point;

	const char* message;

	Timer(const char* _message) : message(_message)
	{
		start_point = ReadTimeStampCount();
	}

	~Timer()
	{
		printf("%s speed ticks: %llu\n", message, ReadTimeStampCount() - start_point);
	}
};