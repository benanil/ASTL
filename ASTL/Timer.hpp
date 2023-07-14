#pragma once
#include <stdio.h>
#include "Common.hpp"

#ifndef NDEBUG
#	define CSTIMER(message) Timer timer = Timer(message);
#else
#   define CSTIMER(message) Timer timer = Timer(message);
#endif

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