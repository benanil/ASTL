// Copyright (c) 2015 Jeff Preshing
// directly coppied from here: 
//https://github.com/cameron314/cpp11-on-multicore/blob/master/common/sema.h

#pragma once

#include <cassert>

#if defined(_WIN32)

extern "C" {
struct _SECURITY_ATTRIBUTES;
__declspec(dllimport) void* __stdcall CreateSemaphoreW(_SECURITY_ATTRIBUTES* lpSemaphoreAttributes, long lInitialCount, long lMAXimumCount, const wchar_t* lpName);
__declspec(dllimport) int __stdcall CloseHandle(void* hObject);
__declspec(dllimport) unsigned long __stdcall WaitForSingleObject(void* hHandle, unsigned long dwMilliseconds);
__declspec(dllimport) int __stdcall ReleaseSemaphore(void* hSemaphore, long lReleaseCount, long* lpPreviousCount);
}

class Semaphore
{
private:
	void* m_hSema;

	Semaphore(const Semaphore& other) = delete;

	Semaphore& operator=(const Semaphore& other) = delete;

public:
	Semaphore(int initialCount = 0)
	{
		assert(initialCount >= 0);
		const long maxLong = 0x7fffffff;
		m_hSema = CreateSemaphoreW(nullptr, initialCount, maxLong, nullptr);
	}

	~Semaphore()
	{
		CloseHandle(m_hSema);
	}

	void wait()
	{
		const unsigned long infinite = 0xffffffff;
		WaitForSingleObject(m_hSema, infinite);
	}

	void signal(int count = 1)
	{
		ReleaseSemaphore(m_hSema, count, nullptr);
	}
};

#elif defined(__unix__)
//---------------------------------------------------------
// Semaphore (POSIX, Linux)
//---------------------------------------------------------

#include <semaphore.h>

class Semaphore
{
private:
	sem_t m_sema;

	Semaphore(const Semaphore& other) = delete;
	Semaphore& operator=(const Semaphore& other) = delete;

public:
	Semaphore(int initialCount = 0)
	{
		assert(initialCount >= 0);
		sem_init(&m_sema, 0, initialCount);
	}

	~Semaphore()
	{
		sem_destroy(&m_sema);
	}

	void wait()
	{
		// http://stackoverflow.com/questions/2013181/gdb-causes-sem-wait-to-fail-with-eintr-error
		int rc;
		do
		{
			rc = sem_wait(&m_sema);
		} while (rc == -1 && errno == EINTR);
	}

	void signal()
	{
		sem_post(&m_sema);
	}

	void signal(int count)
	{
		while (count-- > 0)
		{
			sem_post(&m_sema);
		}
	}
};
#endif