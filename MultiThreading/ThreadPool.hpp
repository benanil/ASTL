#pragma once
#ifdef ASTL_MULTI_THREADING

#include "Common.hpp"
#include <functional>
#include <atomic>

template<int NumJobs = 16>
class LockFreeRingBuffer
{
public:
	LockFreeRingBuffer() : m_NumAdded(0), m_NumProcessed(0)  { }

	bool Push(const std::function<void()>& value)
	{
		// maximum jobs reached
		if (m_NumAdded.load() - m_NumProcessed.load() == NumJobs)
			return false;

		m_Data[m_NumAdded++ % NumJobs] = value;
		return true;
	}

	bool Pop(std::function<void()>& value)
	{
		// if we don't have any job return false
		if (!(m_NumAdded.load() > m_NumProcessed.load())) return false;
		value = m_Data[m_NumProcessed++ % NumJobs];
		return true;
	}

private:
	std::function<void()> m_Data[NumJobs] = {};
	std::atomic<int> m_NumProcessed;
	std::atomic<int> m_NumAdded;
public:
	// we don't need to make this atomic because we will write once and
	// most of the time this will be false, if we miss this data in producer thread, this will be not a problem because it will yield and read this value properly
	bool m_Done = false;
};

class ThreadPool
{
public:
	static const uint32 MAXThreads = 16u;
	static const uint32 MAXJobsPerThread = 16u;

	ThreadPool();
	explicit ThreadPool(uint32 threadCount = 4u);
	~ThreadPool();

	void Destroy();

	// use std::bind for arguments
	void PushJob(std::function<void()> job);

	int GetNumThreads() const { return m_NumThreads; }

private:
	void Initialize();
private:
	int m_NumThreads = 4;
	int m_CurretWorkerIndex = 0;
	LockFreeRingBuffer<MAXJobsPerThread> m_Jobs[MAXThreads];
};
#endif