#pragma once
#ifdef ASTL_MULTI_THREADING

#include <atomic>
#include <functional>
#include "Common.hpp"
#include "../Array.hpp"

template<typename T>
struct ParallelExecutable : public IExecutable
{
	void Execute() override // overrided function
	{
		for (int i = 0; i < m_Len; ++i)
			m_Func(m_Data[i]);
	}

	void Set(T* data, int len, void(* func)(T&))
	{
		m_Data = data;
		m_Len = len;
		m_Func = func;
	}

	T* m_Data = nullptr;
	int m_Len = 0;

	void (* m_Func)(T&) = nullptr;
};

enum eParallelState
{
	eParallelState_None, eParallelState_Done = 1u, eParallelState_Requested = 2u
};

struct ParallelThread
{
	IExecutable* fun = nullptr;
	std::atomic<unsigned> state;
};

// this class spawns n threads and you can use this on whatever data you want
// lets say we are using a parallel for each frame we can use this class for that
// this will wait until we need this other frame that way we will be avoiding thread generation each frame
class ParallelFor
{
public:
	static constexpr unsigned MaxThreads = 16u;
	static constexpr unsigned MaxJobsPerThread = 16u;

	explicit ParallelFor(int numThreads);

	~ParallelFor();

	// spawns n threads and executes job parallely on data
	template<typename T>
	static void ExecuteOnce(const int numThreads, T* data, const int size, std::function<void(T&)> f, bool waitUntilFinish = true)
	{
		const int itemPerThread = size / (numThreads + waitUntilFinish);
		const T* end = data + size;

		auto threadJob = [f, itemPerThread](T* data) -> void
		{
			for (int i = 0; i < itemPerThread; ++i)
				f(data[i]);
		};

		std::thread threads[MaxThreads];

		for (int i = 0; i < numThreads; ++i)
		{
			new(threads + i)std::thread(threadJob, data);
			data += itemPerThread;
		}
		if (waitUntilFinish)
		{
			while (data < end)
			{ // work in this thread too
				f(*data++);
			}
			// wait until all threads finished
			for (int i = 0; i < numThreads; ++i)
			{
				threads[i].join();
			}
		}
	}

	template<typename T>
	void Execute(T* data, int size, std::function<void(T&)> f, bool waitUntilFinish = true)
	{
		const int itemPerThread = size / (m_NumThreads + waitUntilFinish);
		const T* end = data + size;

		ParallelExecutable<T> tasks[MaxThreads];

		for (int i = 0; i < m_NumThreads; ++i)
		{
			// spin lock if thread is bussy
			// todo: look at other threads instead of spinning here
			while (m_Threads[i].state.load() & eParallelState_Requested)
				std::this_thread::yield();

			tasks[i].Set(data, itemPerThread, f);

			m_Threads[i].fun = dynamic_cast<IExecutable*>(&tasks[i]);
			m_Threads[i].state |= eParallelState_Requested;
			data += itemPerThread;
		}

		if (waitUntilFinish)
		{
			// work in this thread too
			while (data < end)
			{
				f(*data++);
			}

			// wait until all threads finished
			int numFinishedTasks = 0;
			while (numFinishedTasks < m_NumThreads)
			{
				numFinishedTasks = 0;
				for (int i = 0; i < m_NumThreads; ++i)
				{
					numFinishedTasks += m_Threads[i].state.load() != eParallelState_Requested;
				}
				std::this_thread::yield();
			}
		}
	}

	const uint32 GetNumThreads() const
	{ return m_NumThreads; }

private:
	ParallelThread m_Threads[MaxThreads] = {};
	unsigned m_NumThreads;
};

#endif