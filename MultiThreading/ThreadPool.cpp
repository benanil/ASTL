#ifdef ASTL_MULTI_THREADING

#include "ThreadPool.hpp"
#include "../Math/Math.hpp"

static void AXThreadProc(LockFreeRingBuffer<16>* jobs)
{
	while (!jobs->m_Done)
	{
		std::function<void()> job; 
		while (jobs->Pop(job)) {
			job();
		}
		std::this_thread::yield();
	}
}

void ThreadPool::Destroy()
{
	for (int i = 0; i < m_NumThreads; ++i) {
		m_Jobs[i].m_Done = true;
	}
}

void ThreadPool::PushJob(std::function<void()> job)
{
	// if thread is full try to find available thread and push to it
	while (!m_Jobs[m_CurretWorkerIndex++ % m_NumThreads].Push(job))
		std::this_thread::yield();
}

ThreadPool::~ThreadPool() {

}

ThreadPool::ThreadPool()
: m_NumThreads(Clamp(std::thread::hardware_concurrency(), 1u, MaxThreads))
{
	Initialize();
}

ThreadPool::ThreadPool(uint32 threadCount)
: m_NumThreads(Clamp(threadCount, 1u, MaxThreads))
{
	Initialize();
}

void ThreadPool::Initialize()
{
	for (int i = 0; i < m_NumThreads; ++i)
	{
		std::thread worker(AXThreadProc, &m_Jobs[i]);
		SetToUniqueCore(worker, i);
		worker.detach();// forget about this thread, let it do it's job in the infinite loop that we created above
	}
}
#endif