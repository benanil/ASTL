#ifdef ASTL_MULTI_THREADING

#include "ParallelFor.hpp"
#include <iostream>

static void AXParallelWorkerProc(ParallelThread* parallelThread)
{
	while ((parallelThread->state & eParallelState_Done) == 0)
	{
		if (parallelThread->state & eParallelState_Requested)
		{
			parallelThread->fun->Execute();
			parallelThread->state &= ~eParallelState_Requested;
		}
		std::this_thread::yield();
	}
}

ParallelFor::ParallelFor(int numThreads) : m_NumThreads(numThreads)
{
	for (int i = 0; i < numThreads; ++i)
	{
		std::thread worker(AXParallelWorkerProc, &m_Threads[i]);
		SetToUniqueCore(worker, i);
		worker.detach();
	}
}

ParallelFor::~ParallelFor() {
	for (unsigned i = 0; i < m_NumThreads; ++i)
		m_Threads[i].state = eParallelState_Done;
}
#endif