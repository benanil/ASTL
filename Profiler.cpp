
#include "Profiler.hpp"
#include <stdio.h>

#if _WIN32


#define WIN32_LEAN_AND_MEAN 
#define NOMINMAX
#define VC_EXTRALEAN
#include <intrin.h>
#include <windows.h>

static uint64_t  GetOSTimerFreq(void)
{
    LARGE_INTEGER Freq;
    QueryPerformanceFrequency(&Freq);
    return Freq.QuadPart;
}

static uint64_t  ReadOSTimer(void)
{
    LARGE_INTEGER Value;
    QueryPerformanceCounter(&Value);
    return Value.QuadPart;
}
#else
#include <x86intrin.h>
#include <sys/time.h>

static uint64_t  GetOSTimerFreq(void)
{
    return 1000000;
}

static uint64_t  ReadOSTimer(void)
{
    // NOTE(casey): The "struct" keyword is not necessary here when compiling in C++,
    // but just in case anyone is using this file from C, I include it.
    struct timeval Value;
    gettimeofday(&Value, 0);
    uint64_t  Result = GetOSTimerFreq() * (uint64_t )Value.tv_sec + (uint64_t )Value.tv_usec;
    return Result;
}
#endif

/* NOTE(casey): This does not need to be "inline", it could just be "static"
   because compilers will inline it anyway. But compilers will warn about
   static functions that aren't used. So "inline" is just the simplest way
   to tell them to stop complaining about that. */
inline uint64_t  ReadCPUTimer(void)
{
    // NOTE(casey): If you were on ARM, you would need to replace __rdtsc
    // with one of their performance counter read instructions, depending
    // on which ones are available on your platform.
    return __rdtsc();
}

static uint64_t  EstimateCPUTimerFreq(void)
{
    uint64_t  MillisecondsToWait = 100;
    uint64_t  OSFreq = GetOSTimerFreq();

    uint64_t  CPUStart = ReadCPUTimer();
    uint64_t  OSStart = ReadOSTimer();
    uint64_t  OSEnd = 0;
    uint64_t  OSElapsed = 0;
    uint64_t  OSWaitTime = OSFreq * MillisecondsToWait / 1000;
    while (OSElapsed < OSWaitTime)
    {
        OSEnd = ReadOSTimer();
        OSElapsed = OSEnd - OSStart;
    }

    uint64_t  CPUEnd = ReadCPUTimer();
    uint64_t  CPUElapsed = CPUEnd - CPUStart;

    uint64_t  CPUFreq = 0;
    if (OSElapsed)
    {
        CPUFreq = OSFreq * CPUElapsed / OSElapsed;
    }

    return CPUFreq;
}

// extern variables
profiler GlobalProfiler{};
uint32_t GlobalProfilerParent{};

profile_block::profile_block(char const* Label_, uint32_t AnchorIndex_)
{
    ParentIndex = GlobalProfilerParent;

    AnchorIndex = AnchorIndex_;
    Label = Label_;

    profile_anchor* Anchor = GlobalProfiler.Anchors + AnchorIndex;
    OldTSCElapsedInclusive = Anchor->TSCElapsedInclusive;

    GlobalProfilerParent = AnchorIndex;
    StartTSC = ReadCPUTimer();
}

profile_block::~profile_block(void)
{
    uint64_t  Elapsed = ReadCPUTimer() - StartTSC;
    GlobalProfilerParent = ParentIndex;

    profile_anchor* Parent = GlobalProfiler.Anchors + ParentIndex;
    profile_anchor* Anchor = GlobalProfiler.Anchors + AnchorIndex;

    Parent->TSCElapsedExclusive -= Elapsed;
    Anchor->TSCElapsedExclusive += Elapsed;
    Anchor->TSCElapsedInclusive = OldTSCElapsedInclusive + Elapsed;
    ++Anchor->HitCount;

    /* NOTE(casey): This write happens every time solely because there is no
       straightforward way in C++ to have the same ease-of-use. In a better programming
       language, it would be simple to have the anchor points gathered and labeled at compile
       time, and this repetative write would be eliminated. */
    Anchor->Label = Label;
}

#ifndef AX_PROFILER_DISABLE

void PrintTimeElapsed(uint64_t  TotalTSCElapsed, profile_anchor* Anchor)
{
    double Percent = 100.0 * ((double)Anchor->TSCElapsedExclusive / (double)TotalTSCElapsed);
    printf("  %s[%llu]: %llu (%.2f%%", Anchor->Label, Anchor->HitCount, Anchor->TSCElapsedExclusive, Percent);
    if (Anchor->TSCElapsedInclusive != Anchor->TSCElapsedExclusive)
    {
        double PercentWithChildren = 100.0 * ((double)Anchor->TSCElapsedInclusive / (double)TotalTSCElapsed);
        printf(", %.2f%% w/children", PercentWithChildren);
    }
    printf(")\n");
}

void BeginProfile()
{
    GlobalProfiler.StartTSC = ReadCPUTimer();
}

void EndAndPrintProfile()
{
    GlobalProfiler.EndTSC = ReadCPUTimer();
    uint64_t  CPUFreq = EstimateCPUTimerFreq();

    uint64_t  TotalCPUElapsed = GlobalProfiler.EndTSC - GlobalProfiler.StartTSC;

    if (CPUFreq)
    {
        printf("\nTotal time: %0.4fms (CPU freq %llu)\n", 1000.0 * (double)TotalCPUElapsed / (double)CPUFreq, CPUFreq);
    }

    const uint32_t numAnchors = sizeof(GlobalProfiler.Anchors) / sizeof(profile_anchor);

    for (uint32_t AnchorIndex = 0; AnchorIndex < numAnchors; ++AnchorIndex)
    {
        profile_anchor* Anchor = GlobalProfiler.Anchors + AnchorIndex;
        if (Anchor->TSCElapsedInclusive)
        {
            PrintTimeElapsed(TotalCPUElapsed, Anchor);
        }
    }
    memset(&GlobalProfiler, 0, sizeof(profiler));
    GlobalProfilerParent = 0;
}
#else


void PrintTimeElapsed(uint64_t  TotalTSCElapsed, profile_anchor* Anchor) {}


void BeginProfile() {}

void EndAndPrintProfile() {}

#endif

