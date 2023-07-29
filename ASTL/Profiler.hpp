#pragma once

// from Casey Muratori's Computer Enhance course

// #define AX_PROFILER_DISABLE

#include "Common.hpp"

#ifdef AX_PROFILER_DISABLE
#    define TimeBlock(Name) 
#    define TimeFunction 
#else
#    define NameConcat2(A, B) A##B
#    define NameConcat(A, B) NameConcat2(A, B)
#    define TimeBlock(Name) profile_block NameConcat(Block, __LINE__)(Name, __COUNTER__ + 1);
#    define TimeFunction TimeBlock(__func__)
#endif

struct profile_anchor
{
    uint64_t  TSCElapsedExclusive; // NOTE(casey): Does NOT include children
    uint64_t  TSCElapsedInclusive; // NOTE(casey): DOES include children
    uint64_t  HitCount;
    char const* Label;
};

struct profiler
{
    profile_anchor Anchors[2048];
    uint64_t  StartTSC;
    uint64_t  EndTSC;
};

extern profiler GlobalProfiler;
extern uint32_t GlobalProfilerParent;

struct profile_block
{
    profile_block(char const* Label_, uint32_t AnchorIndex_);
    ~profile_block(void);

    char const* Label;
    uint64_t  OldTSCElapsedInclusive;
    uint64_t  StartTSC;
    uint32_t ParentIndex;
    uint32_t AnchorIndex;
};


void PrintTimeElapsed(uint64_t  TotalTSCElapsed, profile_anchor* Anchor);

void BeginProfile(void);

void EndAndPrintProfile();
