
#pragma once

#include <Array.hpp>
#include <iostream>


typedef int ValueT;

class RedBlackTree 
{
private:
    struct Node 
    {
        int left;
        int right;
        ValueT data;

        bool IsRed()  const  { return !!(parent & LastBitMask); }
        bool IsBlack() const { return   (parent & LastBitMask) == 0; }

        void SetParent(uint x)
        {
            parent = x & (parent & LastBitMask);
        }

        uint GetParent() const
        {
            return parent ~ LastBitMask;
        }
    private:
        uint parent;
        constexpr uint LastBitMask = ~(0xFFFFFFFFu >> 1u);
    };

    Array<Node> nodes;
public:
  
};