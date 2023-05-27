///////////////////////// ankerl::unordered_dense::{map, set} /////////////////////////

// A fast & densely stored hashmap and hashset based on robin-hood backward shift deletion.
// Version 4.0.0
// https://github.com/martinus/unordered_dense
//
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.
// SPDX-License-Identifier: MIT
// Copyright (c) 2022-2023 Martin Leitner-Ankerl <martin.ankerl@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include "Random.hpp"
#include "Array.hpp"

/*
* example custom hasher
template<> struct Hasher<int> 
{
    static FINLINE uint64 Hash(int obj) {
        return uint64(Random::WangHash(x)) * 0x9ddfea08eb382d69ull;
    }
};
*/
template<typename KeyT,
         typename HasherT = Hasher<KeyT>,
         typename AllocatorT = Allocator<KeyT>>
class HashSet
{
    using Iterator = KeyT*;
    using ConstIterator = const KeyT*;
    
    struct Bucket
    {
        static constexpr uint32_t DistInc = 1U << 8U;             // skip 1 byte fingerprint
        static constexpr uint32_t FingerprintMask = DistInc - 1u; // mask for 1 byte of fingerprint
        uint32_t distAndFingerprint; // upper 3 byte: distance to original bucket. lower byte: fingerprint from hash
        uint32_t valueIdx;            // index into the m_values vector.
    };

    static constexpr uint8_t initial_shifts = 64 - 3; // 2^(64-m_shift) number of buckets
    static constexpr float default_max_load_factor = 0.8F;

    // todo use custom allocator and our vector
    Array<KeyT, AllocatorT> m_keys{};
    Array<Bucket, MallocAllocator<Bucket>> m_buckets(16);
    uint32 m_num_buckets = 0;
    uint32 m_max_bucket_capacity = 0;
    float m_max_load_factor = default_max_load_factor;
    uint8_t m_shifts = initial_shifts;

private:
    uint32 Next(uint bucketIdx) const {
        return AX_UNLIKELY(bucketIdx + 1u == m_num_buckets) ? 0u : bucketIdx + 1u;
    }

    uint32 DistInc(uint32 x) const { return x + Bucket::DistInc; }

    uint32 DistDec(uint32 x) const { return x - Bucket::DistInc; }

    uint32 DistAndFingerprintFromHash(uint64 hash) const {
        return Bucket::DistInc | (hash & Bucket::FingerprintMask);
    }

    Bucket& BucketAt(uint32 index) { return m_buckets[index]; }
    const Bucket& BucketAt(uint32 index) const { return m_buckets[index]; }

    uint32 BucketIdxFromHash(uint64 hash) const { return uint32(hash >> m_shifts); }
    
    Bucket NextWhileLess(const KeyT& key) const
    {
        uint64 hash = HasherT::Hash(key);
        uint32 distAndFingerprint = DistAndFingerprintFromHash(hash);
        uint32 bucketIdx= BucketIdxFromHash(hash);

        while (distAndFingerprint < BucketAt(bucketIdx).distAndFingerprint)
        {
            distAndFingerprint = DistInc(distAndFingerprint);
            bucketIdx = Next(bucketIdx);
        }
        return { distAndFingerprint, bucketIdx };
    }

    void PlaceAndShiftUp(Bucket bucket, uint32 place)
    {
        while (0 != BucketAt(place).distAndFingerprint)
        {
            bucket = Exchange(BucketAt(place), bucket);
            bucket.distAndFingerprint = DistInc(bucket.distAndFingerprint);
            place = Next(place);
        }
        m_buckets[place] = bucket;
    }

    float LoadFactor() const {
        return m_num_buckets ? float(m_keys.Size()) / float(m_num_buckets) : 0.0F;
    }

    constexpr uint32 MaxSize() const { return 1ull << (sizeof(uint) * 8u - 1u); }

    constexpr uint32 CalcNumBuckets(uint8 shifts)
    {
        return Min(MaxSize(), 1u << (64u - shifts));
    }

    constexpr uint8 CalcShiftsForSize(uint32 s)
    {
        uint8 shifts = initial_shifts;

        while (shifts > 0 && uint32(float(CalcNumBuckets(shifts)) * m_max_load_factor) < s)
            --shifts;
        
        return shifts;
    }

    void MaxLoadFactor(float ml) {
        m_max_load_factor = ml;
        if (m_num_buckets != MaxSize()) {
            m_max_bucket_capacity = uint32(float(m_num_buckets) * m_max_load_factor);
        }
    }

    void CopyBuckets(const HashSet& other)
    {
        if (!IsEmpty()) return;

        m_shifts = other.m_shifts;
        ReallocateBuckets(CalcNumBuckets(m_shifts));
        MemCpy(&m_buckets[0], &other.m_buckets[0], m_num_buckets * sizeof(Bucket));
    }

    void ReallocateBuckets(uint32 numBuckets)
    {
        m_num_buckets = numBuckets;
        m_buckets.Resize(m_num_buckets);
        
        if (m_num_buckets == MaxSize()) {
            m_max_bucket_capacity = MaxSize();
        }
        else {
            m_max_bucket_capacity = uint32(float(m_num_buckets) * m_max_load_factor);
        }
    }

    void ClearAndFillBucketsFromValues() 
    {
        MemSet(m_buckets.Data(), 0, m_num_buckets * sizeof(Bucket));

        for (uint32 value_idx = 0u,
            end_idx = uint32(m_keys.Size());
            value_idx < end_idx;
          ++value_idx) 
        {
            const KeyT&  key = m_keys[value_idx];
            Bucket bucket = NextWhileLess(key);

            // we know for certain that key has not yet been inserted, so no need to check it.
            PlaceAndShiftUp({bucket.distAndFingerprint, value_idx}, bucket.valueIdx);
        }
    }

    void IncreaseSize()
    {
        ax_assert(AX_UNLIKELY(m_max_bucket_capacity) == MaxSize() && "bucket overflow");
        --m_shifts;
        ReallocateBuckets(CalcNumBuckets(m_shifts)); // DeallocateBuckets(); AllocateBuffersFromShift();
        ClearAndFillBucketsFromValues();
    }

    void DoErase(uint32 bucketIdx)
    {
        uint32 valueIdxToRemove = m_buckets[bucketIdx].valueIdx;
        uint32 nextBucketIdx = Next(bucketIdx);

        while (BucketAt(nextBucketIdx).distAndFingerprint >= Bucket::DistInc * 2u)
        {
            const Bucket& nextBucket = BucketAt(nextBucketIdx);
            BucketAt(bucketIdx) = { DistDec(nextBucket.distAndFingerprint), nextBucket.valueIdx };
            bucketIdx = Exchange(nextBucketIdx, Next(nextBucketIdx));
        }
        BucketAt(bucketIdx) = {};
        
        if (valueIdxToRemove != m_keys.Size()-1)
        {
            KeyT& key = m_keys[valueIdxToRemove];
            key = Move(m_keys.Back());

            uint64 mh = HasherT::Hash(key);
            bucketIdx = BucketIdxFromHash(mh);

            const uint32 valuesIdxBack = uint32(m_keys.Size()) - 1;
            while (valuesIdxBack != BucketAt(bucketIdx).valueIdx)
            {
                bucketIdx = Next(bucketIdx);
            }
            BucketAt(bucketIdx).valueIdx = valueIdxToRemove;
        }
        m_keys.PopBack();
    }

    template<typename... Args>
    Pair<Iterator, bool> DoTryEmplace(Args&&... args)
    {
        if (AX_UNLIKELY(IsFull())) 
            IncreaseSize();

        KeyT key(Forward<Args>(args)...);
        uint64 hash = HasherT::Hash(key);
        uint32 distAndFootprint = DistAndFingerprintFromHash(hash);
        uint32 bucketIdx = BucketIdxFromHash(hash);
        Bucket bucket = BucketAt(bucketIdx);

        while (distAndFootprint <= bucket.distAndFootprint)
        {
            if (distAndFootprint == bucket.distAndFootprint && 
                key == m_keys[bucket.valueIdx])
            {
                return { begin() + bucket.valueIdx, false };
            }
            distAndFootprint = DistInc(distAndFootprint);
            bucketIdx = Next(bucketIdx);
        }
        m_keys.EmplaceBack(Forward<KeyT>(key));
        uint32 valueIdx = uint32(m_keys.Size()) - 1;
        PlaceAndShiftUp({distAndFootprint, valueIdx}, bucketIdx);
        return { begin() + valueIdx, true };
    }
    
    template<typename K>
    Pair<Iterator, bool> DoTryInsert(K&& key)
    {
        if (AX_UNLIKELY(IsFull())) 
            IncreaseSize();
    
        uint64 hash = HasherT::Hash(key);
        uint32 distAndFootprint = DistAndFingerprintFromHash(hash);
        uint32 bucketIdx = BucketIdxFromHash(hash);
    
        while (true)
        {
            Bucket bucket = BucketAt(bucketIdx);
            if (distAndFootprint == bucket.distAndFingerprint && 
                             key == m_keys[bucket.valueIdx])
            {
                return { begin() + bucket.valueIdx, false };
            }
            else if (distAndFootprint > bucket.distAndFingerprint)
            {
                m_keys.EmplaceBack(Forward<K>(key));
                uint32 valueIdx = uint32(m_keys.Size()) - 1;
                PlaceAndShiftUp({distAndFootprint, valueIdx}, bucketIdx);
                return { begin() + valueIdx, true };
            }
            distAndFootprint = DistInc(distAndFootprint);
            bucketIdx = Next(bucketIdx);
        }
    }

    ConstIterator DoFind(KeyT key) const
    {
        if (AX_UNLIKELY(IsEmpty()))
            return cend();

        uint64 mh = HasherT::Hash(key);
        uint32 distAndFingerPrint = DistAndFingerprintFromHash(mh);
        uint32 bucketIdx = BucketIdxFromHash(mh);
        
        const Bucket* bucket = &BucketAt(bucketIdx);
        if (distAndFingerPrint == bucket->distAndFingerprint &&
                                  key == m_keys[bucket->valueIdx]) {
          return cbegin() + bucket->valueIdx;
        }  

        distAndFingerPrint = DistInc(distAndFingerPrint);
        bucketIdx = Next(bucketIdx);
        bucket = &BucketAt(bucketIdx);

        if (distAndFingerPrint == bucket->distAndFingerprint &&
                                  key == m_keys[bucket->valueIdx]) {
          return cbegin() + bucket->valueIdx;
        }

        distAndFingerPrint = DistInc(distAndFingerPrint);
        bucketIdx = Next(bucketIdx);
        bucket = &BucketAt(bucketIdx);

        while (true) {
            if (distAndFingerPrint == bucket->distAndFingerprint &&
                               key == m_keys[bucket->valueIdx]) {
                return cbegin() + bucket->valueIdx;
            }
            else if (distAndFingerPrint > bucket->distAndFingerprint) {
                return cend();
            }
            distAndFingerPrint = DistInc(distAndFingerPrint);
            bucketIdx = Next(bucketIdx);
            bucket = &BucketAt(bucketIdx);
        }
    }

public:

    ConstIterator Find(KeyT const& key) const {
        return DoFind(key);
    }

    Iterator Insert(const KeyT& key) {
      return DoTryInsert(key).first;
    }
    
    bool Contains(KeyT const& key) const 
    {
      return DoFind(key) != cend();
    }

    HashSet() : HashSet(0ull) {}
 
    HashSet(uint32 bucketCount)
    {
        if (bucketCount != 0) 
            Reserve(bucketCount);
    }

    HashSet(const KeyT* begin, const KeyT* end) 
    : HashSet(0) {
        Insert(begin, end);
    }

    HashSet(const KeyT *begin, uint32 bucketCount)
    : HashSet(0) {
        Insert(begin, begin + bucketCount);
    }

    HashSet& operator = (HashSet const& other) {
        if (&other != this) {
            ReallocateBuckets(other.m_num_buckets); // deallocate before m_values is set (might have another allocator)
            m_keys = other.m_keys;
            m_max_load_factor = other.m_max_load_factor;
            m_shifts = initial_shifts;
            CopyBuckets(other);
        }
        return *this;
    }

    HashSet& operator = (HashSet&& other) noexcept 
    {
        if (&other != this) {
            ReallocateBuckets(other.m_num_buckets); // deallocate before m_values is set (might have another allocator)
            m_keys = Move(other.m_keys);
            m_buckets = Move(other.m_buckets);
            m_num_buckets = other.m_num_buckets;
            m_max_bucket_capacity = other.m_max_bucket_capacity;
            m_max_load_factor = other.m_max_load_factor;
            m_shifts = other.m_shifts;
            other.Clear();
        }
        return *this;
    }

    ConstIterator  cbegin() const { return m_keys.cbegin(); }
    ConstIterator  cend()   const { return m_keys.cend(); }
    ConstIterator  end()    const { return m_keys.cend(); }
    ConstIterator  begin()  const { return m_keys.cbegin(); }
    Iterator  begin()             { return m_keys.begin(); }
    Iterator  end()               { return m_keys.end(); }

    bool IsFull()  const { return Size() >= m_max_bucket_capacity; }
    bool IsEmpty() const { return m_keys.Size() == 0; }
    uint32 Size()  const { return m_keys.Size(); }
    uint32 size()  const { return m_keys.Size(); }

    void Clear() { m_keys.Clear(); m_buckets.Resize(16); }

    void Insert(ConstIterator first, ConstIterator last) 
    {
        while (first != last) {
            DoTryInsert(*first);
            ++first;
        }
    }
    
    template<class M>
    Iterator DoInsertOrAssign(KeyT&& key, M&& mapped)
    {
        Pair<Iterator, bool> inserted = TryEmplace(Forward<KeyT>(key), Forward<M>(mapped));
        if (inserted.second) {
            inserted.first->value = Forward<M>(mapped);
        }
        return inserted.first;
    }

    ConstIterator Erase(ConstIterator it)
    {
        uint64 hash = HasherT::Hash(it->key);
        uint32 bucketIdx = BucketIdxFromHash(hash);
        uint32 valueIdxToRemove = uint32(PointerDistance(cbegin(), it));

        while (BucketAt(bucketIdx).valueIdx != valueIdxToRemove) {
            bucketIdx = Next(bucketIdx);
        }
        DoErase(bucketIdx);
        return cbegin() + valueIdxToRemove;
    }

    uint32 Erase(const KeyT& key)
    {
        if (IsEmpty()) return 0u;
    
        const Bucket bucket = NextWhileLess(key);
    
        uint32 distAndFingerprint = bucket.distAndFingerprint;
        uint32 bucketIdx = bucket.valueIdx;
    
        while (distAndFingerprint == BucketAt(bucketIdx).distAndFingerprint &&
                              key != m_keys[BucketAt(bucketIdx).valueIdx])
        {
            distAndFingerprint = DistInc(distAndFingerprint);
            bucketIdx = Next(bucketIdx);
        }
    
        if (distAndFingerprint != BucketAt(bucketIdx).distAndFingerprint)
        {
            return 0u;
        }
        DoErase(bucketIdx);
        return 1u;
    }

    template<typename Pred>
    uint32 EraseIf(Pred pred)
    {
        uint32 oldSize = Size();
        uint32 idx = oldSize;
        while (idx--)
        {
            Iterator it = begin() + idx;
            if (pred(*it)) Erase(it);
        }
        
        return Size() - oldSize;
    }

    void ReHash(uint32 count)
    {
        count = Min(count, MaxSize());
        uint8 shifts = CalcShiftsForSize(Size());
        if (shifts != m_shifts)
        {
            m_shifts = shifts;
            ReallocateBuckets();
            ClearAndFillBucketsFromValues();
        }
    }
    
    void Reserve(uint32 capacity) 
    {
        capacity = Min((uint32)capacity, MaxSize());
        m_keys.Resize(capacity);
        uint8 shifts = CalcShiftsForSize(Max(capacity, Size()));
        
        if (0 == m_num_buckets || shifts < m_shifts) {
            m_shifts = shifts;
            ReallocateBuckets();
            ClearAndFillBucketsFromValues();
        }
    }

    template<typename... Args>
    Iterator TryEmplace(Args&&... args)
    {
        return DoTryEmplace(Move<Args>(args)...).first;
    }

    friend bool operator == (const HashSet& a, const HashSet& b)
    {
        if (&a == &b) return true; 
        if (a.Size() != b.Size())  return false; 
        
        for (const KeyT& b_entry : b) 
        {
            ConstIterator it = a.Find(b_entry.key);
            // map: check that key is here.
            if (a.cend() == it) {
                return false;
            }
        }
        return true;
    }

    friend bool operator != (const HashSet& a, const HashSet& b) {
        return !(a == b);
    }
#ifdef ASTL_STL_COMPATIBLE
    // stl compatibility
    bool contains(KeyT const& key) const
    {
      return DoFind(key) != cend();
    }

    Iterator insert(const KeyT& key) {
      return DoTryInsert(key).first;
    }
#endif
};