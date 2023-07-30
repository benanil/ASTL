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

AX_NAMESPACE 

/* example custom hasher
template<> struct Hasher<int> 
{
    static __forceinline uint64_t Hash(int obj)  {
        return uint64(WangHash(x)) * 0x9ddfea08eb382d69ull;
    }
};
*/
template<typename KeyT,
    typename ValueT,
    typename HasherT = Hasher<KeyT>,
    typename AllocatorT = Allocator<KeyValuePair<KeyT, ValueT>>>
class HashMap
{
    using Iterator = KeyValuePair<KeyT, ValueT>*;
    using ConstIterator = const KeyValuePair<KeyT, ValueT>*;
    
    struct Bucket
    {
        static const uint32_t DistInc = 1u << 8u;             // skip 1 byte fingerprint
        static const uint32_t FingerprintMask = DistInc - 1u; // mask for 1 byte of fingerprint
        uint32_t distAndFingerprint; // upper 3 byte: distance to original bucket. lower byte: fingerprint from hash
        uint32_t valueIdx;           // index into the m_values vector.
    };

    static const uint8 initial_shifts = 64u - 3u; // 2^(64-m_shift) number of buckets

    Array<KeyValuePair<KeyT, ValueT>, AllocatorT> m_values{};
    Array<Bucket, MallocAllocator<Bucket>> m_buckets{};

    uint32 m_num_buckets         = 0u;
    uint32 m_max_bucket_capacity = 0u;
    float m_max_load_factor      = 0.8f; // max load factor
    uint8 m_shifts               = initial_shifts;

private:
    uint32 Next(uint bucketIdx) const {
        return AX_UNLIKELY(bucketIdx + 1u == m_num_buckets) ? 0u : bucketIdx + 1u;
    }

    uint32 DistInc(uint32 x) const { return x + Bucket::DistInc; }

    uint32 DistDec(uint32 x) const { return x - Bucket::DistInc; }

    uint32 DistAndFingerprintFromHash(uint64_t hash) const {
        return Bucket::DistInc | (hash & Bucket::FingerprintMask);
    }

    Bucket& BucketAt(uint32 index) { return m_buckets[index]; }
    const Bucket& BucketAt(uint32 index) const { return m_buckets[index]; }

    uint32 BucketIdxFromHash(uint64_t hash) const { return uint32(hash >> m_shifts); }
    
    Bucket NextWhileLess(const KeyT& key) const
    {
        uint64_t hash               = HasherT::Hash(key);
        uint32 distAndFingerprint = DistAndFingerprintFromHash(hash);
        uint32 bucketIdx          = BucketIdxFromHash(hash);

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
        return m_num_buckets ? float(m_values.Size()) / float(m_num_buckets) : 0.0f;
    }

    __constexpr uint32 MAXSize() const { return 1u << (sizeof(uint) * 8u - 1u); }

    uint32 CalcNumBuckets(uint8 shifts) const
    {
        return MIN(MAXSize(), 1u << (64u - shifts));
    }

    __constexpr uint8 CalcShiftsForSize(uint32 s)
    {
        uint8 shifts = initial_shifts;

        while (shifts > 0 && uint32(float(CalcNumBuckets(shifts)) * m_max_load_factor) < s)
            --shifts;
        
        return shifts;
    }

    void MaxLoadFactor(float ml) {
        m_max_load_factor = ml;
        if (m_num_buckets != MAXSize()) {
            m_max_bucket_capacity = uint32(float(m_num_buckets) * m_max_load_factor);
        }
    }

    void CopyBuckets(const HashMap& other)
    {
        if (!Empty()) return;
        m_shifts = other.m_shifts;
        ReallocateBuckets(CalcNumBuckets(m_shifts)); // AllocateBuffersFromShift();
        MemCpy<alignof(Bucket)>(&m_buckets[0], &other.m_buckets[0], m_num_buckets * sizeof(Bucket));
    }

    void ReallocateBuckets(uint32 numBuckets)
    {
        m_num_buckets = numBuckets;
        m_buckets.Resize(m_num_buckets);

        if (AX_UNLIKELY(m_num_buckets == MAXSize())) {
            m_max_bucket_capacity = MAXSize();
        }
        else {
            m_max_bucket_capacity = uint32(float(m_num_buckets) * m_max_load_factor);
        }
    }

    void ClearAndFillBucketsFromValues() 
    {
        MemSet<alignof(Bucket)>(m_buckets.Data(), 0, m_num_buckets * sizeof(Bucket));

        for (uint32 value_idx = 0u,
            end_idx = uint32(m_values.Size());
            value_idx < end_idx;
          ++value_idx) 
        {
            const KeyT&  key = m_values[value_idx].key;
            Bucket bucket = NextWhileLess(key);
            // we know for certain that key has not yet been inserted, so no need to check it.
            PlaceAndShiftUp({bucket.distAndFingerprint, value_idx}, bucket.valueIdx);
        }
    }

    void IncreaseSize()
    {
        ASSERT(m_max_bucket_capacity != MAXSize()); 
        --m_shifts;
        ReallocateBuckets(CalcNumBuckets(m_shifts));
        ClearAndFillBucketsFromValues();
    }

    void DoErase(uint32 bucketIdx)
    {
        uint32 valueIdxToRemove = m_buckets[bucketIdx].valueIdx;
        uint32 nextBucketIdx    = Next(bucketIdx);

        while (BucketAt(nextBucketIdx).distAndFingerprint >= Bucket::DistInc * 2u)
        {
            const Bucket& nextBucket = BucketAt(nextBucketIdx);
            BucketAt(bucketIdx) = { DistDec(nextBucket.distAndFingerprint), nextBucket.valueIdx };
            bucketIdx = Exchange(nextBucketIdx, Next(nextBucketIdx));
        }
        BucketAt(bucketIdx) = {};
        
        if (valueIdxToRemove != m_values.Size()-1)
        {
            KeyValuePair<KeyT, ValueT>& val = m_values[valueIdxToRemove];
            val = Move(m_values.Back());

            uint64_t mh = HasherT::Hash(val.key);
            bucketIdx = BucketIdxFromHash(mh);

            const uint32 valuesIdxBack = uint32(m_values.Size()) - 1;
            while (valuesIdxBack != BucketAt(bucketIdx).valueIdx)
            {
                bucketIdx = Next(bucketIdx);
            }
            BucketAt(bucketIdx).valueIdx = valueIdxToRemove;
        }
        m_values.PopBack();
    }

    template<typename K, typename ...Args>
    Pair<Iterator, bool> DoTryEmplace(K&& key, Args&&... args)
    {
        if (AX_UNLIKELY(IsFull())) 
            IncreaseSize();

        uint64_t hash             = HasherT::Hash(key);
        uint32 distAndFootprint = DistAndFingerprintFromHash(hash);
        uint32 bucketIdx        = BucketIdxFromHash(hash);

        while (true)
        {
            Bucket bucket = BucketAt(bucketIdx);
            if (distAndFootprint == bucket.distAndFingerprint)
            {
                if (key == m_values[bucket.valueIdx].key)
                {
                    return { begin() + bucket.valueIdx, false };
                }
            }
            else if (distAndFootprint > bucket.distAndFingerprint)
            {
                ValueT val(Forward<Args>(args)...);
                KeyValuePair<KeyT, ValueT> p((KeyT&&)key, (ValueT&&)val);
                m_values.EmplaceBack(Forward<KeyValuePair<KeyT, ValueT>>(p));
                uint32 valueIdx = uint32(m_values.Size()) - 1;
                PlaceAndShiftUp({distAndFootprint, valueIdx}, bucketIdx);
                return { begin() + valueIdx, true };
            }
            distAndFootprint = DistInc(distAndFootprint);
            bucketIdx = Next(bucketIdx);
        }
    }
    
    Pair<Iterator, bool> DoTryInsert(const KeyT& key, const ValueT& value)
    {
        if (AX_UNLIKELY(IsFull())) 
            IncreaseSize();
    
        uint64_t hash             = HasherT::Hash(key);
        uint32 distAndFootprint = DistAndFingerprintFromHash(hash);
        uint32 bucketIdx        = BucketIdxFromHash(hash);
    
        while (true)
        {
            Bucket bucket = BucketAt(bucketIdx);
            if (distAndFootprint == bucket.distAndFingerprint)
            {
                if (key == m_values[bucket.valueIdx].key)
                {
                    return { begin() + bucket.valueIdx, false };
                }
            }
            else if (distAndFootprint > bucket.distAndFingerprint)
            {
                KeyValuePair<KeyT, ValueT> p(key, value);
                m_values.EmplaceBack((KeyValuePair<KeyT, ValueT>&&)p);
                uint32 valueIdx = uint32(m_values.Size()) - 1;
                PlaceAndShiftUp({distAndFootprint, valueIdx}, bucketIdx);
                return { begin() + valueIdx, true };
            }
            distAndFootprint = DistInc(distAndFootprint);
            bucketIdx = Next(bucketIdx);
        }
    }

    ConstIterator DoFind(const KeyT& key) const
    {
        if (AX_UNLIKELY(Empty()))
            return cend();

        uint64_t mh                 = HasherT::Hash(key);
        uint32 distAndFingerPrint = DistAndFingerprintFromHash(mh);
        uint32 bucketIdx          = BucketIdxFromHash(mh);
        // first check two times without while loop. (unrolling for performance)
        const Bucket* bucket      = &BucketAt(bucketIdx);

        if (distAndFingerPrint == bucket->distAndFingerprint &&
                                  key == m_values[bucket->valueIdx].key) {
          return cbegin() + bucket->valueIdx;
        }
        
        distAndFingerPrint = DistInc(distAndFingerPrint);
        bucketIdx = Next(bucketIdx);
        bucket = &BucketAt(bucketIdx);

        if (distAndFingerPrint == bucket->distAndFingerprint &&
                                  key == m_values[bucket->valueIdx].key) {
          return cbegin() + bucket->valueIdx;
        }
        
        distAndFingerPrint = DistInc(distAndFingerPrint);
        bucketIdx = Next(bucketIdx);
        bucket = &BucketAt(bucketIdx);

        while (true) {
            if (distAndFingerPrint == bucket->distAndFingerprint &&
                               key == m_values[bucket->valueIdx].key) {
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
        return (ConstIterator)DoFind(key);
    }

    bool Contains(KeyT const& key) const 
    {
        return DoFind(key) != cend();
    }

    ValueT& At(KeyT const& key)
    {
        Iterator it = (Iterator)Find(key);
        if (AX_LIKELY(end() != it))
        {
            return it->value;
        }
        ASSERT(true); // key is not exist in array
        return m_values[0].value;
    }

    const ValueT& At(const KeyT& key) const
    {
        return const_cast<HashMap*>(this)->At(key);
    }

    HashMap() : HashMap(0ull) {}
 
    HashMap(uint32 bucketCount)
    {
        if (bucketCount != 0) 
            Reserve(bucketCount);
    }

    HashMap(HashMap const& other)
      : HashMap(other.m_values.Data(), other.m_values.Size()) {}

    HashMap(const KeyValuePair<KeyT, ValueT>* begin, 
            const KeyValuePair<KeyT, ValueT>* end)
    : HashMap(0) {
        Insert(begin, end);
    }

    HashMap(const KeyValuePair<KeyT, ValueT>* begin, size_t bucketCount)
    : HashMap(0) {
        Insert(begin, begin + bucketCount);
    }

    HashMap& operator = (HashMap const& other) {
        if (&other != this) {
            ReallocateBuckets(other.m_num_buckets); 
            m_values          = other.m_values;
            m_max_load_factor = other.m_max_load_factor;
            m_shifts          = initial_shifts;
            CopyBuckets(other);
        }
        return *this;
    }

    HashMap& operator = (HashMap&& other) noexcept 
    {
        if (&other != this) {
            ReallocateBuckets(other.m_num_buckets); 
            m_values              = Move(other.m_values);
            m_buckets             = Move(other.m_buckets);
            m_num_buckets         = other.m_num_buckets;
            m_max_bucket_capacity = other.m_max_bucket_capacity;
            m_max_load_factor     = other.m_max_load_factor;
            m_shifts              = other.m_shifts;
            other.Clear();
        }
        return *this;
    }

    ConstIterator  cbegin() const { return m_values.cbegin(); }
    ConstIterator  cend()   const { return m_values.cend(); }
    ConstIterator  end()    const { return m_values.cend(); }
    ConstIterator  begin()  const { return m_values.cbegin(); }
    Iterator  begin()             { return m_values.begin(); }
    Iterator  end()               { return m_values.end(); }

    bool IsFull()  const { return Size() >= m_max_bucket_capacity; }
    bool Empty() const { return m_values.Size() == 0; }
    uint32 Size()  const { return m_values.Size(); }

    void Clear() { m_values.Clear(); m_buckets.Resize(16); }

    Iterator Insert(const KeyT& key, const ValueT& val) {
        return DoTryInsert(key, val).first;
    }

    void Insert(const KeyValuePair<KeyT, ValueT>* first, 
                const KeyValuePair<KeyT, ValueT>* last) 
    {
        while (first != last) {
            DoTryInsert(first->key, first->value);
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
        uint64_t hash             = HasherT::Hash(it->key);
        uint32 bucketIdx        = BucketIdxFromHash(hash);
        uint32 valueIdxToRemove = uint32(PointerDistance(cbegin(), it));

        while (BucketAt(bucketIdx).valueIdx != valueIdxToRemove) {
            bucketIdx = Next(bucketIdx);
        }
        DoErase(bucketIdx);
        return cbegin() + valueIdxToRemove;
    }

    uint32 Erase(const KeyT& key)
    {
        if (Empty()) 
            return 0u;
    
        const Bucket bucket       = NextWhileLess(key);
        uint32 distAndFingerprint = bucket.distAndFingerprint;
        uint32 bucketIdx          = bucket.valueIdx;
    
        while (distAndFingerprint == BucketAt(bucketIdx).distAndFingerprint &&
                              key != m_values[BucketAt(bucketIdx).valueIdx].key)
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
        while (idx-- > 0)
        {
            Iterator it = begin() + idx;
            if (pred(*it)) Erase(it);
        }
        
        return Size() - oldSize;
    }

    void ReHash(uint32 count)
    {
        count = MIN(count, MAXSize());
        uint8 shifts = CalcShiftsForSize(Size());
        if (shifts != m_shifts)
        {
            m_shifts = shifts;
            ReallocateBuckets(CalcNumBuckets(m_shifts));
            ClearAndFillBucketsFromValues();
        }
    }
    
    void Reserve(uint32 capacity) 
    {
        capacity = MIN(capacity, MAXSize());
        m_values.Resize(capacity);
        uint8 shifts = CalcShiftsForSize(MAX(capacity, Size()));
        
        if (0 == m_num_buckets || shifts < m_shifts) {
            m_shifts = shifts;
            ReallocateBuckets(CalcNumBuckets(m_shifts));
            ClearAndFillBucketsFromValues();
        }
    }

    template<typename... Args>
    Pair<Iterator, bool> TryEmplace(KeyT const& key, Args&& ... args)
    {
        return DoTryEmplace(key, Forward<Args>(args)...);
    }

    template<typename... Args>
    Pair<Iterator, bool> TryEmplace(KeyT&& key, Args&&... args)
    {
        return DoTryEmplace(Move(key), Forward<Args>(args)...);
    }

    template<typename... Args>
    Iterator TryEmplace(ConstIterator /*hint*/, KeyT const& key, Args&&... args)
    {
        return DoTryEmplace(key, Forward<Args>(args)...).first;
    }

    template<typename... Args>
    Iterator TryEmplace(ConstIterator /*hint*/, KeyT&& key, Args&&... args)
    {
        return DoTryEmplace(Move(key), Forward<Args>(args)...).first;
    }

    ValueT& operator[](const KeyT& key) { 
        return TryEmplace(key).first->value;
    }

    ValueT& operator[](KeyT&& key) {
        return TryEmplace(Move(key)).first->value;
    }

    const ValueT& operator[](const KeyT& key) const {
        return Find(key)->value;
    }

    friend bool operator == (const HashMap& a, const HashMap& b)
    {
        if (&a == &b) return true; 
        if (a.Size() != b.Size())  return false; 
        
        for (const KeyValuePair<KeyT, ValueT>& b_entry : b) 
        {
            ConstIterator it = a.Find(b_entry.key);
            // map: check that key is here, then also check that value is the same
            if (a.cend() == it || !(b_entry.value == it->value)) {
                return false;
            }
        }
        return true;
    }

    friend bool operator != (const HashMap& a, const HashMap& b) {
        return !(a == b);
    }

#ifdef ASTL_STL_COMPATIBLE
    Iterator insert(const KeyT& key, const ValueT& val) {
        return DoTryInsert(key, val).first;
    }

    Iterator insert(const KeyValuePair<KeyT, ValueT>& keyVal) {
        return DoTryInsert(keyVal.key, keyVal.value).first;
    }

    Iterator erase(Iterator it) { return Erase(it); }
    
    bool contains(KeyT const& key) const { return DoFind(key) != cend(); }
    
    ConstIterator find(KeyT const& key) const { return DoFind(key); }

    uint32 erase(const KeyT& key) { return Erase(key); }

    ValueT& at(KeyT const& key) {
        Iterator it = (Iterator)Find(key);
        if (AX_LIKELY(end() != it)) 
            return it->value;
        ASSERT(true); // key is not exist in array
        return m_values[0].value;
    }
#endif
};

AX_END_NAMESPACE 