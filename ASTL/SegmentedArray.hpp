#pragma once
#include "Memory.hpp"

// Very much like std::deque, but faster for indexing (in most cases). As of now this doesn't implement the full std::vector
// API, but merely what's necessary to work as an underlying container for ankerl::unordered_dense::{map, set}.
// It allocates blocks of equal size and puts them into the m_blocks vector. That means it can grow simply by adding a new
// block to the back of m_blocks, and doesn't double its size like an std::vector. The disadvantage is that memory is not
// linear and thus there is one more indirection necessary for indexing.
/*
template<typename T, typename AllocatorT = Allocator<T>, size_t MaxSegmentSizeBytes = 4096>
class SegmentedVector
{
private:
    // todo use our array
    // todo custom allocator
    std::vector<T*> m_blocks{};
    size_t m_size{};

    static constexpr size_t NumBitsClosest(size_t maxVal, size_t s)
    {
        size_t f = 0ull;
        while (s << (f + 1) <= maxVal) {
            f++;
        }
        return f;
    }

    [[nodiscard]] static constexpr size_t CalcNumBlocksForCapacity(size_t capacity) {
        return (capacity + num_elements_in_block - 1U) / num_elements_in_block;
    }

public:
    // using vec_alloc = typename std::allocator_traits<Allocator>::template rebind_alloc<pointer>;
    using SelfT = SegmentedVector<T, AllocatorT, MaxSegmentSizeBytes>;
    static constexpr size_t num_bits = NumBitsClosest(MaxSegmentSizeBytes, sizeof(T));
    static constexpr size_t num_elements_in_block = 1ull << num_bits;
    static constexpr size_t mask = num_elements_in_block - 1ull;

    class iter_t
    {
        T* m_data{};
        size_t m_idx{};
    public:
        constexpr iter_t(iter_t const& other) : m_data(other.m_data) , m_idx(other.m_idx) {}
        constexpr iter_t(T* data, size_t idx) : m_data(data) , m_idx(idx) {}
        constexpr iter_t& operator = (iter_t const& other) { m_data = other.m_data; m_idx = other.m_idx; return *this; }
        constexpr iter_t& operator++() noexcept { ++m_idx; return *this; }
        constexpr iter_t operator + (size_t diff) { return {m_data, size_t(size_t(m_idx) + diff)}; }
        constexpr size_t operator - (iter_t const& other) { return size_t(m_idx) - size_t(other.m_idx); }
        constexpr T& operator*()  const { return m_data[m_idx >> num_bits][m_idx & mask];  }
        constexpr T* operator->() const { return &m_data[m_idx >> num_bits][m_idx & mask]; }
        constexpr bool operator == (iter_t const& o) const { return m_idx == o.m_idx; }
        constexpr bool operator != (iter_t const& o) const { return !(*this == o); }
    };

    class citer_t
    {
        const T* m_data{};
        size_t m_idx{};
    public:
        constexpr citer_t(citer_t const& other) : m_data(other.m_data) , m_idx(other.m_idx) {}
        constexpr citer_t(const T* data, size_t idx) : m_data(data) , m_idx(idx) {}
        constexpr citer_t& operator = (citer_t const& other) { m_data = other.m_data; m_idx = other.m_idx; return *this; }
        constexpr citer_t& operator++() noexcept { ++m_idx; return *this; }
        constexpr citer_t operator+(size_t diff) { return {m_data, size_t(size_t(m_idx) + diff)}; }
        constexpr size_t operator-(citer_t const& other) { return size_t(m_idx) - size_t(other.m_idx); }
        constexpr T& operator*() const { return m_data[m_idx >> num_bits][m_idx & mask]; }
        constexpr T* operator->() const { return &m_data[m_idx >> num_bits][m_idx & mask]; }
        constexpr bool operator==(citer_t const& o) const { return m_idx == o.m_idx; }
        constexpr bool operator!=(citer_t const& o) const { return !(*this == o); }
    };

    void IncreaseCapacity()
    {
        // todo use allocator here
        m_blocks.push_back(new T[num_elements_in_block]);
    }

    void AppendEveryThingFrom(SegmentedVector&& other)
    {
        Reserve(Size() + other.size());
        for (auto&& o : other) EmplaceBack(Move(o));
    }

    void AppendEveryThingFrom(SegmentedVector const& other)
    {
        Reserve(Size() + other.size());
        for (auto const& o : other) EmplaceBack(Move(o));
    }

    void Dealloc()
    {
        for (int i = 0; i < m_blocks.size(); ++i)
            delete[] m_blocks[i]; // todo use allocator.deallocate here
    }

    void Clear()
    {
        for (size_t i = 0, s = Size(); i < s; ++i)
            operator[](i).~T();
        m_size = 0;
    }

    void ShrinkToFit() {
        size_t num_blocks_required = CalcNumBlocksForCapacity(m_size);

        while (m_blocks.size() > num_blocks_required)
        {
            delete[] m_blocks.back(); // todo use allocator.deallocate here
            m_blocks.pop_back();
        }
        m_blocks.shrink_to_fit();
    }

public:

    SegmentedVector() = default;
    // segmented_vector(AllocatorT alloc) : m_blocks(vec_alloc(alloc)) {}
    // segmented_vector(segmented_vector&& other, AllocatorT alloc) : m_blocks(vec_alloc(alloc)) {
    //     if (other.get_allocator() == alloc) *this = Move(other);
    //     else // Oh my, allocator is different so we need to copy everything.
    //         AppendEverythingFrom(Move(other));
    // }
    // segmented_vector(segmented_vector const& other, Allocator alloc) : m_blocks(vec_alloc(alloc))
    // { AppendEverythingFrom(other); }

    SegmentedVector(SegmentedVector&& other)
    : m_blocks(Move(other.m_blocks)), m_size(Exchange(other.m_size, {})) {}

    SegmentedVector(SegmentedVector const& other) { AppendEverythingFrom(other); }

    SegmentedVector& operator=(SegmentedVector const& other) {
        if (this == &other) {
            return *this;
        }
        Clear();
        AppendEverythingFrom(other);
        return *this;
    }

    SegmentedVector& operator=(SegmentedVector&& other) {
        Clear();
        Cealloc();
        m_blocks = Move(other.m_blocks);
        m_size = Exchange(other.m_size, {});
        return *this;
    }

    ~SegmentedVector() {
        Clear();
        Dealloc();
    }

    [[nodiscard]] constexpr size_t Size() const { return m_size; }

    [[nodiscard]] constexpr size_t Capacity() const {
        return m_blocks.size() * num_elements_in_block;
    }

    // Indexing is highly performance critical
    constexpr auto operator[](size_t i) const noexcept -> T const& {
        return m_blocks[i >> num_bits][i & mask];
    }

    constexpr T&operator[](size_t i) noexcept {
        return m_blocks[i >> num_bits][i & mask];
    }

    constexpr const T& back() const  { return operator[](m_size - 1); }
    constexpr       T& back()        { return operator[](m_size - 1); }

    constexpr citer_t cbegin() const { return {m_blocks.data(), 0U}; }
    constexpr citer_t cend() const   { return {m_blocks.data(), m_size}; }
    constexpr citer_t begin() const  { return {m_blocks.data(), 0U}; }
    constexpr citer_t end() const    { return {m_blocks.data(), m_size}; }
    constexpr iter_t  begin()        { return {m_blocks.data(), 0U}; }
    constexpr iter_t  end()          { return {m_blocks.data(), m_size}; }

    void pop_back() { back().~T(); --m_size; }

    [[nodiscard]] bool empty() const { return 0 == m_size; }

    void reserve(size_t new_capacity) {
        m_blocks.reserve(CalcNumBlocksForCapacity(new_capacity));
        while (new_capacity > Capacity()) {
            IncreaseCapacity();
        }
    }
    // [[nodiscard]] allocator_type get_allocator() const {
    //     return allocator_type{m_blocks.get_allocator()};
    // }

    template <class... Args>
    reference emplace_back(Args&&... args)
    {
        if (m_size == Capacity()) {
            IncreaseCapacity();
        }
        void* ptr = static_cast<void*>(&operator[](m_size));
        auto& ref = *new (ptr) T(Forward<Args>(args)...);
        ++m_size;
        return ref;
    }
};
*/
