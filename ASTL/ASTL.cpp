
#include "Queue.hpp"
#include "Stack.hpp"
#include "String.hpp"
#include "Array.hpp"
#include "HashMap.hpp"
#include "HashSet.hpp"
#include "RedBlackTree.hpp"
#include "Math/Vector.hpp"
#include "String.hpp"

#include <stdio.h>
#include "Timer.hpp"

extern void AdventOfCodeTests();


int ints[13]{ 23, 456, 789, 42,
               8675309, 555, 999,
               314159, 271828, 777,
               888, 666, 99999 };

void PassingToConstFunction(HashSet<int>& testSet)
{
    bool allExist = true;
    for (int i = 0; i < 13; ++i)
    {
        testSet.Insert(ints[i]);
    }

    for (int i = 0; i < 13; ++i)
    {
        allExist &= testSet.Contains(ints[i]);
    }
    allExist &= testSet.Contains(213444590);
    ASSERT(allExist);
    printf("ours %s the test\n", allExist ? "passed" : "failed");
}


template <typename MapT> void TestMap(const char* name, MapT& map) {
    printf("\n%s BEGIN\n", name);

    CSTIMER("total speed: ");

    {
        CSTIMER("insert speed: ");
        uint64 pw = 1ull;
        for (uint64 i = 0ull; i < 90000ull; ++i) {
            map.TryEmplace(pw, "123");
            pw += 331ull;
        }
    }

    {
        CSTIMER("contains speed: ");
        uint64 pw = 1ull;
        uint64 numContains = 0ull;

        for (uint64 i = 0ull; i < 90000ull; ++i) {
            numContains += map.Find(pw) != map.cend();
            pw += 331ull;
        }
        printf("numContains: %i\n", numContains);
    }

    {
        CSTIMER("Iteration speed: ");
        uint64 keySum = 0ull;
        uint64 valueSum = 0ull;

        for (const auto& val : map) {
            keySum += val.key;
            valueSum += val.value.Length();
        }
        printf("key sum: %lu, value sum: %lu\n", keySum, valueSum);
    }

    {
        CSTIMER("Erase speed: ");
        uint64 pw = 1ull;
        for (uint64 i = 0ull; i < 90000ull; ++i) {
            map.Erase((int)pw);
            pw += 331ull;
        }
    }

    printf("%s END\n", name);
}

int main()
{
    AdventOfCodeTests();
    String testStr = "floating test: ";
    testStr += 1234.567f;
    testStr.Replace("floating", "integing");
    testStr.Replace("1234", "43210");

    printf("%s\n", testStr.CStr());
    
    String testStr1 = "int test: ";
    testStr1 += 1234567;
    printf("%s\n", testStr1.CStr());
    

    static const KeyValuePair<int, String> testMapInitializer[13] =
    {
        { 23     , "apple"    },
        { 456    , "banana"   },
        { 789    , "cherry"   },
        { 42     , "dog"      },
        { 8675309, "elephant" },
        { 555    , "frog"     },
        { 999    , "guitar"   },
        { 314159 , "house"    },
        { 271828 , "igloo"    },
        { 777    , "jigsaw"   },
        { 888    , "kangaroo" },
        { 666    , "lemon"    },
        { 99999  , "monkey"   }
    };

    static uint64_t arraySize = ArraySize(testMapInitializer);
    HashMap<int, String> testMap(testMapInitializer, arraySize);

    testMap[4] = "no its five";

    String strings[13] = {
        "apple", "banana", "cherry", "dog", "elephant", "frog", "guitar",
        "house", "igloo", "jigsaw", "kangaroo", "lemon", "monkey"
    };

    Array<String> strs(strings, strings+13);

    strs.RemoveAt(3, 3);

    for (int i = 0; i < strs.Size(); ++i)
    {
        printf("%s \n", strs[i].c_str());
    }

    for (const KeyValuePair<int, String>& x : testMap)
    {
        printf("%i, %s\n", x.key, x.value.CStr());
    }
    ASSERT(testMap[99999] == String("monkey"));

    HashSet<int> testSet{};
    testSet.Insert(213444590);

    PassingToConstFunction(testSet);

    HashMap<int, String> ourMap = testMap;

    TestMap("ourMap", ourMap);

    PassingToConstFunction(testSet);

    String num = ourMap[4];
    ASSERT(num == String("no its five"));
    bool finded1 = ourMap.Find(5ull) != ourMap.end();
    ASSERT(finded1 == false);

    String val1 = ourMap.At(4);
    ASSERT(val1 == String("no its five"));

    ourMap.TryEmplace(33, "33ull");

    bool contains1 = ourMap.Contains(33);
    ASSERT(contains1);
    ourMap.Erase(33);

    contains1 = ourMap.Contains(33ull);
    ASSERT(!contains1);

    ourMap.Insert(33, "33ull");
    contains1 = ourMap.Contains(33ull);
    ASSERT(contains1);

    return 0;
}
  
#ifdef __NOTHING

/* Genetic algorithm to explore xorshift-multiply-xorshift hashes.
 */
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "Math/SIMDCommon.hpp"

#define POOL      40
#define THRESHOLD 2.0  // Use exact when estimate is below this
#define DONTCARE  4.0  // Only print tuples with bias below this threshold
#define QUALITY   18   // 2^N iterations of estimate samples
#define RESETMINS 90   // Reset pool after this many minutes of no progress

static uint64_t
rand64(uint64_t s[4])
{
    uint64_t x = s[1] * 5;
    uint64_t r = ((x << 7) | (x >> 57)) * 9;
    uint64_t t = s[1] << 17;
    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = (s[3] << 45) | (s[3] >> 19);
    return r;
}


#define FLAG_SCORED  (1u << 0)
#define FLAG_EXACT   (1u << 1)
#define FLAG_PRINTED (1u << 2)

struct gene {
    double score;
    short s[3];
    uint32_t c[2];
    unsigned flags;
};

static uint32_t hash(const struct gene* g, uint32_t x)
{
    x ^= x >> g->s[0];
    x *= g->c[0];
    x ^= x >> g->s[1];
    x *= g->c[1];
    x ^= x >> g->s[2];
    return x;
}

static __m256i SIMDhash(const struct gene* g, __m256i x)
{
    x = _mm256_xor_si256(x, _mm256_slli_epi32(x, g->s[0]));
    x = _mm256_mul_epi32(x, _mm256_set1_epi32(g->c[0]));
    x = _mm256_xor_si256(x, _mm256_slli_epi32(x, g->s[1]));
    x = _mm256_mul_epi32(x, _mm256_set1_epi32(g->c[1]));
    x = _mm256_xor_si256(x, _mm256_slli_epi32(x, g->s[2]));
    return x;
}

#include "Math/SIMDCommon.hpp"

#include <string.h>

static double estimate_bias32SIMD(const struct gene* g, uint64_t rng[4])
{
    const long n = 1L << QUALITY;
    static const __m256i m1   = _mm256_set1_epi32(0xed5ad4bbu), 
                         m2   = _mm256_set1_epi32(0xac4c1b51u),
                         m3   = _mm256_set1_epi32(0x31848babu),
                         zero  = _mm256_setzero_si256(),
                         one   = _mm256_set1_epi32(1);

    __m256i bins[8][8];
    memset(bins, 0, sizeof(__m256i) * 8 * 8);
    __m256i x = _mm256_loadu_epi64(rng);

    for (long i = 0; i < n / 8; ++i) {
        x = _mm256_xor_si256(x, _mm256_slli_epi32(x, 17));
        x = _mm256_mul_epi32(x, m1);
        x = _mm256_xor_si256(x, _mm256_slli_epi32(x, 11));
        x = _mm256_mul_epi32(x, m2);
        x = _mm256_xor_si256(x, _mm256_slli_epi32(x, 15));
        x = _mm256_mul_epi32(x, m3);
        x = _mm256_xor_si256(x, _mm256_slli_epi32(x, 14));
        x = _mm256_add_epi32(x, _mm256_set_epi32(1, 7, 3, 5, 9, 3, 5, 7));
        __m256i h0 = SIMDhash(g, x);

        for (int j = 0; j < 8; j++) {
            __m256i bit = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
            __m256i h1 = SIMDhash(g, _mm256_xor_si256(x, bit));
            __m256i set = _mm256_xor_si256(h0, h1);
            __m256i k = _mm256_and_epi32(_mm256_srlv_epi32(set, bit), one);
            
            bins[j][0] = _mm256_add_epi64(bins[j][0], _mm256_unpacklo_epi32(k, zero));
            bins[j][1] = _mm256_add_epi64(bins[j][1], _mm256_unpackhi_epi32(k, zero));

            bit = _mm256_set_epi32(8, 9, 10, 11, 12, 13, 14, 15);
            k = _mm256_and_epi32(_mm256_srlv_epi32(set, bit), one);
            bins[j][2] = _mm256_add_epi64(bins[j][2], _mm256_unpacklo_epi32(k, zero));
            bins[j][3] = _mm256_add_epi64(bins[j][3], _mm256_unpackhi_epi32(k, zero));
          
            bit = _mm256_set_epi32(16, 17, 18, 19, 20, 21, 22, 23);
            k = _mm256_and_epi32(_mm256_srlv_epi32(set, bit), one);
            bins[j][4] = _mm256_add_epi64(bins[j][4], _mm256_unpacklo_epi32(k, zero));
            bins[j][5] = _mm256_add_epi64(bins[j][5], _mm256_unpackhi_epi32(k, zero));
          
            bit = _mm256_set_epi32(24, 25, 26, 27, 28, 29, 30, 31);
            k = _mm256_and_epi32(_mm256_srlv_epi32(set, bit), one);
            bins[j][6] = _mm256_add_epi64(bins[j][6], _mm256_unpacklo_epi32(k, zero));
            bins[j][7] = _mm256_add_epi64(bins[j][7], _mm256_unpackhi_epi32(k, zero));
        }
    }
    _mm256_storeu_epi64(rng, x);
    
    static const __m256d t1024 = _mm256_set1_pd(1024.0), two = _mm256_set1_pd(2.0);
    static const __m256i n4 = _mm256_set1_epi64x(n);

    __m256d mean = _mm256_setzero_pd(),
           ndiv2 = _mm256_div_pd(_mm256_castsi256_pd(n4), two);

    for (int j = 0; j < 8; ++j) { 
        for (int i = 0; i < 8; ++i) {
            __m256d d = _mm256_div_pd(_mm256_castsi256_pd(_mm256_sub_epi64(bins[j][i], n4)), ndiv2);
            mean = _mm256_div_pd(_mm256_mul_pd(d, d), t1024);
        }
    }
    mean = _mm256_add_pd(mean, _mm256_permute2f128_pd(mean, mean, 1));
    mean = _mm256_hadd_pd(mean, mean);
    return _mm256_cvtsd_f64(_mm256_hadd_pd(mean, mean));
}

static double
estimate_bias32(const struct gene* g, uint64_t rng[4])
{
    long n = 1L << QUALITY;
    long bins[32][32] = { {0} };
    for (long i = 0; i < n; i++) {
        uint32_t x = rand64(rng);
        uint32_t h0 = hash(g, x);
        for (int j = 0; j < 32; j++) {
            uint32_t bit = UINT32_C(1) << j;
            uint32_t h1 = hash(g, x ^ bit);
            uint32_t set = h0 ^ h1;
            for (int k = 0; k < 32; k++)
                bins[j][k] += (set >> k) & 1;
        }
    }
    double mean = 0;
    for (int j = 0; j < 32; j++) {
        for (int k = 0; k < 32; k++) {
            double diff = (bins[j][k] - n / 2) / (n / 2.0);
            mean += (diff * diff) / (32.0 * 32.0);
        }
    }
    return sqrt(mean) * 1000.0;
}

#define EXACT_SPLIT 32  // must be power of two
static double
exact_bias32(const struct gene* g)
{
    long long bins[32][32] = { {0} };
    static const uint64_t range = (UINT64_C(1) << 32) / EXACT_SPLIT;
#pragma omp parallel for
    for (int i = 0; i < EXACT_SPLIT; i++) {
        long long b[32][32] = { {0} };
        for (uint64_t x = i * range; x < (i + 1) * range; x++) {
            uint32_t h0 = hash(g, x);
            for (int j = 0; j < 32; j++) {
                uint32_t bit = UINT32_C(1) << j;
                uint32_t h1 = hash(g, x ^ bit);
                uint32_t set = h0 ^ h1;
                for (int k = 0; k < 32; k++)
                    b[j][k] += (set >> k) & 1;
            }
        }
#pragma omp critical
        for (int j = 0; j < 32; j++)
            for (int k = 0; k < 32; k++)
                bins[j][k] += b[j][k];
    }
    double mean = 0.0;
    for (int j = 0; j < 32; j++) {
        for (int k = 0; k < 32; k++) {
            double diff = (bins[j][k] - 2147483648L) / 2147483648.0;
            mean += (diff * diff) / (32 * 32);
        }
    }

    return sqrt(mean) * 1000.0;
}

static double
SIMDexact_bias32(const struct gene* g)
{
    __m256i vbins[8][8];
    {
        long long bins[32][32] = { {0} };
        static const uint64_t range = (UINT64_C(1) << 32) / EXACT_SPLIT;
        #pragma omp parallel for
        for (int i = 0; i < EXACT_SPLIT; i++) {
            long long b[32][32] = { {0} };
            for (uint64_t x = i * range; x < (i + 1) * range; x++) {
                uint32_t h0 = hash(g, x);
                for (int j = 0; j < 32; j++) {
                    uint32_t bit = UINT32_C(1) << j;
                    uint32_t h1 = hash(g, x ^ bit);
                    uint32_t set = h0 ^ h1;
                    for (int k = 0; k < 32; k++)
                        b[j][k] += (set >> k) & 1;
                }
            }
            __m256i* vb    = (__m256i*)b;
            __m256i* vbins = (__m256i*)bins;
            // todo prefetch a
            for (int j = 0; j < 64; j += 4)
            {
                vbins[j + 0] = _mm256_add_epi64(vbins[j + 0], vb[j + 0]);
                vbins[j + 1] = _mm256_add_epi64(vbins[j + 1], vb[j + 1]);
                vbins[j + 2] = _mm256_add_epi64(vbins[j + 2], vb[j + 2]);
                vbins[j + 3] = _mm256_add_epi64(vbins[j + 3], vb[j + 3]);
            }
        }
        __m256i* xbins = (__m256i*)bins;
        for (int j = 0; j < 64; ++j)
        {
            vbins[0][j] = *xbins++;
        }
    }
    
    static const __m256d t1024 = _mm256_set1_pd(1024.0);
    static const __m256i two = _mm256_set1_epi32(2);

    __m256i ml = _mm256_set1_epi64x(2147483648L);
    __m256d mean = _mm256_setzero_pd(), md = _mm256_set1_pd(2147483648.0);

    for (int j = 0; j < 8; ++j) { 
        for (int i = 0; i < 8; ++i) {
            __m256d d = _mm256_div_pd(_mm256_castsi256_pd(_mm256_sub_epi64(vbins[j][i], ml)), md);
            d =  _mm256_div_pd(_mm256_mul_pd(d, d), t1024);
            mean = _mm256_add_pd(mean, d);
        }
    }
    mean = _mm256_add_pd(mean, _mm256_permute2f128_pd(mean, mean, 1));
    mean = _mm256_hadd_pd(mean, mean);
    return sqrt(_mm256_cvtsd_f64(_mm256_hadd_pd(mean, mean))) * 1000.0;
}

static void
gene_gen(struct gene* g, uint64_t rng[4])
{
    uint64_t s = rand64(rng);
    uint64_t c = rand64(rng);
    g->s[0] = 10 + (s >> 0) % 10;
    g->s[1] = 10 + (s >> 24) % 10;
    g->s[2] = 10 + (s >> 48) % 10;
    g->c[0] = c | 1u;
    g->c[1] = (c >> 32) | 1u;
    g->flags = 0;
}

static void
gene_print(const struct gene* g, FILE* f)
{
    fprintf(f, "[%2d %08lx %2d %08lx %2d]",
        g->s[0], (unsigned long)g->c[0],
        g->s[1], (unsigned long)g->c[1], g->s[2]);
}

static int
small(uint64_t r)
{
    static const int v[] = { -3, -2, -1, +1, +2, +3 };
    return v[r % 6];
}

static void
gene_mutate(struct gene* g, uint64_t rng[4])
{
    uint64_t r = rand64(rng);
    int s = r % 5;
    r >>= 3;
    switch (s) {
    case 0:
        g->s[0] += small(r);
        break;
    case 1:
        g->s[1] += small(r);
        break;
    case 2:
        g->s[2] += small(r);
        break;
    case 3:
        g->c[0] += (int)(r & 0xffff) - 32768;
        break;
    case 4:
        g->c[1] += (int)(r & 0xffff) - 32768;
        break;
    }
    g->flags = 0;
}

static void
gene_cross(struct gene* g,
    const struct gene* a,
    const struct gene* b,
    uint64_t rng[4])
{
    uint64_t r = rand64(rng);
    *g = *a;
    switch (r & 2) {
    case 0: g->c[0] = b->c[0]; /* FALLTHROUGH */
    case 1: g->s[1] = b->s[1]; /* FALLTHROUGH */
    case 2: g->c[1] = b->c[1]; /* FALLTHROUGH */
    case 3: g->s[2] = b->s[2];
    }
    g->flags = 0;
}

static int
gene_same(const struct gene* a, const struct gene* b)
{
    return a->s[0] == b->s[0] &&
        a->s[1] == b->s[1] &&
        a->s[2] == b->s[2] &&
        a->c[0] == b->c[0] &&
        a->c[1] == b->c[1];
}

#include "Random.hpp"

static void
rng_init(void* p, size_t len)
{
    using namespace Random;
    uint64_t xoro[2];
    Xoroshiro128PlusInit(xoro);
    uint64_t* ip = (uint64_t*)p;

    for (uint64_t i = 0; i < len / sizeof(uint64_t); i++)
    {
        *ip++ = Xoroshiro128Plus(xoro);
    }
}

static int
cmp(const void* pa, const void* pb)
{
    double a = *(double*)pa;
    double b = *(double*)pb;
    if (a < b)
        return -1;
    if (b < a)
        return 1;
    return 0;
}

static void
undup(struct gene* pool, uint64_t rng[4])
{
    for (int i = 0; i < POOL; i++)
        for (int j = i + 1; j < POOL; j++)
            if (gene_same(pool + i, pool + j))
                gene_mutate(pool + j, rng);
}

#include "Timer.hpp"

int
main(void)
{
    int verbose = 1;
    double best = 1000.0;
    time_t best_time = time(0);
    uint64_t rng[POOL][4];
    struct gene pool[POOL];

    rng_init(rng, sizeof(rng));
    for (int i = 0; i < POOL; i++)
        gene_gen(pool + i, rng[0]);
    // 2719.915000
    int q = 3;
    CSTIMER("total time: ");
    while (q-- > 0) {
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < POOL; i++) {
            if (!(pool[i].flags & FLAG_SCORED)) {
                pool[i].score = estimate_bias32SIMD(pool + i, rng[i]);
                pool[i].flags |= FLAG_SCORED;
            }
        }
        for (int i = 0; i < POOL; i++) {
            if (!(pool[i].flags & FLAG_EXACT) && pool[i].score < THRESHOLD) {
                pool[i].score = exact_bias32(pool + i);
                pool[i].flags |= FLAG_EXACT;
            }
        }

        qsort(pool, POOL, sizeof(*pool), cmp);
        if (verbose) {
            for (int i = 0; i < POOL; i++) {
                if (!(pool[i].flags & FLAG_PRINTED) &&
                    pool[i].score < DONTCARE) {
                    gene_print(pool + i, stdout);
                    printf(" = %.17g\n", pool[i].score);
                    pool[i].flags |= FLAG_PRINTED;
                }
            }
        }

        time_t now = time(0);
        if (pool[0].score < best) {
            best = pool[0].score;
            best_time = now;
        }
        else if (now - best_time > RESETMINS * 60) {
            best = 1000.0;
            best_time = now;
            for (int i = 0; i < POOL; i++)
                gene_gen(pool + i, rng[0]);
        }

        int c = POOL / 4;
        for (int a = 0; c < POOL && a < POOL / 4; a++)
            for (int b = a + 1; c < POOL && b < POOL / 4; b++)
                gene_cross(pool + c++, pool + a, pool + b, rng[0]);
        undup(pool, rng[0]);
    }
}
#endif
