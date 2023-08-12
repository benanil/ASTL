
#include "../Queue.hpp"
#include "../Stack.hpp"
#include "../String.hpp"
#include "../Array.hpp"
#include "../HashMap.hpp"
#include "../HashSet.hpp"
#include "../RedBlackTree.hpp"
#include "../Math/Vector.hpp"
#include "../String.hpp"

#include "../Profiler.hpp"

#include <stdio.h>

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

    // {
    //     TimeBlock("Testmap HashMap insert speed: ");
    //     uint64_t pw = 1ull;
    //     for (uint64_t i = 0ull; i < 90000ull; ++i) {
    //         map.TryEmplace(pw, "123");
    //         pw += 331ull;
    //     }
    // }
    // 
    // {
    //     TimeBlock("Testmap HashMap contains speed: ");
    //     uint64_t pw = 1ull;
    //     uint64_t numContains = 0ull;
    // 
    //     for (uint64_t i = 0ull; i < 90000ull; ++i) {
    //         numContains += map.Find(pw) != map.cend();
    //         pw += 331ull;
    //     }
    //     printf("numContains: %i\n", numContains);
    // }
    // 
    // {
    //     TimeBlock("Testmap HashMap Iteration speed: ");
    //     uint64_t keySum = 0ull;
    //     uint64_t valueSum = 0ull;
    // 
    //     for (const auto& val : map) {
    //         keySum += val.key;
    //         valueSum += val.value.Length();
    //     }
    //     printf("key sum: %lu, value sum: %lu\n", keySum, valueSum);
    // }
    // 
    // {
    //     TimeBlock("Testmap HashMap Erase speed: ");
    //     uint64_t pw = 1ull;
    //     for (uint64_t i = 0ull; i < 90000ull; ++i) {
    //         map.Erase((int)pw);
    //         pw += 331ull;
    //     }
    // }

    String num = map[4];
    ASSERT(num == String("no its five"));
    bool finded1 = map.Find(5ull) != map.end();
    ASSERT(finded1 == false);

    // String val1 = map.At(4);
    // ASSERT(val1 == String("no its five"));
    // 
    // map.TryEmplace(33, "33ull");
    // 
    // bool contains1 = map.Contains(33);
    // ASSERT(contains1);
    // map.Erase(33);
    // 
    // contains1 = map.Contains(33ull);
    // ASSERT(!contains1);
    // 
    // map.Insert(33, "33ull");
    // contains1 = map.Contains(33ull);
    // ASSERT(contains1);

    printf("%s END\n", name);
}

#include "../Aditional.hpp"

void StrViewTest(StringView view)
{
  printf("str view test: %s \n", view.ptr);
}

extern void BeginProfile(void);
extern void EndAndPrintProfile();

#include "../Math/Matrix.hpp"

int main()
{   
    AdventOfCodeTests();

    Matrix4 matrix = Matrix4::Identity();
    matrix = matrix * matrix;

    StackHashMap<int, float, 32> stackHashMap;
    float flt = 0.5f;
    Random::PCG pcg;
    Random::PCGInitialize(pcg, 1254567);
    stackHashMap[Random::PCGNext(pcg)] = flt += 1.0f;
    stackHashMap[Random::PCGNext(pcg)] = flt += 1.0f;
    stackHashMap[Random::PCGNext(pcg)] = flt += 1.0f;
    stackHashMap[Random::PCGNext(pcg)] = flt += 1.0f;
    stackHashMap[Random::PCGNext(pcg)] = flt += 1.0f;
    stackHashMap[Random::PCGNext(pcg)] = flt += 1.0f;

    StackArray<int, 5> stackArray;
    stackArray.Add(1);
    stackArray.Add(2);
    stackArray.Add(3);
    stackArray.Add(4);
    stackArray.Add(5);

    for (const KeyValuePair<int, float>& x : stackHashMap)
    {
        printf("stack test: %u, %f \n", x.key, x.value);
    }

    String testStr = "floating test: ";
    testStr += 1234.567f;
    testStr.Replace("floating", "integing");
    testStr.Replace("1234", "43210");
    testStr.Insert(testStr.Length()-2, "long string verry verry long it is huge that I can test it with this string so it will pass and works fine yeah just for test you now I'm sayin");

    printf("%s\n", testStr.CStr());
    
    String testStr1 = "int test: ";
    testStr1 += 1234567;
    StrViewTest(testStr1);
    
    GrowablePoolAllocator<uint64_t > allocator(10);
    
    for (int i = 0; i < 9; i++)
        *allocator.Allocate(1) = i + '0';

    uint64_t * walk = allocator.firstBlock->data;

    for (int i = 0; i < 9; i++)
        printf("walk[i]: %c  ", (char)*walk++);

    allocator.Allocate(1);

    
    BeginProfile();
    {
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

        static uint64_t  arraySize = ArraySize(testMapInitializer);
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

        String doozyuzdoksandooz = testMap[99999];
    
        ASSERT(testMap[99999] == String("monkey"));

        HashMap<int, String> ourMap = testMap;
        TestMap("ourMap", ourMap);
    } // begin profile
    EndAndPrintProfile();

    HashSet<int> testSet{};
    testSet.Insert(213444590);

    PassingToConstFunction(testSet);
    
    getchar();

    return 0;
}
  