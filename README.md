# ASTL

Work in progress. This library is composition of lots of paper, research and development <br>
For compatibility I've used lots of macros in Common.hpp but still we have platform dependent code <br>
Also I've avoided most of the modern C++ features <br>
feel free to contribute, use, modify or oppening PR<br>
I havent used any C++ headers except thread and atomic.<br>
Goal is make this library compile as fast as possible,  easy to use. and easy to read <br>
We have Hash functions and Random number generators in Random.hpp. <br>
Matrix4 and Vector4 uses SIMD extensions(SSE3)<br>
math library is combination of glm and DirectX Math. examples:
```cpp
static float f = 1.0f; f += 0.01f;
constexpr float distance = 3.14159265f; // this is distance from cube but I did use pi anyways 
Vector3f position(Sin(f) * distance, 0.0f, Cos(f) * distance );
float verticalFOV = 65.0f, nearClip = 0.01f, farClip = 500.0f;

Matrix4 projection = Matrix4::PerspectiveFovRH(verticalFOV * DegToRad, m_NativeWindow->GetWidth(), m_NativeWindow->GetHeight(), nearClip, farClip);
Matrix4 viewProjection = projection * Matrix4::LookAtRH(position, -Vector3f::Normalize(position), Vector3f::Up());

Vector3f baryCentrics = Vector3f(1.0f - hitOut.u - hitOut.v, hitOut.u, hitOut.v);
Matrix3 inverseMat3 = Matrix4::ConvertToMatrix3(hitInstance.inverseTransform);
		
Vector3f n0 = Matrix3::Multiply(inverseMat3, ConvertToFloat3(&triangle.normal0x));
Vector3f n1 = Matrix3::Multiply(inverseMat3, ConvertToFloat3(&triangle.normal1x));
Vector3f n2 = Matrix3::Multiply(inverseMat3, ConvertToFloat3(&triangle.normal2x));
	        
record.normal = Vector3f::Normalize((n0 * baryCentrics.x) + (n1 * baryCentrics.y) + (n2 * baryCentrics.z));

Vector2i x = Vector2i(1, 1) + Vector2(2, 2);
float length =  ToVector2f(x).Length();
```
<br>
HashMap uses Ankerl's algorithm (way faster than std::unordered_map and uses contigues memory)<br>
Todo: RedBlackTree, ScopedPtr, and SharedPtr. coming soon.<br>
I have Array<T> instead of std::vector<T> <br>
also you can use Bitset1024, Bitset512, Bitset256(SIMD optimized) or Bitset<1234> instead of std::bitset <br>
Queue and Stack will be fixed soon. Here is examples

```cpp
String testStr = "floating test: ";
testStr += 1234.567f;
printf("%s\n", testStr.c_str());
String testStr1 = "int test: ";
testStr1 += 1234567;
printf("%s\n", testStr1.c_str());
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
static HashMap<int, String> testMap(testMapInitializer, 13);
        
testMap[4] = "no its five";
String strings[13] = {
    "apple", "banana", "cherry", "dog", "elephant", "frog", "guitar",
    "house", "igloo", "jigsaw", "kangaroo", "lemon", "monkey"
};       

Array<String> strs{};

for (int i = 0; i < 13; ++i)
{
    strs.Add(strings[i]);
}

strs.RemoveAt(3, 3);
for (int i = 0; i < strs.Size(); ++i)
{
    printf("%s \n", strs[i].c_str());
}

for (const auto&[key, value] : testMap)
{
    printf("%i, %s\n", key, value.CStr());
}

HashSet<int> testSet{};
testSet.insert(213444590);
HashMap<int, String> ourMap = testMap;
TestMap("ourMap", ourMap);
String num = ourMap[4];

bool finded1 = ourMap.Find(5ull) != ourMap.end();
String val1 = ourMap.At(4);

ourMap.TryEmplace(33, "33ull");
bool contains1 = ourMap.Contains(33);

ourMap.Erase(33);
contains1 = ourMap.Contains(33ull);

ourMap.Insert(33, "33ull");
contains1 = ourMap.Contains(33ull);
```
also you can use BVH algorithm(Jacco Bikker's but optimized with SIMD) for Ray casting to the 3d scene<br>
it will return mesh index, hit position hit color etc.<br>
note: I haven't touched multi threading code for a while