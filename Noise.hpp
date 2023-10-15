
#include "Math/Vector.hpp"
#include "Random.hpp"

// todo make simd

// simplex noise in 2D/3D/4D, source: Stefan Gustavson, https://weber.itn.liu.se/~stegu/simplexnoise/SimplexNoise.java
// https://github.com/ProjectPhysX/FluidX3D/blob/master/src/utilities.hpp#L2417
struct SimplexNoise 
{
    struct float4_32 
    {
        float x, y, z, w;
 
        float4_32(float x, float y, float z, float w)
        {
            this->x = x; this->y = y; this->z = z; this->w = w;
        }

        float4_32(int x, int y, int z, int w = 0) 
        {
            this->x = (float)x; this->y = (float)y; this->z = (float)z; this->w = (float)w;
        }
    };

    static const float4_32 grad3[12] = {
        { 1.0f, 1.0f, 0.0f }, {-1.0f, 1.0f, 0.0f}, {1.0f,-1.0f, 0.0f}, {-1.0f,-1.0f, 0.0f},
        { 1.0f, 0.0f, 1.0f }, {-1.0f, 0.0f, 1.0f}, {1.0f, 0.0f,-1.0f}, {-1.0f, 0.0f,-1.0f},
        { 0.0f, 1.0f, 1.0f }, { 0.0f,-1.0f, 1.0f}, {0.0f, 1.0f,-1.0f}, { 0.0f,-1.0f,-1.0f}
    };

    static const float4_32 grad4[32] = {
        { 0.0f, 1.0f, 1.0f, 1.0f}, { 0.0f, 1.0f, 1.0f, -1.0f}, { 0.0f, 1.0f,-1.0f, 1.0f},  0.0f, 1.0f,-1.0f,-1.0f},
        { 0.0f,-1.0f, 1.0f, 1.0f}, { 0.0f,-1.0f, 1.0f, -1.0f}, { 0.0f,-1.0f,-1.0f, 1.0f},  0.0f,-1.0f,-1.0f,-1.0f},
        { 1.0f, 0.0f, 1.0f, 1.0f}, { 1.0f, 0.0f, 1.0f, -1.0f}, { 1.0f, 0.0f,-1.0f, 1.0f},  1.0f, 0.0f,-1.0f,-1.0f},
        {-1.0f, 0.0f, 1.0f, 1.0f}, {-1.0f, 0.0f, 1.0f, -1.0f}, {-1.0f, 0.0f,-1.0f, 1.0f}, -1.0f, 0.0f,-1.0f,-1.0f},
        { 1.0f, 1.0f, 0.0f, 1.0f}, { 1.0f, 1.0f, 0.0f, -1.0f}, { 1.0f,-1.0f, 0.0f, 1.0f},  1.0f,-1.0f, 0.0f,-1.0f},
        {-1.0f, 1.0f, 0.0f, 1.0f}, {-1.0f, 1.0f, 0.0f, -1.0f}, {-1.0f,-1.0f, 0.0f, 1.0f}, -1.0f,-1.0f, 0.0f,-1.0f},
        { 1.0f, 1.0f, 1.0f, 0.0f}, { 1.0f, 1.0f,-1.0f,  0.0f}, { 1.0f,-1.0f, 1.0f, 0.0f},  1.0f,-1.0f,-1.0f, 0.0f},
        {-1.0f, 1.0f, 1.0f, 0.0f}, {-1.0f, 1.0f,-1.0f,  0.0f}, {-1.0f,-1.0f, 1.0f, 0.0f}, -1.0f,-1.0f,-1.0f, 0.0f}
    };
   
    const float F2 = 0.5f*(1.73205080757f-1.0f), G2=(3.0f-1.73205080757f)/6.0f; // 1.73205080757 is sqrt3
    const float F3 = 1.0f/3.0f, G3=1.0f/6.0f;
    const float F4 = (2.2360679775f-1.0f)*0.25f, G4=(5.0f-2.2360679775f)*0.05f; // 2.2360679775f is sqrt5
    uchar perm[512];
    uchar perm12[512];

    int floor(const float x) const { return (int)x-(x<=0.0f); }
    float dot(const float4_32& g, const float x, const float y) const { return g.x*x+g.y*y; }
    float dot(const float4_32& g, const float x, const float y, const float z) const { return g.x*x+g.y*y+g.z*z; }
    float dot(const float4_32& g, const float x, const float y, const float z, const float w) const { return g.x*x+g.y*y+g.z*z+g.w*w; }
    
    SimplexNoise() 
    {
        const uchar p[256] = {
            151,160,137, 91, 90, 15,131, 13,201, 95, 96, 53,194,233,  7,225,140, 36,103, 30, 69,142,  8, 99, 37,240, 21, 10, 23,190,  6,148,
            247,120,234, 75,  0, 26,197, 62, 94,252,219,203,117, 35, 11, 32, 57,177, 33, 88,237,149, 56, 87,174, 20,125,136,171,168, 68,175,
            74,165, 71,134,139, 48, 27,166, 77,146,158,231, 83,111,229,122, 60,211,133,230,220,105, 92, 41, 55, 46,245, 40,244,102,143, 54,
            65, 25, 63,161,  1,216, 80, 73,209, 76,132,187,208, 89, 18,169,200,196,135,130,116,188,159, 86,164,100,109,198,173,186,  3, 64,
            52,217,226,250,124,123,  5,202, 38,147,118,126,255, 82, 85,212,207,206, 59,227, 47, 16, 58, 17,182,189, 28, 42,223,183,170,213,
            119,248,152,  2, 44,154,163, 70,221,153,101,155,167, 43,172,  9,129, 22, 39,253, 19, 98,108,110,79,113,224,232,178,185, 112,104,
            218,246, 97,228,251, 34,242,193,238,210,144, 12,191,179,162,241, 81, 51,145,235,249, 14,239,107, 49,192,214, 31,181,199,106,157,
            184, 84,204,176,115,121, 50, 45,127,  4,150,254,138,236,205, 93,222,114, 67, 29, 24, 72,243,141,128,195, 78, 66,215, 61,156,180
        };
        for(int i = 0; i < 512; i++) 
        {
            perm[i] = p[i & 255];
            perm12[i] = (uchar)(perm[i] % 12);
        }
    }
    
    // 2D simplex noise
    float noise(float x, float y) const 
    { 
        float n0, n1, n2;
        float s  = (x+y)*F2;
        int   i  = Floor(x+s), j = Floor(y+s);
        float t  = (i+j)*G2;
        float X0 = i - t , Y0  = j - t;
        float x0 = x - X0, y0 = y - Y0;
        int i1, j1;
        
        if (x0 > y0) { i1 = 1; j1 = 0; }
        else { i1 = 0; j1 = 1; }

        float x1 = x0 -   i1 +        G2, y1 = y0 -   j1 +        G2;
        float x2 = x0 - 1.0f + 2.0f * G2, y2 = y0 - 1.0f + 2.0f * G2;
        int ii = i & 255, jj = j & 255;
        int gi0 = perm12[ii      + perm[jj     ]];
        int gi1 = perm12[ii + i1 + perm[jj + j1]];
        int gi2 = perm12[ii + 1  + perm[jj + 1]];
        
        float t0 = 0.5f - x0 * x0 - y0 * y0;
        if(t0 < 0.0f) n0 = 0.0f; else { t0 *= t0; n0 = t0*t0*dot(grad3[gi0], x0, y0); }
        
        float t1 = 0.5f-x1*x1-y1*y1;
        if(t1 < 0.0f) n1 = 0.0f; else { t1 *= t1; n1 = t1*t1*dot(grad3[gi1], x1, y1); }
        
        float t2 = 0.5f - x2 * x2 -y2 * y2;
        if(t2 < 0.0f) n2 = 0.0f; else { t2 *= t2; n2 = t2*t2*dot(grad3[gi2], x2, y2); }
        
        return 70.0f * (n0 + n1 + n2);
    }

    // 3D simplex noise
    float noise(float x, float y, float z) const 
    { 
        float n0, n1, n2, n3;
        float s = (x + y + z) * F3;
        int i = Floor(x+s), j = Floor(y+s), k = Floor(z+s);
        float t = (i+j+k)*G3;
        float X0=i-t, Y0=j-t, Z0=k-t;
        float x0=x-X0, y0=y-Y0, z0=z-Z0;
        int i1, j1, k1, i2, j2, k2;
        if(x0 >= y0) {
            if (y0 >= z0)   { i1=1; j1=0; k1=0; i2=1; j2=1; k2=0; }
            else if(x0>=z0) { i1=1; j1=0; k1=0; i2=1; j2=0; k2=1; }
            else            { i1=0; j1=0; k1=1; i2=1; j2=0; k2=1; }
        } else {
            if (y0 < z0)     { i1=0; j1=0; k1=1; i2=0; j2=1; k2=1; }
            else if(x0 < z0) { i1=0; j1=1; k1=0; i2=0; j2=1; k2=1; }
            else             { i1=0; j1=1; k1=0; i2=1; j2=1; k2=0; }
        }
        
        float x1 = x0 -   i1 +        G3, y1 = y0 -   j1 +        G3, z1 = z0 -   k1 +        G3;
        float x2 = x0 -   i2 + 2.0f * G3, y2 = y0 -   j2 + 2.0f * G3, z2 = z0 -   k2 + 2.0f * G3;
        float x3 = x0 - 1.0f + 3.0f * G3, y3 = y0 - 1.0f + 3.0f * G3, z3 = z0 - 1.0f + 3.0f * G3;
        
        int ii = i & 255, jj = j & 255, kk = k & 255;
        int gi0 = perm12[ii      + perm[jj      + perm[kk     ]]];
        int gi1 = perm12[ii + i1 + perm[jj + j1 + perm[kk + k1]]];
        int gi2 = perm12[ii + i2 + perm[jj + j2 + perm[kk + k2]]];
        int gi3 = perm12[ii +  1 + perm[jj +  1 + perm[kk +  1]]];
        
        float t0 = 0.6f - x0 * x0 - y0 * y0 - z0 * z0;
        if(t0 < 0.0f) n0 = 0.0f; else { t0 *= t0; n0 = t0 * t0 * dot(grad3[gi0], x0, y0, z0); }
        
        float t1 = 0.6f - x1 * x1 - y1 * y1 - z1 * z1;
        if(t1 < 0.0f) n1 = 0.0f; else { t1 *= t1; n1 = t1 * t1 * dot(grad3[gi1], x1, y1, z1); }
        
        float t2 = 0.6f - x2 * x2 - y2 * y2 - z2 * z2;
        if(t2 < 0.0f) n2 = 0.0f; else { t2 *= t2; n2 = t2 * t2 * dot(grad3[gi2], x2, y2, z2); }
        
        float t3 = 0.6f - x3 * x3 - y3 * y3 - z3 * z3;
        if(t3 < 0.0f) n3 = 0.0f; else { t3 *= t3; n3 = t3 * t3 * dot(grad3[gi3], x3, y3, z3); }
        
        return 32.0f * (n0+n1+n2+n3);
    }
};

// https://github.com/Scrawk/GPU-GEMS-Improved-Perlin-Noise
// https://github.com/jbikker/tmpl8/blob/master/template/tmpl8math.cpp
struct PerlinNoise2D
{
    const int SIZE = 256;
    static const float GRADIENT2[8 * 2] = 
    {
         0,  1,
         1,  1,
         1,  0,
         1, -1,
         0, -1,
        -1, -1,
        -1,  0,
        -1,  1,
    };

    int    m_perm[SIZE + SIZE];
    float  m_permutationTable1D[SIZE];
    float2 m_gradient2D[SIZE];
    
    int    octaves = 4;
    float _Frequency = 1.0f, _Lacunarity = 2.0f, _Gain = 0.5f;

    PerlinNoise2D(uint seed = 0xe6546b64u)
    {   
        int i, j, k;
        for (i = 0; i < SIZE; i++)
            m_perm[i] = i;

        while (--i != 0)
        {
            k = m_perm[i];
            j = Abs(Random::PCG2Next(seed));
            m_perm[i] = m_perm[j];
            m_perm[j] = k;
        }

        for (i = 0; i < SIZE; i++)
            m_perm[SIZE + i] = m_perm[i];

        for (int x = 0; x < SIZE; x++)
        {
            m_permutationTable1D[x] = (float)m_perm[x] / (float)(SIZE - 1);
        }

        for (int i = 0; i < 8; i++)
        {
            float R = (GRADIENT2[i * 2 + 0] + 1.0f) * 0.5f;
            float G = (GRADIENT2[i * 2 + 1] + 1.0f) * 0.5f;
            m_gradient2D[i].x = R; 
            m_gradient2D[i].x = G; 
        }
    }

    float2 fade(float2 t)
    {
        return t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f);
    }
			
    float perm(float x)
    {
        return m_permutationTable1D[(int)(x * SIZE)];
    }
			
    float grad(float x, float2 p)
    {
        float2 g = m_gradient2D[(int)(x * 8.0 * SIZE)] * 2.0f;
        g.x -= 1.0f;
        g.y -= 1.0f;
        return g.x * p.x + g.y * p.y;
    }

    float inoise(float2 p)
    {
        float2 P = MakeVec2(FMod(Floor(p.x), 256.0f), FMod(Floor(p.y), 256.0f));	// FIND UNIT SQUARE THAT CONTAINS POINT
        p -= MakeVec2(Floor(p.x), Floor(p.y)); // FIND RELATIVE X,Y OF POINT IN SQUARE.
        float2 f = fade(p);                    // COMPUTE FADE CURVES FOR EACH OF X,Y.
			
        P = P / 256.0f;
        const float one = 1.0f / 256.0f;
				
        // HASH COORDINATES OF THE 4 SQUARE CORNERS
        float A = perm(P.x) + P.y;
        float B = perm(P.x + one) + P.y;
			 
        // AND ADD BLENDED RESULTS FROM 4 CORNERS OF SQUARE
        return Lerp(Lerp(grad(perm(A      ), p ),  
                         grad(perm(B      ), p + MakeVec2(-1.0f,  0.0f)), f.x),
                    Lerp(grad(perm(A + one), p + MakeVec2( 0.0f, -1.0f)),
                         grad(perm(B + one), p + MakeVec2(-1.0f, -1.0f)), f.x), f.y);                
    }
			
    // fractal sum, range -1.0 - 1.0
    float fbm(float2 p)
    {
        float freq = _Frequency, amp = 0.5;
        float sum = 0;	
        for(int i = 0; i < octaves; i++) 
        {
            sum  += inoise(p * freq) * amp;
            freq *= _Lacunarity;
            amp  *= _Gain;
        }
        return sum;
    }

    // fractal abs sum, range 0.0 - 1.0
    float turbulence(float2 p)
    {
        float sum = 0;
        float freq = _Frequency, amp = 1.0;
        for(int i = 0; i < octaves; i++) 
        {
            sum  += Abs(inoise(p*freq))*amp;
            freq *= _Lacunarity;
            amp  *= _Gain;
        }
        return sum;
    }
			
    // Ridged multifractal, range 0.0 - 1.0
    // See "Texturing & Modeling, A Procedural Approach", Chapter 12
    float ridge(float h, float offset)
    {
        h = Abs(h);
        h = offset - h;
        h = h * h;
        return h;
    }
			
    float ridgedmf(float2 p, float offset)
    {
        float sum = 0;
        float freq = _Frequency, amp = 0.5;
        float prev = 1.0;
        for(int i = 0; i < octaves; i++) 
        {
            float n = ridge(inoise(p*freq), offset);
            sum  += n * amp * prev;
            prev  = n;
            freq *= _Lacunarity;
            amp  *= _Gain;
        }
        return sum;
    }
};
