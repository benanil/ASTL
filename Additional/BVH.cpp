
#include "BVH.hpp"

AX_NAMESPACE

// todo make non simd version
#ifdef AX_SUPPORT_SSE

// from Jacco Bikker's tutorial I've optimized with SIMD. the unlicense
// https://jacco.ompf2.com/2022/04/13/how-to-build-a-bvh-part-1-basics/

// use simd more
struct aabb 
{ 
	vec_t bmin, bmax;
	
	aabb() : bmin(VecSet1(1e30f)), bmax(VecSet1(-1e30f)) {}

	void grow(vec_t p)
	{ 
		bmin = VecMin(bmin, p), bmax = VecMax(bmax, p);
	}

	void grow(const Tri* tri) 
	{ 
		bmin = VecMin(bmin, tri->v0);
		bmin = VecMin(bmin, tri->v1); 
		bmin = VecMin(bmin, tri->v2);
		bmax = VecMax(bmax, tri->v0);
		bmax = VecMax(bmax, tri->v1);
		bmax = VecMax(bmax, tri->v2);
	}
	
	void grow(aabb other)
	{
		if (AX_LIKELY(VecGetX(other.bmin) != 1e30f))
		{
			bmin = VecMin(bmin, other.bmin);
			bmax = VecMax(bmax, other.bmin);
			bmin = VecMin(bmin, other.bmax);
			bmax = VecMax(bmax, other.bmax);
		}
	}

	float area() 
	{ 
		vec_t e = VecSub(bmax, bmin); // box extent
		vec_t eSurface = VecMask(VecMul(e, VecSwizzle(e, 1, 2, 0, 0)), VecSelect1110);
		return VecGetX(VecHSum(eSurface));
	}
};

static uint totalNodesUsed = 0;
static BVHNode* nodes;

#define GetCenteroid(tri, axis) VecGetW(*((vec_t*)(tri) + axis)) 

static void UpdateNodeBounds(BVHNode* bvhNode, const Tri* tris, uint nodeIdx)
{
    BVHNode* node = bvhNode + nodeIdx;
    vec_t nodeMin = VecSet1(1e30f), nodeMax = VecSet1(-1e30f);
    const vec_t* leafPtr = (const vec_t*)(tris + node->leftFirst);

    AX_ASSUME(node->triCount > 0);
    for (uint i = 0; i < node->triCount; i++)
    {
    	nodeMin = VecMin(nodeMin, leafPtr[0]);
    	nodeMin = VecMin(nodeMin, leafPtr[1]);
    	nodeMin = VecMin(nodeMin, leafPtr[2]);
    	
    	nodeMax = VecMax(nodeMax, leafPtr[0]);
    	nodeMax = VecMax(nodeMax, leafPtr[1]);
    	nodeMax = VecMax(nodeMax, leafPtr[2]);
    	leafPtr += 5; // +3 for vertexPositions + 1 for (texcoords + material index) + 1 for normals
    }
    
    Vec3Store(&node->aabbMin.x, nodeMin);
    Vec3Store(&node->aabbMax.x, nodeMax);
}

purefn float CalculateCost(const BVHNode* node)
{ 
	vec_t e = VecSub(node->maxv, node->minv); // box extent
	vec_t eSurface = VecMask(VecMul(e, VecSwizzle(e, 1, 2, 0, 0)), VecSelect1110);
	return node->triCount * VecGetX(VecHSum(eSurface));
}

static float EvaluateSAH(const BVHNode* node, Tri* tris, int axis, float pos)
{
	aabb leftBox, rightBox;
	int leftCount = 0, rightCount = 0;
	
	AX_ASSUME(node->triCount > 0);
	for (uint i = 0; i < node->triCount; ++i)
	{
		Tri* triangle = tris + node->leftFirst + i;
		if (GetCenteroid(triangle, axis) < pos) {
			leftCount++;
			leftBox.grow(triangle);
		}
		else {
			rightCount++;
			rightBox.grow(triangle);
		}
	}
	float cost = leftCount * leftBox.area() + rightCount * rightBox.area();
	return cost > 0.0f ? cost : 1e30f;
}

static float FindBestSplitPlane(const BVHNode* node, Tri* tris, int* outAxis, float* splitPos)
{
	float bestCost = 1e30f;
	const uint triCount = node->triCount, leftFirst = node->leftFirst;
	AX_ASSUME(triCount > 0);

	for (int axis = 0; axis < 3; ++axis)
	{
		float boundsMin = 1e30f, boundsMax = -1e30f;
		
		for (uint i = 0; i < triCount; ++i)
		{
			Tri* triangle = tris + leftFirst + i;
			float val = GetCenteroid(triangle, axis);
			boundsMin = MIN(boundsMin, val);
			boundsMax = MAX(boundsMax, val);
		}

		if (boundsMax == boundsMin) continue;

		const int BINS = 8;
		struct Bin { aabb bounds; uint triCount = 0; };
		Bin bin[BINS] = {};
		float scale = float(BINS) / (boundsMax - boundsMin);
		for (uint i = 0; i < triCount; i++)
		{
			Tri* triangle = tris + leftFirst + i;
			float centeroid = GetCenteroid(triangle, axis);
			int binIdx = MIN(BINS - 1, (int)((centeroid - boundsMin) * scale));
			bin[binIdx].triCount++;
			bin[binIdx].bounds.grow(triangle);
		}
		
		float leftArea[BINS - 1], rightArea[BINS - 1];
		int leftCount[BINS - 1], rightCount[BINS - 1];

		int leftSum = 0, rightSum = 0;
		aabb leftBox{}, rightBox{};

		for (int i = 0; i < BINS - 1; i++)
		{
			leftSum += bin[i].triCount;
			leftCount[i] = leftSum;
			leftBox.grow( bin[i].bounds );
			leftArea[i] = leftBox.area();
			rightSum += bin[BINS - 1 - i].triCount;
			rightCount[BINS - 2 - i] = rightSum;
			rightBox.grow( bin[BINS - 1 - i].bounds );
			rightArea[BINS - 2 - i] = rightBox.area();
		}

		scale = (boundsMax - boundsMin) / float(BINS);
		for (int i = 0; i < BINS - 1; i++)
		{
			float planeCost = leftCount[i] * leftArea[i] + rightCount[i] * rightArea[i];
			if (planeCost < bestCost)
				*splitPos = boundsMin + scale * (i + 1),
				*outAxis = axis, bestCost = planeCost;
		}
	}
	return bestCost;
}

static void SubdivideBVH(BVHNode* bvhNode, Tri* tris, uint nodeIdx)
{
	// terminate recursion
	BVHNode* node = bvhNode + nodeIdx;
	uint leftFirst = node->leftFirst, triCount = node->triCount;
	// determine split axis and position
	int axis;
	float splitPos;
	float splitCost = FindBestSplitPlane(node, tris, &axis, &splitPos);
	float nosplitCost = CalculateCost(node);
	
	if (splitCost >= nosplitCost) return;

	// in-place partition
	int i = leftFirst;
	int j = i + triCount - 1;
	while (i <= j)
	{
		if (GetCenteroid(tris + i, axis) < splitPos)
			i++;
		else {
			// swap elements with simd 2x faster
			vec_t* a = (vec_t*)(tris + i);
			vec_t* b = (vec_t*)(tris + j);
			
			vec_t t = *a;
			*a++ = *b, *b++ = t, t = *a; // swap a[0], b[0] vertex0
			*a++ = *b, *b++ = t, t = *a; // swap a[1], b[1] vertex1
			*a++ = *b, *b++ = t, t = *a; // swap a[2], b[2] vertex2
			*a++ = *b, *b++ = t, t = *a; // swap a[3], b[3] uv's + material index + normal0x
			*a = *b, *b = t;             // swap a[5], b[5] normals
			j--;
		}
	}
	// abort split if one of the sides is empty
	int leftCount = i - leftFirst;
	if (leftCount == 0 || leftCount == triCount) return;
	// create child nodes
	int leftChildIdx = totalNodesUsed++;
	int rightChildIdx = totalNodesUsed++;
	bvhNode[leftChildIdx].leftFirst = leftFirst;
	bvhNode[leftChildIdx].triCount = leftCount;
	bvhNode[rightChildIdx].leftFirst = i;
	bvhNode[rightChildIdx].triCount = triCount - leftCount;
	node->leftFirst = leftChildIdx;
	node->triCount = 0;
	UpdateNodeBounds(bvhNode, tris, leftChildIdx);
	UpdateNodeBounds(bvhNode, tris, rightChildIdx);
	// recurse
	SubdivideBVH(bvhNode, tris, leftChildIdx);
	SubdivideBVH(bvhNode, tris, rightChildIdx);
}

uint BuildBVH(Tri* tris, MeshInfo* meshes, int numMeshes, BVHNode* nodes, uint* bvhIndices)
{
	// 1239.74ms SIMD
	// 556.51ms  SIMD with custom swap
	// 6511.79ms withut
	int numTriangles = 0;
	for (int i = 0; i < numMeshes; ++i) {
		numTriangles += meshes[i].numTriangles;
	}

	// calculate triangle centroids for partitioning
	for (int i = 0; i < numTriangles; i++) { // this loop will automaticly vectorized by compiler
		// set centeroids
		Tri* tri = tris + i;
		tri->centeroidx = (tri->vertex0.x + tri->vertex1.x + tri->vertex2.x) * 0.333333f;
		tri->centeroidy = (tri->vertex0.y + tri->vertex1.y + tri->vertex2.y) * 0.333333f;
		tri->centeroidz = (tri->vertex0.z + tri->vertex1.z + tri->vertex2.z) * 0.333333f;
	}
	
	uint nodesUsedStart = totalNodesUsed;
	int currTriangle = 0;
	for (int i = 0; i < numMeshes; ++i) {
		// assign all triangles to root node
		uint rootNodeIndex = totalNodesUsed++;

		bvhIndices[i] = rootNodeIndex;
		
		BVHNode& root = nodes[rootNodeIndex];
		root.leftFirst = currTriangle, root.triCount = meshes[i].numTriangles;

		UpdateNodeBounds(nodes, tris, rootNodeIndex);
		// subdivide recursively
		SubdivideBVH(nodes, tris, rootNodeIndex);	
		currTriangle += meshes[i].numTriangles;
	}
	
	return totalNodesUsed - nodesUsedStart;
}

#endif // AX_SUPPORT_SSE

AX_END_NAMESPACE