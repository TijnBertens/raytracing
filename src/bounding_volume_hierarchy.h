//
// Created by Tijn Bertens on 11-2-2021.
//

#ifndef RAYTRACING_BOUNDING_VOLUME_HIERARCHY_H
#define RAYTRACING_BOUNDING_VOLUME_HIERARCHY_H

#include "types.h"
#include "math_utils.h"
#include "geometries.h"
#include "scene.h"

/**
 * Axis aligned bounding box with an associated centroid.
 */
struct BVH_AABB : AABB {
    Vec3f centroid;

    BVH_AABB() = default;
    BVH_AABB(Vec3f min, Vec3f max, Vec3f centroid) : AABB() {
        this->min = min;
        this->max = max;
        this->centroid = centroid;
    }
};

/**
 * Keeps track of a triangle and its associated bounding box.
 */
struct BVH_Triangle {
    BVH_AABB boundingBox;
    TriangleObject *triangle;
};

/**
 * A node in a bounding volume hierarchy.
 */
struct BVH_Node {
    BVH_AABB boundingBox;

    u32 child1;
    u32 child2;

    u32 numTriangles;
    u32 triangleOffset;
};

/**
 * Represents a bounding volume hierarchy.
 */
struct BVH {
    u32 numTriangles;
    BVH_Triangle *triangles;        // list of triangles (with associated bounding boxes)
    u32 *triangleIDs;

    u32 numNodes;
    BVH_Node *nodes;                // list of nodes
};

/**
 * Default empty box. Is used as starting point for a growing bounding box.
 */
const BVH_AABB EMPTY_BOX({INFINITY, INFINITY, INFINITY},
                         {-INFINITY, -INFINITY, -INFINITY},
                         {0,0,0});

/**
 * Calculate the centroid of an AABB.
 */
Vec3f calculateCentroid(const BVH_AABB *a) {
    return a->min + 0.5f * (a->max - a->min);
}

/**
 * Creates a new AABB, that is the given AABB grown to fit a given point.
 * Does not update the centroid, used for performance reasons.
 */
BVH_AABB combineNoCentroid(const BVH_AABB *a, const Vec3f p) {
    BVH_AABB result = {};

    result.min = Vec3f::min(a->min, p);
    result.max = Vec3f::max(a->max, p);

    return result;
}

/**
 *  Creates a new AABB, that is the given AABB grown to fit the other given AABB.
 *  Does not update the centroid, used for performance reasons.
 */
BVH_AABB combineNoCentroid(const BVH_AABB *a, const BVH_AABB *b) {
    BVH_AABB result = {};

    result.min = Vec3f::min(a->min, b->min);
    result.max = Vec3f::max(a->max, b->max);

    return result;
}

/**
 * Calculates the surface area of a given AABB.
 */
f32 SA(const BVH_AABB *box) {
    Vec3f size = {box->max.x - box->min.x, box->max.y - box->min.y, box->max.z - box->min.z};
    return (2 * size.x * size.y) + (2 * size.y * size.x) + (2 * size.z * size.y);
}

/**
 * Recursively constructs a BVH.
 * I made this implementation based on the technique proposed by Ingo Wald in 2007:
 * https://doi.org/10.1109/RT.2007.4342588
 */
static u32 recursivelyConstruct(BVH *bvh, u32 tBegin, u32 tEnd, BVH_AABB vb, BVH_AABB cb) {

      // uncomment for extra logging
//    printf("I was created with: \n");
//    printf("tBegin: %u, tEnd %u \n", tBegin, tEnd);
//    printf("vb: {%f, %f, %f}  {%f, %f, %f}\n", vb.min.x, vb.min.y, vb.min.z, vb.max.x, vb.max.y, vb.max.z);
//    printf("cb: {%f, %f, %f}  {%f, %f, %f}\n", cb.min.x, cb.min.y, cb.min.z, cb.max.x, cb.max.y, cb.max.z);
//    printf("\n");

    // base case, 4 or fewer triangles remain in the bounding volume
    if(tEnd - tBegin <= 4) {
        u32 nodeIndex = bvh->numNodes++;
        BVH_Node *node = &bvh->nodes[nodeIndex];
        node->boundingBox = vb;
        node->numTriangles = tEnd -tBegin;
        node->triangleOffset = tBegin;

        return nodeIndex;
    }

    // find dominant axis to split on

    u32 dominantAxis = 0; // 0 = x, 1 = y, 2 = z
    Vec3f cbSize = {cb.max.x - cb.min.x, cb.max.y - cb.min.y, cb.max.z - cb.min.z};

    if(cbSize.x >= cbSize.y && cbSize.x >= cbSize.z) {
        dominantAxis = 0;
    } else if(cbSize.y >= cbSize.x && cbSize.y >= cbSize.z) {
        dominantAxis = 1;
    } else {
        dominantAxis = 2;
    }

    // pre-compute binning constants

    const u32 K = 8;
    const f32 eps = 0.0001f;
    f32 k1 = (K*(1.0f - eps)) / (cb.max.values[dominantAxis] - cb.min.values[dominantAxis]);
    f32 k0 = cb.min.values[dominantAxis];

    // find the number of triangles and the bounding boxes for each bin

    u32 numTrianglesInBin[K] = {};
    BVH_AABB binBoundingBox[K] = {};

    for(auto & i : binBoundingBox) {
        i = EMPTY_BOX;
    }

    for(u32 i = tBegin; i < tEnd; i++) {
        BVH_AABB *bb = &bvh->triangles[bvh->triangleIDs[i]].boundingBox;
        u32 bin = k1*(bb->centroid.values[dominantAxis] - k0);

        numTrianglesInBin[bin]++;
        binBoundingBox[bin] = combineNoCentroid(&binBoundingBox[bin], bb);
    }

    // calculate costs

    f32 costL[K-1];
    f32 costR[K-1];

    // cost left of split

    u32 triangleCount = 0;
    BVH_AABB growingBB = EMPTY_BOX;

    for(u32 i = 0; i < K - 1; i++) {
        triangleCount += numTrianglesInBin[i];
        growingBB = combineNoCentroid(&growingBB, &binBoundingBox[i]);
        costL[i] = triangleCount * SA(&growingBB);
    }

    // cost right of split

    triangleCount = 0;
    growingBB = EMPTY_BOX;

    for(s32 i = K - 2; i >= 0; i--) {
        triangleCount += numTrianglesInBin[i + 1];
        growingBB = combineNoCentroid(&growingBB, &binBoundingBox[i + 1]);
        costR[i] = triangleCount * SA(&growingBB);
    }

    // find optimal split

    s32 splitNumber = -1;
    f32 optimalSplitCost = INFINITY;
    for(u32 i = 0; i < K - 1; i++) {
        f32 cost = costL[i] + costR[i];
        if(cost <= optimalSplitCost) {
            splitNumber = i;
            optimalSplitCost = cost;
        }
    }

    // repartition list

    //TODO: could be done without the copy using two iterators, may be worth it?

    u32 numIDs = tEnd - tBegin;
    u32 *tempList = (u32 *) malloc(sizeof(u32) * numIDs);
    memcpy(tempList, bvh->triangleIDs + tBegin, sizeof(u32) * numIDs);

    BVH_AABB leftVb = EMPTY_BOX;
    BVH_AABB leftCb = EMPTY_BOX;

    BVH_AABB rightVb = EMPTY_BOX;
    BVH_AABB rightCb = EMPTY_BOX;

    u32 leftCounter = tBegin;
    u32 rightCounter = tEnd - 1;
    for(u32 i = 0; i < numIDs; i++) {
        BVH_AABB *bb = &bvh->triangles[tempList[i]].boundingBox;
        u32 bin = k1*(bb->centroid.values[dominantAxis] - k0);

        if(bin <= splitNumber) {
            bvh->triangleIDs[leftCounter++] = tempList[i];
            leftVb = combineNoCentroid(&leftVb, bb);
            leftCb = combineNoCentroid(&leftCb, bb->centroid);
        } else {
            bvh->triangleIDs[rightCounter--] = tempList[i];
            rightVb = combineNoCentroid(&rightVb, bb);
            rightCb = combineNoCentroid(&rightCb, bb->centroid);
        }
    }

    // post-calculate centroids

    leftVb.centroid = calculateCentroid(&leftVb);
    leftCb.centroid = calculateCentroid(&leftCb);
    rightVb.centroid = calculateCentroid(&rightVb);
    rightCb.centroid = calculateCentroid(&rightCb);

    u32 mid = leftCounter;

    free(tempList);

    // recurse

    u32 nodeIndex = bvh->numNodes++;
    BVH_Node *node = &bvh->nodes[nodeIndex];
    node->boundingBox = vb;
    node->child1 = recursivelyConstruct(bvh, tBegin, mid, leftVb, leftCb);
    node->child2 = recursivelyConstruct(bvh, mid, tEnd, rightVb, rightCb);
    node->numTriangles = 0;
    node->triangleOffset = 0;

    return nodeIndex;
}
/**
 * Construct a BVH on a given list of triangles.
 */
BVH constructBVH(TriangleObject *triangles, u32 numTriangles) {
    BVH result = {};

    // Allocate space for triangles, ids and nodes

    result.numTriangles = numTriangles;
    result.triangles = (BVH_Triangle *) malloc(sizeof(BVH_Triangle) * numTriangles);
    result.triangleIDs = (u32 *) malloc(sizeof(u32) * numTriangles);

    // Can have at most 2N-1 nodes, so we can pre-allocate node space
    result.nodes = (BVH_Node *) malloc(sizeof(BVH_Node) * (numTriangles * 2 - 1));

    for(u32 i = 0; i < numTriangles; i++) {
        result.triangleIDs[i] = i;
    }

    // Find bounding boxes of triangle AABBs (vb) and centroids (cb)

    BVH_AABB vb = EMPTY_BOX;
    BVH_AABB cb = EMPTY_BOX;

    for(u32 i = 0; i < numTriangles; i++) {
        BVH_Triangle *t = &result.triangles[i];

        t->triangle = &triangles[i];
        t->boundingBox.min = Vec3f::min(Vec3f::min(triangles[i].triangle.A, triangles[i].triangle.B), triangles[i].triangle.C);
        t->boundingBox.max = Vec3f::max(Vec3f::max(triangles[i].triangle.A, triangles[i].triangle.B), triangles[i].triangle.C);
        t->boundingBox.centroid = calculateCentroid(&t->boundingBox);

        vb = combineNoCentroid(&vb, &t->boundingBox);
        cb = combineNoCentroid(&cb, t->boundingBox.centroid);
    }

    // Actually construct BVH
    recursivelyConstruct(&result, 0, numTriangles, vb, cb);

    return result;
}

/**
 * Free all memory associated with a given BVH.
 */
void BVH_free(BVH *bvh) {
    free(bvh->triangles);
    free(bvh->triangleIDs);
    free(bvh->nodes);
}

#endif //RAYTRACING_BOUNDING_VOLUME_HIERARCHY_H
