//
// Created by s150818 on 14-2-2019.
//

#ifndef RAYTRACING_GEOMETRY_H
#define RAYTRACING_GEOMETRY_H

#include "math_utils.h"
#include "file_io.h"

/**
 * 3D triangle, defined by three points
 */
struct Triangle {
    Vec3f A;
    Vec3f B;
    Vec3f C;
};

/**
 * 3D sphere, defined by center and radius
 */
struct Sphere {
    Vec3f position;      // center
    f32 radius;
};

/**
 * 3D mesh data. List of vertex coordinates and a list of indices that form triangles.
 */
struct Mesh {
    u32 numTriangles;

    u32 numVertices;
    Vec3f *vertices;

    u32 numIndices;
    u32 *indices;
};

/**
 * Load a mesh from a given file.
 */
bool loadMesh(Mesh *mesh, const char *file) {
    ObjContent objContent;
    if(!readInObj(&objContent, file)) {
        printf("Could not load mesh from file %s\n", file);
        return false;
    }

    // Copy vertex data

    mesh->numVertices = objContent.attrib.vertices.size() / 3;
    mesh->vertices = (Vec3f *) malloc(sizeof(Vec3f) * mesh->numVertices);
    // NOTE: kind of banking on the fact that real_t is a 32 bit float in the tiny obj library
    memcpy(mesh->vertices, objContent.attrib.vertices.data(), objContent.attrib.vertices.size() * sizeof(f32));

    // Count number of indices over all shapes

    mesh->numIndices = 0;
    for(const auto &shape : objContent.shapes) {
        mesh->numIndices += shape.mesh.indices.size();
    }
    mesh->numTriangles = mesh->numIndices / 3;

    // Copy index data

    mesh->indices = (u32 *) malloc(sizeof(u32) * mesh->numIndices);
    u32 i = 0;
    for(const auto &shape : objContent.shapes) {
        for(const auto &index : shape.mesh.indices) {
            mesh->indices[i++] = index.vertex_index;
        }
    }

    return true;
}

struct BVH_AABB {
    Vec3f min;
    Vec3f max;
    Vec3f centroid;
};

struct BVH_Triangle {
    BVH_AABB boundingBox;
    Triangle *triangle;
};

struct BVH_Node {
    BVH_AABB boundingBox;

    u32 child1;
    u32 child2;

    u32 numTriangles;
    u32 triangleOffset;
};

struct BVH {
    u32 numTriangles;
    BVH_Triangle *triangles;
    u32 *triangleIDs;

    u32 numNodes;
    BVH_Node *nodes;
};

const BVH_AABB EMPTY_BOX = {.min = {INFINITY, INFINITY, INFINITY},
                            .max = {-INFINITY, -INFINITY, -INFINITY},
                            .centroid = {0,0,0}};

Vec3f calculateCentroid(const BVH_AABB *a) {
    return a->min + 0.5f * (a->max - a->min);
}

BVH_AABB combineNoCentroid(const BVH_AABB *a, const Vec3f p) {
    BVH_AABB result = {};

    result.min = Vec3f::min(a->min, p);
    result.max = Vec3f::max(a->max, p);

    return result;
}

BVH_AABB combineNoCentroid(const BVH_AABB *a, const BVH_AABB *b) {
    BVH_AABB result = {};

    result.min = Vec3f::min(a->min, b->min);
    result.max = Vec3f::max(a->max, b->max);

    return result;
}

f32 SA(const BVH_AABB *box) {
    Vec3f size = {box->max.x - box->min.x, box->max.y - box->min.y, box->max.z - box->min.z};
    return (2 * size.x * size.y) + (2 * size.y * size.x) + (2 * size.z * size.y);
}

static u32 recursivelyConstruct(BVH *bvh, u32 tBegin, u32 tEnd, BVH_AABB vb, BVH_AABB cb) {

//    printf("I was created with: \n");
//    printf("tBegin: %u, tEnd %u \n", tBegin, tEnd);
//    printf("vb: {%f, %f, %f}  {%f, %f, %f}\n", vb.min.x, vb.min.y, vb.min.z, vb.max.x, vb.max.y, vb.max.z);
//    printf("cb: {%f, %f, %f}  {%f, %f, %f}\n", cb.min.x, cb.min.y, cb.min.z, cb.max.x, cb.max.y, cb.max.z);
//    printf("\n");

    // base case
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
    if(cbSize.x > cbSize.y && cbSize.x > cbSize.z) {
        dominantAxis = 0;
    } else if(cbSize.y > cbSize.x && cbSize.y > cbSize.z) {
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

BVH constructBVH(Triangle *triangles, u32 numTriangles) {
    BVH result = {};

    result.numTriangles = numTriangles;
    result.triangles = (BVH_Triangle *) malloc(sizeof(BVH_Triangle) * numTriangles);
    result.triangleIDs = (u32 *) malloc(sizeof(u32) * numTriangles);

    // Can have at most 2N-1 nodes, so e can pre-allocate node space
    result.nodes = (BVH_Node *) malloc(sizeof(BVH_Node) * (numTriangles * 2 - 1));

    for(u32 i = 0; i < numTriangles; i++) {
        result.triangleIDs[i] = i;
    }

    BVH_AABB vb = EMPTY_BOX;
    BVH_AABB cb = EMPTY_BOX;

    for(u32 i = 0; i < numTriangles; i++) {
        BVH_Triangle *t = &result.triangles[i];

        t->triangle = &triangles[i];
        t->boundingBox.min = Vec3f::min(Vec3f::min(triangles[i].A, triangles[i].B), triangles[i].C);
        t->boundingBox.max = Vec3f::max(Vec3f::max(triangles[i].A, triangles[i].B), triangles[i].C);
        t->boundingBox.centroid = calculateCentroid(&t->boundingBox);

        vb = combineNoCentroid(&vb, &t->boundingBox);
        cb = combineNoCentroid(&cb, t->boundingBox.centroid);
    }

    recursivelyConstruct(&result, 0, numTriangles, vb, cb);

    return result;
}

void BVH_free(BVH *bvh) {
    free(bvh->triangles);
    free(bvh->triangleIDs);
    free(bvh->nodes);
}


#endif //RAYTRACING_GEOMETRY_H

