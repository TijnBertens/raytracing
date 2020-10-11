//
// Created by Tijn Bertens on 3-10-2020.
//

#ifndef RAYTRACING_RAY_TRACER_H
#define RAYTRACING_RAY_TRACER_H

#include "math_utils.h"
#include "ray.h"
#include "scene.h"

//  ------------------------------------------------------------------------------
// The following functions are implementations of Physically Based Shading functions
// all credits to: https://learnopengl.com/PBR/Lighting
//  ------------------------------------------------------------------------------

/**
 * Credits to:
 * https://learnopengl.com/PBR/Lighting
 */
Color fresnelSchlick(f32 cosTheta, Color F0)
{
    Color one = {1, 1, 1, 1};
    return F0 + (one - F0) * pow(1.0 - cosTheta, 5.0);
}

Color fresnelSchlickRoughness(float cosTheta, Color F0, f32 roughness)
{
    f32 invR = 1.0f - roughness;
    Color invRc = {invR, invR, invR, invR};
    return F0 + (cmax(invRc, F0) - F0) * pow(1.0 - cosTheta, 5.0);
}

float DistributionGGX(Vec3f N, Vec3f H, f32 roughness)
{
    f32 a      = roughness*roughness;
    f32 a2     = a*a;
    f32 NdotH  = fmax(N.dot(H), 0.0);
    f32 NdotH2 = NdotH*NdotH;

    f32 num   = a2;
    f32 denom = (NdotH2 * (a2 - 1.0) + 1.0);
    denom = PI_32 * denom * denom;

    return num / denom;
}

f32 GeometrySchlickGGX(f32 NdotV, f32 roughness)
{
    f32 r = (roughness + 1.0);
    f32 k = (r*r) / 8.0;

    f32 num   = NdotV;
    f32 denom = NdotV * (1.0 - k) + k;

    return num / denom;
}
f32 GeometrySmith(Vec3f N, Vec3f V, Vec3f L, f32 roughness)
{
    f32 NdotV = fmax(N.dot(V), 0.0);
    f32 NdotL = fmax(N.dot(L), 0.0);
    f32 ggx2  = GeometrySchlickGGX(NdotV, roughness);
    f32 ggx1  = GeometrySchlickGGX(NdotL, roughness);

    return ggx1 * ggx2;
}

/**
 * Evaluates the radiance at a surface point for a single incoming light direction.
 * All input directions must be normalized!
 */
Color calculateRadiancePBR(Vec3f surfaceNormal, PBMaterial surface, Vec3f viewDirection, Vec3f lightDirection, Color lightRadiance) {
    Vec3f V = viewDirection;
    Vec3f N = surfaceNormal;

    Color F0 = {0.3, 0.3, 0.3, 1};
    F0 = mix(F0, surface.albedo, surface.metallic);

    // calculate per-light radiance
    Vec3f L = lightDirection;
    Vec3f H = (V + L).normalized();
    //f32 distance    = (light.position - surfacePosition).length();
    //f32 attenuation = 1.0 / (distance * distance); //todo: different attenuation

    Color radiance = lightRadiance;

    // cook-torrance brdf
    f32 NDF = DistributionGGX(N, H, surface.roughness);
    f32 G = GeometrySmith(N, V, L, surface.roughness);
    Color F = fresnelSchlickRoughness(fmax(H.dot(V), 0.0), F0, surface.roughness);

    Color kS = F;
    Color kD = {1 - kS.r, 1 - kS.g, 1 - kS.b, 1}; //todo: what about alpha here?  //vec3(1.0) - kS; // here
    kD = kD * (1.0 - surface.metallic);

    Color numerator = NDF * G * F;
    f32 denominator = 4.0 * fmax(N.dot(V), 0.0) * fmax(N.dot(L), 0.0);
    Color specular = numerator / fmax(denominator, 0.001);

    // return outgoing radiance
    f32 NdotL = fmax(N.dot(L), 0.0);
    return (kD * surface.albedo / PI_32 + specular) * radiance * NdotL;
}

//  ----------------------------------------------------
//    Bounding Volume Hierarchy acceleration structure
//  ----------------------------------------------------

struct BVH_AABB : AABB {
    Vec3f centroid;

    BVH_AABB() = default;
    BVH_AABB(Vec3f min, Vec3f max, Vec3f centroid) : AABB() {
        this->min = min;
        this->max = max;
        this->centroid = centroid;
    }
};

struct BVH_Triangle {
    BVH_AABB boundingBox;
    TriangleObject *triangle;
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

const BVH_AABB EMPTY_BOX({INFINITY, INFINITY, INFINITY},
        {-INFINITY, -INFINITY, -INFINITY},
        {0,0,0});

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

BVH constructBVH(TriangleObject *triangles, u32 numTriangles) {
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
        t->boundingBox.min = Vec3f::min(Vec3f::min(triangles[i].triangle.A, triangles[i].triangle.B), triangles[i].triangle.C);
        t->boundingBox.max = Vec3f::max(Vec3f::max(triangles[i].triangle.A, triangles[i].triangle.B), triangles[i].triangle.C);
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

//  -------------
//    Ray tracer
//  -------------

struct RayTracer {
    Scene *scene;

    TriangleObject *triangles;
    u32 numTriangles;

    BVH bvh;
};

/**
 * Result of intersecting a ray with a scene.
 */
struct SceneIntersectReport {
    RayHit hit;
    PBMaterial hitMaterial;
};

TriangleObject *buildTriangleList(Scene *scene, u32 *numTriangles) {
    u32 totalNumTriangles = scene->numTriangles;

    for(u32 i = 0; i < scene->numModels; i++) {
        totalNumTriangles += scene->models[i].mesh->numTriangles;
    }

    TriangleObject *allTriangles = (TriangleObject *) malloc(sizeof(TriangleObject) * totalNumTriangles);

    u32 tIndex = 0;
    for(u32 i = 0; i < scene->numModels; i++) {
        Mesh *mesh = scene->models[i].mesh;
        for(u32 j = 0; j < mesh->numTriangles; j++) {
            Vec3f v1 = mesh->vertices[mesh->indices[j*3+0]];
            Vec3f v2 = mesh->vertices[mesh->indices[j*3+1]];
            Vec3f v3 = mesh->vertices[mesh->indices[j*3+2]];

            Vec3f min = Vec3f::min(Vec3f::min(v1, v2), v3);
            Vec3f max = Vec3f::max(Vec3f::max(v1, v2), v3);

            // If two (or more) dimensions of the AABB are 0, then the triangle has no surface area
            // and won't be rendered.
            if((min.x == max.x && min.y == max.y) ||
               (min.x == max.x && min.z == max.z) ||
               (min.y == max.y && min.z == max.z)) {
                continue;
            }

            Vec4f v1h = Vec3f::toHomogeneous(mesh->vertices[mesh->indices[j*3+0]], false);
            Vec4f v2h = Vec3f::toHomogeneous(mesh->vertices[mesh->indices[j*3+1]], false);
            Vec4f v3h = Vec3f::toHomogeneous(mesh->vertices[mesh->indices[j*3+2]], false);

            Vec3f v1t = Vec3f::fromHomogeneous(scene->models[i].transform * v1h);
            Vec3f v2t = Vec3f::fromHomogeneous(scene->models[i].transform * v2h);
            Vec3f v3t = Vec3f::fromHomogeneous(scene->models[i].transform * v3h);

            allTriangles[tIndex].triangle.A = v1t;
            allTriangles[tIndex].triangle.B = v2t;
            allTriangles[tIndex].triangle.C = v3t;
            allTriangles[tIndex].material = PBM_SMOOTH_OFF_WHITE; //TODO: load materials for meshes

            tIndex++;
        }
    }

    for(u32 i = 0; i < scene->numTriangles; i++) {
        allTriangles[tIndex++] = scene->triangles[i];
    }

    u32 usedNumberOfTriangles = tIndex;
    allTriangles = (TriangleObject *) realloc(allTriangles, sizeof(TriangleObject) * usedNumberOfTriangles);

    *numTriangles = tIndex;
    return allTriangles;
}

RayTracer createRayTracer(Scene *scene) {
    RayTracer result = {};

    result.scene = scene;
    result.triangles = buildTriangleList(scene, &result.numTriangles);
    printf("LOG: built triangle list for scene! %u triangles in total.\n", result.numTriangles);
    result.bvh = constructBVH(result.triangles, result.numTriangles);
    printf("LOG: BVH constructed.\n", result.numTriangles);

    return result;
}

SceneIntersectReport BVHNodeIntersect(const BVH *bvh, const BVH_Node *node, Ray ray, f32 t = INFINITY) {
    // TODO: numTriangles should maybe be s32 so we can use -1 as an invalid value rather than 0
    // TODO: pass down closest TOI to prune subtrees that will never result in a closer hit
    if(node->numTriangles != 0) {
        f32 closestTOI = INFINITY;
        SceneIntersectReport closestHit = {};
        closestHit.hit.hit = false;

        for(u32 i = 0; i < node->numTriangles; i++) {
            u32 idx = i + node->triangleOffset;
            BVH_Triangle *triangle = &bvh->triangles[bvh->triangleIDs[idx]];

            RayHit hit = intersect(ray, triangle->triangle->triangle);
            if(hit.hit && hit.TOI < closestTOI) {
                closestTOI = hit.TOI;
                closestHit.hit = hit;
                closestHit.hitMaterial = triangle->triangle->material;
            }
        }

        return closestHit;
    } else {
        f32 closestTOI = t;
        SceneIntersectReport closestHit = {};
        closestHit.hit.hit = false;

        HitCheck h1 = hitCheck(ray, bvh->nodes[node->child1].boundingBox);
        if(h1.hit && h1.TOI <= closestTOI) {
            SceneIntersectReport c1 = BVHNodeIntersect(bvh, &bvh->nodes[node->child1], ray, closestTOI);

            if(c1.hit.hit) {
                closestTOI = c1.hit.TOI;
                closestHit = c1;
            }
        }

        HitCheck h2 = hitCheck(ray, bvh->nodes[node->child2].boundingBox);
        if(h2.hit && h2.TOI <= closestTOI) {
            SceneIntersectReport c2 = BVHNodeIntersect(bvh, &bvh->nodes[node->child2], ray, closestTOI);

            if(c2.hit.hit && c2.hit.TOI < closestTOI) {
                closestHit = c2;
            }
        }

        return closestHit;
    }
}

SceneIntersectReport intersectBVH(Ray ray, BVH bvh) {
    SceneIntersectReport result = {};

    result = BVHNodeIntersect(&bvh, &bvh.nodes[0], ray);

    return result;
}

/**
 * Use a RayTracer to trace a ray through a scene. Returns the color of the queried ray.
 * @param traceDepth How many levels of recursive reflections will be sampled.
 */
Color traceRay(RayTracer *tracer, Ray ray, u32 traceDepth = 5) {
    if(traceDepth == 0) return tracer->scene->backgroundColor;

    SceneIntersectReport sceneIntersect = intersectBVH(ray, tracer->bvh);

    // Did not hit anything, so we can return the background color
    if(!sceneIntersect.hit.hit) {
        return tracer->scene->backgroundColor;
    }

    Color surfaceColor = {0, 0, 0, 1};
    Color ambientTerm = 0.03 *sceneIntersect.hitMaterial.albedo * sceneIntersect.hitMaterial.ao;
    Vec3f viewDirection = (ray.start - sceneIntersect.hit.hitPosition).normalized();

    // Apply lighting
    for(u32 i = 0; i < tracer->scene->numLights; i++) {
        Vec3f hitToLight = (tracer->scene->lights[i].position - sceneIntersect.hit.hitPosition);

        Ray shadowRay = {};
        shadowRay.start = sceneIntersect.hit.hitPosition;
        shadowRay.direction = hitToLight.normalized();

        SceneIntersectReport shadowTrace = intersectBVH(shadowRay, tracer->bvh);

        if(!(shadowTrace.hit.hit && shadowTrace.hit.TOI <= hitToLight.length())) {
            Vec3f lightDirection = (tracer->scene->lights[i].position - sceneIntersect.hit.hitPosition).normalized();
            surfaceColor = surfaceColor + calculateRadiancePBR(sceneIntersect.hit.hitNormal, sceneIntersect.hitMaterial, viewDirection, lightDirection, tracer->scene->lights[i].color);
        }
    }

    //reflection

    if(sceneIntersect.hitMaterial.roughness < 0.8f) {
        Ray reflectionRay = {};
        reflectionRay.start = sceneIntersect.hit.hitPosition;
        reflectionRay.direction = ray.direction.reflectIn(sceneIntersect.hit.hitNormal);
        Color reflectionColor = traceRay(tracer, reflectionRay, traceDepth - 1);

        surfaceColor = surfaceColor + calculateRadiancePBR(sceneIntersect.hit.hitNormal, sceneIntersect.hitMaterial, viewDirection, reflectionRay.direction, reflectionColor);
    }

    surfaceColor = surfaceColor + ambientTerm;

    //todo: refraction
    surfaceColor.a = 1; //todo: handle alpha
    return surfaceColor;
}



#endif //RAYTRACING_RAY_TRACER_H
