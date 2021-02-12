//
// Created by Tijn Bertens on 3-10-2020.
//

#ifndef RAYTRACING_RAY_TRACER_H
#define RAYTRACING_RAY_TRACER_H

#include "math_utils.h"
#include "ray.h"
#include "scene.h"
#include "bounding_volume_hierarchy.h"

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
 * Evaluates the radiance at a surface point for a single light.
 * All input directions must be normalized!
 */
Color calculateRadiancePBR(Vec3f surfacePosition, Vec3f surfaceNormal, PBMaterial surface, Vec3f viewDirection, PointLight light) {
    Vec3f V = viewDirection;
    Vec3f N = surfaceNormal;

    Vec3f lightVector = light.position - surfacePosition;
    Vec3f lightDirection = lightVector.normalized();

    Color F0 = {0.3, 0.3, 0.3, 1};
    F0 = mix(F0, surface.albedo, surface.metallic);

    // calculate per-light radiance
    Vec3f L = lightDirection;
    Vec3f H = (V + L).normalized();
    f32 distance    = lightVector.length() / 5.0;
    f32 attenuation = 1.0 / (distance * distance); //todo: different attenuation?

    Color radiance = light.color * attenuation;

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

/**
 * Hack to incorporate reflections. Computes added radiance from a reflection ray with given direction and color.
 */
Color calculateReflectionInfluence(Vec3f surfaceNormal, PBMaterial surface, Vec3f viewDirection, Vec3f reflectionDirection, Color reflectionColor) {
    Vec3f V = viewDirection;

    Color F0 = {0.3, 0.3, 0.3, 1};
    F0 = mix(F0, surface.albedo, surface.metallic);

    // calculate per-light radiance
    Vec3f L = reflectionDirection;
    Vec3f H = (V + L).normalized();

    Color radiance= reflectionColor;
    Color F = fresnelSchlickRoughness(fmax(H.dot(V), 0.0), F0, surface.roughness);
    f32 NdotL = fmax(surfaceNormal.dot(L), 0.0); //TODO: does this improve the situation?

    return F * radiance * NdotL;
}

//  -------------
//    Ray tracer
//  -------------

/**
 * Contains all data used for raytracing a scene, including the scene, a full list of triangles, and a BVH.
 */
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
    PBMaterial hitMaterial;         // TODO: if we every switch to more advanced materials, this should be a pointer
};

/**
 * Result of tracing a ray through a scene.
 */
struct TraceReport {
    SceneIntersectReport firstHit;
    Color color;
};

/**
 * Build a full list of all triangles in a given scene.
 * This implies instantiating all models in the scene at their respective positions.
 */
TriangleObject *buildTriangleList(Scene *scene, u32 *numTriangles) {
    // Count total number of triangles
    u32 totalNumTriangles = scene->numTriangles;

    for(u32 i = 0; i < scene->numModels; i++) {
        totalNumTriangles += scene->models[i].mesh->numTriangles;
    }

    // Allocate memory for all triangles
    TriangleObject *allTriangles = (TriangleObject *) malloc(sizeof(TriangleObject) * totalNumTriangles);

    // Create triangles for all models in the scene
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

            // Create homogeneous coordinates to find world-space positions

            Vec4f v1h = Vec3f::toHomogeneous(mesh->vertices[mesh->indices[j*3+0]], false);
            Vec4f v2h = Vec3f::toHomogeneous(mesh->vertices[mesh->indices[j*3+1]], false);
            Vec4f v3h = Vec3f::toHomogeneous(mesh->vertices[mesh->indices[j*3+2]], false);

            Vec3f v1t = Vec3f::fromHomogeneous(scene->models[i].transform * v1h);
            Vec3f v2t = Vec3f::fromHomogeneous(scene->models[i].transform * v2h);
            Vec3f v3t = Vec3f::fromHomogeneous(scene->models[i].transform * v3h);

            allTriangles[tIndex].triangle.A = v1t;
            allTriangles[tIndex].triangle.B = v2t;
            allTriangles[tIndex].triangle.C = v3t;
            allTriangles[tIndex].material = &mesh->materials[mesh->materialIndices[j]];

            tIndex++;
        }
    }

    // Copy non-model triangles from the scene
    for(u32 i = 0; i < scene->numTriangles; i++) {
        allTriangles[tIndex++] = scene->triangles[i];
    }

    // Trim non-renderable triangles from the list
    u32 usedNumberOfTriangles = tIndex;
    allTriangles = (TriangleObject *) realloc(allTriangles, sizeof(TriangleObject) * usedNumberOfTriangles);

    *numTriangles = tIndex;

    return allTriangles;
}

/**
 * Create a RayTracer for a given scene.
 */
RayTracer createRayTracer(Scene *scene) {
    RayTracer result = {};

    result.scene = scene;
    result.triangles = buildTriangleList(scene, &result.numTriangles);
    printf("LOG: built triangle list for scene! %u triangles in total.\n", result.numTriangles);
    result.bvh = constructBVH(result.triangles, result.numTriangles);
    printf("LOG: BVH constructed.\n", result.numTriangles);

    return result;
}

/**
 * Intersect a ray with a BVH node. Recursively intersects with all children of that node.
 * @param t the current closest time of intersection found. Used to prune sub-trees that cant produce closer hits.
 */
SceneIntersectReport BVHNodeIntersect(const BVH *bvh, const BVH_Node *node, Ray ray, f32 t = INFINITY) {
    if(node->numTriangles != 0) {
        // If this is a leaf node

        f32 closestTOI = INFINITY;
        SceneIntersectReport closestHit = {};
        closestHit.hit.hit = false;

        // Loop over triangles in this leaf node, find closest intersection
        for(u32 i = 0; i < node->numTriangles; i++) {
            u32 idx = i + node->triangleOffset;
            BVH_Triangle *triangle = &bvh->triangles[bvh->triangleIDs[idx]];

            RayHit hit = intersect(ray, triangle->triangle->triangle);
            if(hit.hit && hit.TOI < closestTOI) {
                closestTOI = hit.TOI;
                closestHit.hit = hit;
                closestHit.hitMaterial = *triangle->triangle->material;
            }
        }

        return closestHit;
    } else {
        // If this is not a leaf node

        f32 closestTOI = t;
        SceneIntersectReport closestHit = {};
        closestHit.hit.hit = false;

        // Recursively intersect first child node
        HitCheck h1 = hitCheck(ray, bvh->nodes[node->child1].boundingBox);
        if(h1.hit && h1.TOI <= closestTOI) {
            SceneIntersectReport c1 = BVHNodeIntersect(bvh, &bvh->nodes[node->child1], ray, closestTOI);

            if(c1.hit.hit) {
                closestTOI = c1.hit.TOI;
                closestHit = c1;
            }
        }

        // Recursively intersect second child node
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

/**
 * Intersect a ray with a BVH.
 */
SceneIntersectReport intersectBVH(Ray ray, BVH bvh) {
    // Simply start a recursive intersect from the root node
    SceneIntersectReport result = BVHNodeIntersect(&bvh, &bvh.nodes[0], ray);
    return result;
}

/**
 *  Intersects a ray with the scene belonging to a given raytracer.
 */
SceneIntersectReport intersectScene(RayTracer *tracer, Ray ray) {
    // First intersect the BVH, and take that as closest hit
    SceneIntersectReport closestHit = intersectBVH(ray, tracer->bvh);

    // Go over all spheres separately, and check for closer hits
    for(u32 i = 0; i < tracer->scene->numSpheres; i++) {
        RayHit sphereHit = intersect(ray, tracer->scene->spheres[i].sphere);

        bool isFirstHit = (!closestHit.hit.hit && sphereHit.hit);
        bool isCloserHit = (sphereHit.hit && (sphereHit.TOI < closestHit.hit.TOI));

        if(isFirstHit || isCloserHit) {
            closestHit.hit = sphereHit;
            closestHit.hitMaterial = *tracer->scene->spheres[i].material;
        }
    }

    return closestHit;
}

/**
 * Use a RayTracer to trace a ray through a scene. Returns the color of the queried ray.
 * @param traceDepth How many levels of recursive reflections will be sampled.
 */
TraceReport traceRay(RayTracer *tracer, Ray ray, u32 traceDepth = 5) {
    // Base case for recursion
    if(traceDepth == 0) {
        TraceReport result = {};
        result.firstHit.hit.hit = false;
        result.color = tracer->scene->backgroundColor;
        return result;
    }

    // Intersect ray with scene
    SceneIntersectReport sceneIntersect = intersectScene(tracer, ray);

    // Did not hit anything, so we can return the background color
    if(!sceneIntersect.hit.hit) {
        TraceReport result = {};
        result.firstHit = sceneIntersect;
        result.color = tracer->scene->backgroundColor;
        return result;
    }

    Color surfaceColor = {0, 0, 0, 1};
    Color ambientTerm = 0.03 * sceneIntersect.hitMaterial.albedo * sceneIntersect.hitMaterial.ao;
    Vec3f viewDirection = (ray.start - sceneIntersect.hit.hitPosition).normalized();

    // Apply lighting
    for(u32 i = 0; i < tracer->scene->numLights; i++) {
        Vec3f hitToLight = (tracer->scene->lights[i].position - sceneIntersect.hit.hitPosition);

        Ray shadowRay = {};
        shadowRay.start = sceneIntersect.hit.hitPosition;
        shadowRay.direction = hitToLight.normalized();

        SceneIntersectReport shadowTrace = intersectScene(tracer, shadowRay);

        if(!(shadowTrace.hit.hit && shadowTrace.hit.TOI <= hitToLight.length())) {
            surfaceColor = surfaceColor + calculateRadiancePBR(
                    sceneIntersect.hit.hitPosition,
                    sceneIntersect.hit.hitNormal,
                    sceneIntersect.hitMaterial,
                    viewDirection,
                    tracer->scene->lights[i]);
        }
    }

    //reflection

    if(sceneIntersect.hitMaterial.roughness <= 0.3f) {
        Ray reflectionRay = {};
        reflectionRay.start = sceneIntersect.hit.hitPosition;
        reflectionRay.direction = ray.direction.reflectIn(sceneIntersect.hit.hitNormal);
        TraceReport reflectionReport = traceRay(tracer, reflectionRay, traceDepth - 1);

        surfaceColor = surfaceColor + calculateReflectionInfluence(
                sceneIntersect.hit.hitNormal,
                sceneIntersect.hitMaterial,
                viewDirection,
                reflectionRay.direction,
                reflectionReport.color);
    }

    surfaceColor = surfaceColor + ambientTerm;

    //todo: refraction
    surfaceColor.a = 1; //todo: handle alpha

    // Create result report
    TraceReport result = {};
    result.firstHit = sceneIntersect;
    result.color = surfaceColor;

    return result;
}



#endif //RAYTRACING_RAY_TRACER_H
