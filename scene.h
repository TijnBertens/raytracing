//
// Created by s150818 on 14-2-2019.
//

#ifndef RAYTRACING_SCENE_H
#define RAYTRACING_SCENE_H

#include <random>
#include "math_utils.h"
#include "color.h"
#include "geometries.h"
#include "ray.h"
#include "material.h"

//  ----------------------
//      Scene objects
//  ----------------------

/**
 * 3D sphere.
 */
struct SphereObject {
    Sphere sphere;              // sphere geometry
    PBMaterial material;        // physically based material
};

/**
 * 3D triangle.
 */
struct TriangleObject {
    Triangle triangle;          // triangle geometry
    PBMaterial material;        // physically based material
};

/**
 * Point light.
 */
struct PointLight {
    Vec3f position;          // world space position
    Color color;                // color of emitted light
};

/**
 * Representation of a camera object.
 */
struct Camera {
    Vec3f position;          // world space position
    Vec3f viewDirection;     // world space view direction (normalized)
    Vec3f upVector;          // world space up vector (normalized)

    f32 fovy;                   // vertical field of view
    f32 nearClippingDistance;   // distance from camera to screen
    f32 aspectRatio;            // width over height
};

/**
 * Calculates the world space position of a pixel on the screen of a camera.
 */
Vec3f pixelToWorldSpace(Camera camera, u32 x, u32 y, u32 width, u32 height) {
    Vec3f screenCenter = camera.position + (camera.viewDirection * camera.nearClippingDistance);

    Vec3f screenVerticalDirection = camera.upVector;
    Vec3f screenHorizontalDirection = camera.upVector.cross(camera.viewDirection).normalized();

    f32 screenHeight = 2 * (f32) tan((camera.fovy / 2.0) * PI_32 / 180.0) * camera.nearClippingDistance;
    f32 screenWidth =  screenHeight * camera.aspectRatio;

    Vec3f botLeft =
            screenCenter
            - ((screenWidth / 2) * screenHorizontalDirection)
            - ((screenHeight / 2) * screenVerticalDirection);

    Vec3f botRight =
            screenCenter
            + ((screenWidth / 2) * screenHorizontalDirection)
            - ((screenHeight / 2) * screenVerticalDirection);

    Vec3f topLeft =
            screenCenter
            - ((screenWidth / 2) * screenHorizontalDirection)
            + ((screenHeight / 2) * screenVerticalDirection);

    f32 u = (f32) x / (f32) width;
    f32 v = (f32) y / (f32) height;

    Vec3f result =
            botLeft + (u * (botRight - botLeft)) + (v * (topLeft - botLeft));

    return result;
}

Vec3f pixelToWorldSpaceRand(Camera camera, u32 x, u32 y, u32 width, u32 height) {
    Vec3f screenCenter = camera.position + (camera.viewDirection * camera.nearClippingDistance);

    Vec3f screenVerticalDirection = camera.upVector;
    Vec3f screenHorizontalDirection = camera.upVector.cross(camera.viewDirection).normalized();

    f32 screenHeight = 2 * (f32) tan((camera.fovy / 2.0) * PI_32 / 180.0) * camera.nearClippingDistance;
    f32 screenWidth =  screenHeight * camera.aspectRatio;

    Vec3f botLeft =
            screenCenter
            - ((screenWidth / 2) * screenHorizontalDirection)
            - ((screenHeight / 2) * screenVerticalDirection);

    Vec3f botRight =
            screenCenter
            + ((screenWidth / 2) * screenHorizontalDirection)
            - ((screenHeight / 2) * screenVerticalDirection);

    Vec3f topLeft =
            screenCenter
            - ((screenWidth / 2) * screenHorizontalDirection)
            + ((screenHeight / 2) * screenVerticalDirection);

    f32 u = (f32) x / (f32) width;
    f32 v = (f32) y / (f32) height;

    Vec3f result =
            botLeft + (u * (botRight - botLeft)) + (v * (topLeft - botLeft));

    f32 r1 = (float) rand() / RAND_MAX;
    f32 r2 = (float) rand() / RAND_MAX;

    f32 pixWidth = screenWidth / width;
    f32 pixHeight = screenHeight / height;

    result = result + r1 * pixWidth * screenHorizontalDirection;
    result = result + r2 * pixHeight * screenVerticalDirection;

    return result;
}

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

//  ----------------------
//    Scene intersection
//  ----------------------

/**
 * Container that stores objects in a scene. Thus representing an environment filled with 3D objects.
 */
struct Scene {
    Color backgroundColor;

    SphereObject *spheres;
    u32 numSpheres;

    TriangleObject *triangles;
    u32 numTriangles;

    PointLight *lights;
    u32 numLights;
};

/**
 * Result of intersecting a ray with a scene.
 */
struct SceneIntersectReport {
    RayHit hit;
    PBMaterial hitMaterial;
};

/**
 * Traces a ray through a scene of spheres and triangles.
 * Returns the closest hit (or an empty hit with hit.hit = false if there is no hit) and the color at that point.
 */
SceneIntersectReport intersectScene(Ray ray, Scene scene) {
    SceneIntersectReport result = {};

    f32 closestTOI = FP_INFINITE;
    RayHit closestHit = {};

    for(u32 i = 0; i < scene.numSpheres; i++) {
        Sphere sphere = scene.spheres[i].sphere;
        RayHit hit = intersect(ray, sphere);

        if(hit.hit && hit.TOI < closestTOI) {
            closestHit = hit;
            closestTOI = hit.TOI;
            result.hitMaterial = scene.spheres[i].material;
        }
    }

    for(u32 i = 0; i < scene.numTriangles; i++) {
        Triangle triangle = scene.triangles[i].triangle;
        RayHit hit = intersect(ray, triangle);

        if(hit.hit && hit.TOI < closestTOI) {
            closestHit = hit;
            closestTOI = hit.TOI;
            result.hitMaterial = scene.triangles[i].material;
        }
    }

    result.hit = closestHit;

    return result;
}

/**
 * Trace a ray through a scene. Returns the color at the intersection point.
 * @param traceDepth How many levels of recursive reflections will be sampled.
 */
Color traceThroughScene(Ray ray, Scene scene, u32 traceDepth = 5) {
    if(traceDepth == 0) return scene.backgroundColor;

    SceneIntersectReport sceneIntersect = intersectScene(ray, scene);

    // Did not hit anything, so we can return the background color
    if(!sceneIntersect.hit.hit) {
        return scene.backgroundColor;
    }

    Color surfaceColor = {0, 0, 0, 1};
    Color ambientTerm = 0.03 *sceneIntersect.hitMaterial.albedo * sceneIntersect.hitMaterial.ao;
    Vec3f viewDirection = (ray.start - sceneIntersect.hit.hitPosition).normalized();

    // Apply lighting
    for(u32 i = 0; i < scene.numLights; i++) {
        Vec3f hitToLight = (scene.lights[i].position - sceneIntersect.hit.hitPosition);

        Ray shadowRay = {};
        shadowRay.start = sceneIntersect.hit.hitPosition;
        shadowRay.direction = hitToLight.normalized();

        SceneIntersectReport shadowTrace = intersectScene(shadowRay, scene);

        if(!(shadowTrace.hit.hit && shadowTrace.hit.TOI <= hitToLight.length())) {
            Vec3f lightDirection = (scene.lights[i].position - sceneIntersect.hit.hitPosition).normalized();
            surfaceColor = surfaceColor + calculateRadiancePBR(sceneIntersect.hit.hitNormal, sceneIntersect.hitMaterial, viewDirection, lightDirection, scene.lights[i].color);
        }
    }

    //reflection

    if(sceneIntersect.hitMaterial.roughness < 0.8f) {
        Ray reflectionRay = {};
        reflectionRay.start = sceneIntersect.hit.hitPosition;
        reflectionRay.direction = ray.direction.reflectIn(sceneIntersect.hit.hitNormal);
        Color reflectionColor = traceThroughScene(reflectionRay, scene, traceDepth - 1);

        surfaceColor = surfaceColor + calculateRadiancePBR(sceneIntersect.hit.hitNormal, sceneIntersect.hitMaterial, viewDirection, reflectionRay.direction, reflectionColor);
    }

    surfaceColor = surfaceColor + ambientTerm;

    //todo: refraction
    surfaceColor.a = 1; //todo: handle alpha
    return surfaceColor;
}


//  ------------------------------------------------------------------------------
//  ------------------------------------------------------------------------------

#endif //RAYTRACING_SCENE_H
