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
 * 3D model
 */
struct Model {
    Mesh *mesh;                 // geometry and materials of the model
    Matrix4x4 transform;        // model to world space transform
};

/**
 * Point light.
 */
struct PointLight {
    Vec3f position;          // world space position
    Color color;                // color of emitted light
};

/**
 * Container that stores objects in a scene. Thus representing an environment filled with 3D objects.
 */
struct Scene {
    Color backgroundColor;

    SphereObject *spheres;
    u32 numSpheres;

    TriangleObject *triangles;
    u32 numTriangles;

    Model *models;
    u32 numModels;

    PointLight *lights;
    u32 numLights;
};


//  ------------------------------------------------------------------------------
//  ------------------------------------------------------------------------------

#endif //RAYTRACING_SCENE_H
