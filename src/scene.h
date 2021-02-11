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
#include "camera.h"

/**
 * 3D sphere.
 */
struct SphereObject {
    Sphere sphere;              // sphere geometry
    const PBMaterial *material; // physically based material
};

/**
 * 3D triangle.
 */
struct TriangleObject {
    Triangle triangle;          // triangle geometry
    const PBMaterial *material; // physically based material
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
    Vec3f position;             // world space position
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


#endif //RAYTRACING_SCENE_H
