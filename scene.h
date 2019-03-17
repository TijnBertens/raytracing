//
// Created by s150818 on 14-2-2019.
//

#ifndef RAYTRACING_SCENE_H
#define RAYTRACING_SCENE_H

#include "math_utils.h"
#include "color.h"
#include "geometries.h"

struct Material {
    Color color;
    f32 shininess;

    f32 specularIntensity = 0.1f;
};

//  ----------------------
// Physically based material
//  ----------------------

struct PBMaterial {
    Color albedo;
    f32 metallic;
    f32 roughness;
    f32 ao;
};

static const PBMaterial PBM_ROUGH_RED = {
        {1, 0, 0, 1},
        0,
        0.8,
        1
};

static const PBMaterial PBM_METALLIC_GREEN = {
        {0, 1, 0, 1},
        1,
        0.2,
        1
};

static const PBMaterial PBM_SMOOTH_BLUE = {
        {0, 0, 1, 1},
        0,
        0.1,
        1
};

static const PBMaterial PBM_GRAY = {
        {0.5, 0.5, 0.5, 1},
        0,
        0.4,
        1
};


//  ----------------------
//  ----------------------


struct SphereObject {
    Sphere sphere;
    PBMaterial material;
};

struct TriangleObject {
    Triangle triangle;
    PBMaterial material;
};

struct PointLight {
    Vector3D position;

    f32 ambientIntensity;
};


struct Scene {
    SphereObject *spheres;
    u32 numSpheres;

    TriangleObject *triangles;
    u32 numTriangles;

    PointLight *lights;
    u32 numLights;
};

#include "scene.cpp"

#endif //RAYTRACING_SCENE_H
