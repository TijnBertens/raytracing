//
// Created by s150818 on 14-2-2019.
//

#ifndef RAYTRACING_SCENE_H
#define RAYTRACING_SCENE_H

#include "math_utils.h"

struct Material {
    //Color ambient;
    //Color diffuse;
    //Color specular;
    //f32 shininess;

    Color color;
    f32 shininess;

    f32 specularIntensity = 0.1f;
};

struct SphereObject {
    Sphere sphere;
    Material material;
};

struct TriangleObject {
    Triangle triangle;
    Material material;
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
