//
// Created by s150818 on 14-2-2019.
//

#ifndef RAYTRACING_GEOMETRY_H
#define RAYTRACING_GEOMETRY_H

#include "math_utils.h"
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

#endif //RAYTRACING_GEOMETRY_H
