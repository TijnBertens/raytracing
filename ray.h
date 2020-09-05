//
// Created by s150818 on 14-2-2019.
//

#ifndef RAYTRACING_RAY_H
#define RAYTRACING_RAY_H

#include "color.h"
#include "math_utils.h"
#include "geometries.h"

#define HIT_EPSILON 0.0001f     // correction factor to prevent hit positions getting too close to objects

/**
 * Representation of a ray as a half line. Starts at start, goes in direction.
 */
struct Ray {
    Vec3f start;         // start position
    Vec3f direction;     // direction of the ray
};

/**
 * Contains data representing a ray hit, aka, details about the intersection of a ray and some shape.
 */
struct RayHit {
    bool hit;                   // was there a hit?

    f32 TOI;                    // time of intersection
    Vec3f hitPosition;          // world space position of intersection
    Vec3f hitNormal;            // world space normal at the hit position//
};


/**
 * Intersect a ray with a sphere.
 */
RayHit intersect(Ray ray, Sphere sphere) {
    RayHit result = {};

    // sphere formula:
    // || point - center ||^2 - r^2 = 0

    // b = direction dot (start - center)
    f32 b = ray.direction.dot(ray.start - sphere.position);

    // || start - center ||^2 - r^2
    f32 c = (ray.start - sphere.position).lengthSq() - (sphere.radius * sphere.radius);

    // no roots, no solution
    if(b*b - c < 0) {
        result.hit = false;
        return result;
    }

    // intersection times
    float t1 = -b - sqrt(b*b - c);
    float t2 = -b + sqrt(b*b - c);

    // no hit if the intersection is behind the start
    if(t1 <= 0 && t2 <= 0) {
        result.hit = false;
        return result;
    }

    // Get closest hit
    float t = fmin(t1, t2) - HIT_EPSILON;

    result.hit = true;
    result.TOI = t;

    Vec3f hitPosition = ray.start + (t) * ray.direction;

    result.hitPosition = hitPosition;
    result.hitNormal = (hitPosition - sphere.position).normalized();

    return result;
}

/**
 * Intersect a ray with a triangle.
 */
RayHit intersect(Ray ray, Triangle triangle) {
    RayHit result = {};

    // plain equation, normal N, point on triangle X:
    // N dot (X - P) = 0

    // normal of the plain of the triangle
    Vec3f planeNormal = (triangle.B - triangle.A).cross(triangle.C - triangle.B).normalized();

    // time of intersection
    // t = N dot (P - start) / (N dot direction)
    f32 t = (planeNormal.dot(triangle.A - ray.start)) / (planeNormal.dot(ray.direction));

    // no hit if the intersection is behind the start
    if(t <= 0) {
        result.hit = false;
        return result;
    }

    // apply correction
    t = t - HIT_EPSILON;

    Vec3f hitPosition = ray.start + t * ray.direction;

    Vec3f normalAB = (triangle.B - triangle.A).cross(planeNormal);
    Vec3f normalBC = (triangle.C - triangle.B).cross(planeNormal);
    Vec3f normalCA = (triangle.A - triangle.C).cross(planeNormal);

    Vec3f AH = hitPosition - triangle.A;
    Vec3f BH = hitPosition - triangle.B;
    Vec3f CH = hitPosition - triangle.C;

    // check if the hitposition is inside the bounds of the triangle
    if(normalAB.dot(AH) <=0
       && normalBC.dot(BH) <= 0
       && normalCA.dot(CH) <= 0) {

        // hit

        result.hit = true;
        result.hitPosition = hitPosition;
        result.TOI = t;
        result.hitNormal = planeNormal;
    } else {

        // no hit

        result.hit = false;
    }

    return result;
}



#endif //RAYTRACING_RAY_H
