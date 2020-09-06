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
    Vec3f hitNormal;            // world space normal at the hit position
};

/**
 * A more stripped down version of a RayHit, used for simply checking whether and when a hit happened.
 */
struct HitCheck {
    bool hit;                   // was there a hit
    f32 TOI;                    // time of intersection
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
    if (b * b - c < 0) {
        result.hit = false;
        return result;
    }

    // intersection times
    float t1 = -b - sqrt(b * b - c);
    float t2 = -b + sqrt(b * b - c);

    // no hit if the intersection is behind the start
    if (t1 <= 0 && t2 <= 0) {
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
    if (t <= 0) {
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
    if (normalAB.dot(AH) <= 0
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

#define NUMDIM    3
#define RIGHT    0
#define LEFT    1
#define MIDDLE    2

/**
 * Checks when and if a ray will hit an AABB.
 * Adapted from Andrew Woo:
 * https://web.archive.org/web/20090803054252/http://tog.acm.org/resources/GraphicsGems/gems/RayBox.c
 */
HitCheck hitCheck(Ray ray, BVH_AABB aabb) {
    HitCheck result = {};

    bool inside = true;
    char quadrant[NUMDIM];
    s32 whichPlane;
    Vec3f maxT = {};
    Vec3f candidatePlane = {};

    /* Find candidate planes; this loop can be avoided if
       rays cast all from the eye(assume perpsective view) */
    for (s32 i = 0; i < NUMDIM; i++)
        if (ray.start.values[i] < aabb.min.values[i]) {
            quadrant[i] = LEFT;
            candidatePlane.values[i] = aabb.min.values[i];
            inside = false;
        } else if (ray.start.values[i] > aabb.max.values[i]) {
            quadrant[i] = RIGHT;
            candidatePlane.values[i] = aabb.max.values[i];
            inside = false;
        } else {
            quadrant[i] = MIDDLE;
        }

    /* Ray origin inside bounding box */
    if (inside) {
        result.hit = true;
        result.TOI = 0;
        return result;
    }


    /* Calculate T distances to candidate planes */
    for (s32 i = 0; i < NUMDIM; i++)
        if (quadrant[i] != MIDDLE && ray.direction.values[i] != 0.)
            maxT.values[i] = (candidatePlane.values[i] - ray.start.values[i]) / ray.direction.values[i];
        else
            maxT.values[i] = -1.f;

    /* Get largest of the maxT's for final choice of intersection */
    whichPlane = 0;
    for (s32 i = 1; i < NUMDIM; i++)
        if (maxT.values[whichPlane] < maxT.values[i])
            whichPlane = i;

    /* Check final candidate actually inside box */
    if (maxT.values[whichPlane] < 0.) {
        result.hit = false;
        result.TOI = -1.f;
        return result;
    }
    Vec3f hitPos = {};
    for (s32 i = 0; i < NUMDIM; i++) {
        if (whichPlane != i) {
            hitPos.values[i] = ray.start.values[i] + maxT.values[whichPlane] * ray.direction.values[i];
            if (hitPos.values[i] < aabb.min.values[i] || hitPos.values[i] > aabb.max.values[i]) {
                result.hit = false;
                result.TOI = -1.f;
                return result;
            }
        } else {
            hitPos.values[i] = candidatePlane.values[i];
        }
    }

    result.hit = true;
    result.TOI = maxT.values[whichPlane];

    return result;                /* ray hits box */
}


#endif //RAYTRACING_RAY_H
