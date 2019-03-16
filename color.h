//
// Created by s150818 on 13-2-2019.
//

#ifndef RAYTRACING_COLOR_H
#define RAYTRACING_COLOR_H

#include "types.h"
#include "math_utils.h"


/**
 * RGBA formatted color, components are floats.
 */
struct Color {
    union {
        Vector4D v;
        struct {
            f32 r;
            f32 g;
            f32 b;
            f32 a;
        };
    };
};

static const Color RED      =   {1, 0, 0, 1};
static const Color GREEN    =   {0, 1, 0, 1};
static const Color BLUE     =   {0, 0, 1, 1};


#pragma pack(push, 1)
/**
 * ARGB formatted color, components are 8bit unsigned ints.
 */
struct IntColor {
    union {
        u32 color;
        struct {
            u8 b;
            u8 g;
            u8 r;
            u8 a;
        };
    };
};
#pragma pack(pop)


//----------------------
//-------Color----------
//----------------------

/**
 * Convert to 32bit IntColor.
 */
inline IntColor toIntColor(Color c) {
    IntColor result;

    result.r = (u32) (c.r * 255);
    result.g = (u32) (c.g * 255);
    result.b = (u32) (c.b * 255);
    result.a = (u32) (c.a * 255);

    return result;
}


/**
 *  Mixes two colors by averaging their rgba values.
 */
inline Color mix(Color a, Color b) {
    Color result = {};

    result.r = (a.r + b.r) / 2;
    result.g = (a.g + b.g) / 2;
    result.b = (a.b + b.b) / 2;
    result.a = (a.a + b.a) / 2;

    return result;
}

/**
 * Computes the sum of colors a and b.
 */
inline Color operator+(Color a, Color b) {
    Color result = {};

    result.r = fmin(1.0, a.r + b.r);
    result.g = fmin(1.0, a.g + b.g);
    result.b = fmin(1.0, a.b + b.b);
    result.a = fmin(1.0, a.a + b.a);

    return result;
}

/**
 * Scalar multiplication for colors.
 */
inline Color operator*(f32 s, Color c) {
    Color result = {};

    result.r = fmin(1.0, s * c.r);
    result.g = fmin(1.0, s * c.g);
    result.b = fmin(1.0, s * c.b);
    result.a = fmin(1.0, s * c.a);

    return result;
}

/**
 * Scalar multiplication for colors.
 */
inline Color operator*(Color c, f32 s) {
    Color result = {};

    result.r = fmin(1.0, s * c.r);
    result.g = fmin(1.0, s * c.g);
    result.b = fmin(1.0, s * c.b);
    result.a = fmin(1.0, s * c.a);

    return result;
}

//----------------------
//------IntColor--------
//----------------------

/**
 *  Mixes two colors by averaging their rgba values.
 */
inline IntColor mix(IntColor a, IntColor b) {
    IntColor result = {};

    result.r = ((u32) a.r + b.r) / 2;
    result.g = ((u32) a.g + b.g) / 2;
    result.b = ((u32) a.b + b.b) / 2;
    result.a = ((u32) a.a + b.a) / 2;

    return result;
}

/**
 * Computes the sum of colors a and b.
 */
inline IntColor operator+(IntColor a, IntColor b) {
    IntColor result = {};

    result.r = fmin(255, (u32)a.r + b.r);
    result.g = fmin(255, (u32)a.g + b.g);
    result.b = fmin(255, (u32)a.b + b.b);
    result.a = fmin(255, (u32)a.a + b.a);

    return result;
}

/**
 * Scalar multiplication for colors.
 */
inline IntColor operator*(f32 s, IntColor c) {
    IntColor result = {};

    result.r = fmin(255, (u8) (s * (f32)c.r));
    result.g = fmin(255, (u8) (s * (f32)c.g));
    result.b = fmin(255, (u8) (s * (f32)c.b));
    result.a = fmin(255, (u8) (s * (f32)c.a));

    return result;
}

/**
 * Scalar multiplication for colors.
 */
inline IntColor operator*(IntColor c, f32 s) {
    IntColor result = {};

    result.r = fmin(255, (u8) (s * (f32)c.r));
    result.g = fmin(255, (u8) (s * (f32)c.g));
    result.b = fmin(255, (u8) (s * (f32)c.b));
    result.a = fmin(255, (u8) (s * (f32)c.a));

    return result;
}

#endif //RAYTRACING_COLOR_H
