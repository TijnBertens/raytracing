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
        Vec4f v;
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
 *  Linear interpolation between color a and b, by some mixFactor [0.0, 1.0].
 */
inline Color mix(Color a, Color b, f32 mixFactor) {
    Color result = {};

    result.r = a.r + mixFactor * (b.r - a.r);
    result.g = a.g + mixFactor * (b.g - a.g);
    result.b = a.b + mixFactor * (b.b - a.b);
    result.a = a.a + mixFactor * (b.a - a.a);

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
 * Computes the difference of colors a and b.
 */
inline Color operator-(Color a, Color b) {
    Color result = {};

    result.r = fmax(0, a.r - b.r);
    result.g = fmax(0, a.g - b.g);
    result.b = fmax(0, a.b - b.b);
    result.a = fmax(0, a.a - b.a);

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

/**
 * Component wise color-color multiplication for colors.
 */
inline Color operator*(Color a, Color b) {
    Color result = {};

    result.r = a.r * b.r;
    result.g = a.g * b.g;
    result.b = a.b * b.b;
    result.a = a.a * b.a;

    return result;
}

/**
 * Scalar multiplication for colors.
 */
inline Color operator/(f32 s, Color c) {
    Color result = {};

    result.r = clamp(c.r / s, 0, 1);
    result.g = clamp(c.g / s, 0, 1);
    result.b = clamp(c.b / s, 0, 1);
    result.a = clamp(c.a / s, 0, 1);

    return result;
}

/**
 * Scalar multiplication for colors.
 */
inline Color operator/(Color c, f32 s) {
    Color result = {};

    result.r = clamp(c.r / s, 0, 1);
    result.g = clamp(c.g / s, 0, 1);
    result.b = clamp(c.b / s, 0, 1);
    result.a = clamp(c.a / s, 0, 1);

    return result;
}

/**
 * For each component of a color, take the max between that component and n.
 */
inline Color cmax(Color c, f32 n) {
    Color result = {};

    result.r = fmax(c.r, n);
    result.g = fmax(c.g, n);
    result.b = fmax(c.b, n);
    result.a = fmax(c.a, n);

    return result;
}

/**
 * For each component of a color, take the min between that component and n.
 */
inline Color cmin(Color c, f32 n) {
    Color result = {};

    result.r = fmin(c.r, n);
    result.g = fmin(c.g, n);
    result.b = fmin(c.b, n);
    result.a = fmin(c.a, n);

    return result;
}

/**
 * For each component of a color, take the max between that component and n.
 */
inline Color cmax(Color c, Color n) {
    Color result = {};

    result.r = fmax(c.r, n.r);
    result.g = fmax(c.g, n.g);
    result.b = fmax(c.b, n.b);
    result.a = fmax(c.a, n.a);

    return result;
}

/**
 * For each component of a color, take the min between that component and n.
 */
inline Color cmin(Color c, Color n) {
    Color result = {};

    result.r = fmin(c.r, n.r);
    result.g = fmin(c.g, n.g);
    result.b = fmin(c.b, n.b);
    result.a = fmin(c.a, n.a);

    return result;
}

//----------------------
//------IntColor--------
//----------------------

/**
 * Convert to 32bit IntColor.
 */
inline Color toColor(IntColor c) {
    Color result;

    result.r = c.r / 255.0;
    result.g = c.g / 255.0;
    result.b = c.b / 255.0;
    result.a = c.a / 255.0;

    return result;
}

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
