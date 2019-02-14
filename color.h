//
// Created by s150818 on 13-2-2019.
//

#ifndef RAYTRACING_COLOR_H
#define RAYTRACING_COLOR_H

#pragma pack(push, 1)
/**
 * ARGB formatted color.
 */
struct Color {
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

/**
 *  Mixes two colors by averaging their rgba values.
 */
inline Color mix(Color a, Color b) {
    Color result = {};

    result.r = ((u32) a.r + b.r) / 2;
    result.g = ((u32) a.g + b.g) / 2;
    result.b = ((u32) a.b + b.b) / 2;
    result.a = ((u32) a.a + b.a) / 2;

    return result;
}

/**
 * Computes the sum of colors a and b.
 */
inline Color operator+(Color a, Color b) {
    Color result = {};

    result.r = fmin(255, (u32)a.r + b.r);
    result.g = fmin(255, (u32)a.g + b.g);
    result.b = fmin(255, (u32)a.b + b.b);
    result.a = fmin(255, (u32)a.a + b.a);

    return result;
}

/**
 * Scalar multiplication for colors.
 */
inline Color operator*(f32 s, Color c) {
    Color result = {};

    result.r = fmin(255, (u8) (s * (f32)c.r));
    result.g = fmin(255, (u8) (s * (f32)c.g));
    result.b = fmin(255, (u8) (s * (f32)c.b));
    result.a = fmin(255, (u8) (s * (f32)c.a));

    return result;
}

/**
 * Scalar multiplication for colors.
 */
inline Color operator*(Color c, f32 s) {
    Color result = {};

    result.r = fmin(255, (u8) (s * (f32)c.r));
    result.g = fmin(255, (u8) (s * (f32)c.g));
    result.b = fmin(255, (u8) (s * (f32)c.b));
    result.a = fmin(255, (u8) (s * (f32)c.a));

    return result;
}

#endif //RAYTRACING_COLOR_H
