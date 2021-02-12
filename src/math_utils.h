//
// Created by Tijn Bertens on 31-10-2018.
//

#ifndef CLIONPROJECT_MATH_UTILS_H
#define CLIONPROJECT_MATH_UTILS_H

#include "math.h"
#include "types.h"

/*
 * DEFINITIONS
 */

#define PI_32 3.14159265359

struct Vec2f;
struct Vec3f;
struct Vec4f;

struct Vec2f {
    union {
        struct {
            f32 x;
            f32 y;
        };
        f32 values[2];
    };

public:
    inline f32 dot(Vec2f b);
    inline f32 length();

    inline Vec2f normalized();
};

struct Vec3f {
    union {
        struct {
            f32 x;
            f32 y;
            f32 z;
        };
        f32 values[3];
    };

public:
    inline f32 dot(Vec3f b);
    inline f32 length();
    inline f32 lengthSq();

    inline Vec3f cross(Vec3f b);
    inline Vec3f normalized();
    inline Vec3f reflectIn(Vec3f b);

    static inline Vec3f min(Vec3f a, Vec3f b);
    static inline Vec3f max(Vec3f a, Vec3f b);

    static inline Vec4f toHomogeneous(Vec3f v, bool isDirection);
    static inline Vec3f fromHomogeneous(Vec4f v);
};

struct Vec4f {
    union {
        struct {
            f32 x;
            f32 y;
            f32 z;
            f32 w;
        };
        f32 values[4];
    };

public:
    inline f32 dot(Vec4f b);
    inline f32 length();

    inline Vec4f normalized();
};

struct Matrix4x4 {
    union {
        f32 val[4][4];

        struct {
            f32 m00;
            f32 m01;
            f32 m02;
            f32 m03;

            f32 m10;
            f32 m11;
            f32 m12;
            f32 m13;

            f32 m20;
            f32 m21;
            f32 m22;
            f32 m23;

            f32 m30;
            f32 m31;
            f32 m32;
            f32 m33;
        };
    };

public:
    inline Matrix4x4 transpose();

    inline static Matrix4x4 identity();

    inline static Matrix4x4 scale(Vec3f s);
    inline static Matrix4x4 translation(Vec3f t);

    inline static Matrix4x4 rotationX(f32 a);
    inline static Matrix4x4 rotationY(f32 a);
    inline static Matrix4x4 rotationZ(f32 a);

    static Matrix4x4 perspective(f32 aspectRatio, f32 fovy, f32 zNear, f32 zFar);
};

struct Quaternion {
    f32 w;
    union {
        struct {
            f32 x;
            f32 y;
            f32 z;
        };
        Vec3f v;
    };

public:
    inline Quaternion inverse();
    inline Vec3f rotateVector(Vec3f);

    inline static Quaternion rotation(f32 angle, Vec3f axis);
    inline static Quaternion slerp(Quaternion start, Quaternion end, f32 t);
};

/*
 * ------------------------------
 * -------- Vector 2D -----------
 * ------------------------------
 */

/*
 * Methods
 */

inline f32 Vec2f::dot(Vec2f b) {
    return x*b.x + y*b.y;
}


inline f32 Vec2f::length() {
    return sqrt(x*x + y*y);
}

inline Vec2f Vec2f::normalized() {
    Vec2f result = {};

    f32 length = this->length();
    result.x = this->x / length;
    result.y = this->y / length;

    return result;
}

/*
 * OPERATORS
 */

inline Vec2f operator+(Vec2f a, Vec2f b) {
    Vec2f result = {};

    result.x = a.x + b.x;
    result.y = a.y + b.y;

    return result;
}

inline Vec2f operator-(Vec2f a, Vec2f b) {
    Vec2f result = {};

    result.x = a.x - b.x;
    result.y = a.y - b.y;

    return result;
}

inline Vec2f operator-(Vec2f a) {
    Vec2f result = {};

    result.x = - a.x;
    result.y = - a.y;

    return result;
}

inline Vec2f operator*(f32 f, Vec2f v) {
    Vec2f result = {};

    result.x = v.x * f;
    result.y = v.y * f;

    return result;
}

inline Vec2f operator*(Vec2f v, f32 f) {
    Vec2f result = {};

    result.x = v.x * f;
    result.y = v.y * f;

    return result;
}

inline Vec2f operator/(f32 f, Vec2f v) {
    Vec2f result = {};

    result.x = v.x / f;
    result.y = v.y / f;

    return result;
}

inline Vec2f operator/(Vec2f v, f32 f) {
    Vec2f result = {};

    result.x = v.x / f;
    result.y = v.y / f;

    return result;
}

/*
 * ------------------------------
 * -------- Vector 3D -----------
 * ------------------------------
 */

/*
 * OPERATORS
 */

inline Vec3f operator+(Vec3f a, Vec3f b) {
    Vec3f result = {};

    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;

    return result;
}

inline Vec3f operator-(Vec3f a, Vec3f b) {
    Vec3f result = {};

    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;

    return result;
}

inline Vec3f operator-(Vec3f a) {
    Vec3f result = {};

    result.x = - a.x;
    result.y = - a.y;
    result.z = - a.z;

    return result;
}

inline Vec3f operator*(f32 f, Vec3f v) {
    Vec3f result = {};

    result.x = v.x * f;
    result.y = v.y * f;
    result.z = v.z * f;

    return result;
}

inline Vec3f operator*(Vec3f v, f32 f) {
    Vec3f result = {};

    result.x = v.x * f;
    result.y = v.y * f;
    result.z = v.z * f;

    return result;
}

inline Vec3f operator/(f32 f, Vec3f v) {
    Vec3f result = {};

    result.x = v.x / f;
    result.y = v.y / f;
    result.z = v.z / f;

    return result;
}

inline Vec3f operator/(Vec3f v, f32 f) {
    Vec3f result = {};

    result.x = v.x / f;
    result.y = v.y / f;
    result.z = v.z / f;

    return result;
}

/*
 * Methods
 */

inline f32 Vec3f::dot(Vec3f b) {
    return x*b.x + y*b.y + z*b.z;
}

inline Vec3f Vec3f::cross(Vec3f b) {
    Vec3f result = {};

    result.x = y*b.z - z*b.y;
    result.y = z*b.x - x*b.z;
    result.z = x*b.y - y*b.x;

    return result;
}

inline f32 Vec3f::length() {
    return sqrt(x*x + y*y + z*z);
}

inline f32 Vec3f::lengthSq() {
    return x*x + y*y + z*z;
}

inline Vec3f Vec3f::normalized() {
    Vec3f result = {};

    f32 length = this->length();
    result.x = this->x / length;
    result.y = this->y / length;
    result.z = this->z / length;

    return result;
}

inline Vec3f Vec3f::reflectIn(Vec3f b) {
    Vec3f result = {};

    result = *this - 2.0f*(this->dot(b))*b;

    return result;
}

inline Vec3f Vec3f::min(Vec3f a, Vec3f b) {
    Vec3f result = {};

    result.x = fminf(a.x, b.x);
    result.y = fminf(a.y, b.y);
    result.z = fminf(a.z, b.z);

    return result;
}

inline Vec3f Vec3f::max(Vec3f a, Vec3f b) {
    Vec3f result = {};

    result.x = fmaxf(a.x, b.x);
    result.y = fmaxf(a.y, b.y);
    result.z = fmaxf(a.z, b.z);

    return result;
}

inline Vec4f Vec3f::toHomogeneous(const Vec3f v, const bool isDirection) {
    Vec4f result = {v.x, v.y, v.z, isDirection ? 0.f : 1.f};
    return result;
}

inline Vec3f Vec3f::fromHomogeneous(Vec4f v) {
    Vec3f result = {v.x, v.y, v.z};
    return result;
}

/*
 * ------------------------------
 * -------- Vector 4D -----------
 * ------------------------------
 */

/*
 * Methods
 */

inline f32 Vec4f::dot(Vec4f b) {
    return x*b.x + y*b.y + z*b.z + w*b.w;
}

inline f32 Vec4f::length() {
    return sqrt(x*x + y*y + z*z + w*w);
}

inline Vec4f Vec4f::normalized() {
    Vec4f result = {};

    f32 length = this->length();
    result.x = this->x / length;
    result.y = this->y / length;
    result.z = this->z / length;
    result.w = this->w / length;

    return result;
}

/*
 * OPERATORS
 */

inline Vec4f operator+(Vec4f a, Vec4f b) {
    Vec4f result = {};

    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;
    result.w = a.w + b.w;

    return result;
}

inline Vec4f operator-(Vec4f a, Vec4f b) {
    Vec4f result = {};

    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;
    result.w = a.w - b.w;

    return result;
}

inline Vec4f operator-(Vec4f a) {
    Vec4f result = {};

    result.x = -a.x;
    result.y = -a.y;
    result.z = -a.z;
    result.w = -a.w;

    return result;
}

inline Vec4f operator*(f32 f, Vec4f v) {
    Vec4f result = {};

    result.x = v.x * f;
    result.y = v.y * f;
    result.z = v.z * f;
    result.w = v.w * f;

    return result;
}

inline Vec4f operator*(Vec4f v, f32 f) {
    Vec4f result = {};

    result.x = v.x * f;
    result.y = v.y * f;
    result.z = v.z * f;
    result.w = v.w * f;

    return result;
}

inline Vec4f operator/(f32 f, Vec4f v) {
    Vec4f result = {};

    result.x = v.x / f;
    result.y = v.y / f;
    result.z = v.z / f;
    result.w = v.w / f;

    return result;
}

inline Vec4f operator/(Vec4f v, f32 f) {
    Vec4f result = {};

    result.x = v.x / f;
    result.y = v.y / f;
    result.z = v.z / f;
    result.w = v.w / f;

    return result;
}

/*
 * ------------------------------
 * -------- Matrix 4D -----------
 * ------------------------------
 */

/*
 * Methods
 */

inline Matrix4x4 Matrix4x4::transpose() {
    Matrix4x4 result = {};
    for(u8 row = 0; row < 4; row++) {
        for(u8 column = 0; column < 4; column++) {
            result.val[row][column] = this->val[column][row];
        }
    }
    return result;
}

inline Matrix4x4 Matrix4x4::identity() {
    Matrix4x4 result = {
            {{
                     { 1, 0, 0, 0},
                     { 0, 1, 0, 0},
                     { 0, 0, 1, 0},
                     { 0, 0, 0, 1}}
            }};
    return result;
}

inline Matrix4x4 Matrix4x4::scale(Vec3f s) {
    f32 a = s.x;
    f32 b = s.y;
    f32 c = s.z;
    Matrix4x4 result = {
           {{
              { a, 0, 0, 0},
              { 0, b, 0, 0},
              { 0, 0, c, 0},
              { 0, 0, 0, 1}}
           }};
    return result;
}

inline Matrix4x4 Matrix4x4::translation(Vec3f t) {
    f32 x = t.x;
    f32 y = t.y;
    f32 z = t.z;
    Matrix4x4 result = {
            {{
                     { 1, 0, 0, 0},
                     { 0, 1, 0, 0},
                     { 0, 0, 1, 0},
                     { x, y, z, 1}}
            }};
    return result;
}

inline Matrix4x4 Matrix4x4::rotationX(const f32 a) {
    f32 rad = a * PI_32 / 180;
    f32 s = sin(rad);
    f32 c = cos(rad);
    Matrix4x4 result = {
            {{
                     {1, 0, 0, 0},
                     {0, c, s, 0},
                     {0, -s, c, 0},
                     {0, 0, 0, 1}}
            }};
    return result;
}

inline Matrix4x4 Matrix4x4::rotationY(const f32 a) {
    f32 rad = a * PI_32 / 180;
    f32 s = sin(rad);
    f32 c = cos(rad);
    Matrix4x4 result = {
            {{
                     {c, 0, -s, 0},
                     {0, 1, 0, 0},
                     {s, 0, c, 0},
                     {0, 0, 0, 1}}
            }};
    return result;
}

inline Matrix4x4 Matrix4x4::rotationZ(const f32 a) {
    f32 rad = a * PI_32 / 180;
    f32 s = sin(rad);
    f32 c = cos(rad);
    Matrix4x4 result = {
            {{
                     {c, s, 0, 0},
                     {-s, c, 0, 0},
                     {0, 0, 1, 0},
                     {0, 0, 0, 1}}
            }};
    return result;
}

Matrix4x4 Matrix4x4::perspective(f32 aspectRatio, f32 fovy, f32 zNear, f32 zFar) {
    f32 aa = 1.0 / tan((fovy*PI_32 / 180) / 2);
    f32 a = aa / aspectRatio;
    f32 b = aa;
    f32 c = (zFar + zNear) / (zNear - zFar);
    f32 d = (2 * zFar * zNear) / (zNear - zFar);
    Matrix4x4 result = {
            {{
                {a, 0, 0, 0},
                {0, b, 0, 0},
                {0, 0, c,-1},
                {0, 0, d, 0}
            }}
    };
    return result;
}

/*
 * Operators
 */

inline Matrix4x4 operator-(Matrix4x4 A, Matrix4x4 B) {
    Matrix4x4 result = {};
    for(u8 column = 0; column < 4; column++) {
        for(u8 row = 0; row < 4; row++) {
            result.val[column][row] = A.val[column][row] - B.val[column][row];
        }
    }
    return result;
}

inline Matrix4x4 operator+(Matrix4x4 A, Matrix4x4 B) {
    Matrix4x4 result = {};
    for(u8 column = 0; column < 4; column++) {
        for(u8 row = 0; row < 4; row++) {
            result.val[column][row] = A.val[column][row] + B.val[column][row];
        }
    }
    return result;
}

inline Matrix4x4 operator*(Matrix4x4 A, Matrix4x4 B) {
    Matrix4x4 result = {};
    for(u8 column = 0; column < 4; column++) {
        for(u8 row = 0; row < 4; row++) {
            for(u8 i = 0; i < 4; i++) {
                result.val[column][row] += A.val[i][row] * B.val[column][i];
            }
        }
    }
    return result;
}

inline Vec4f operator*(Matrix4x4 A, Vec4f B) {
    Vec4f result = {};
    for(u8 column = 0; column < 4; column++) {
        for(u8 row = 0; row < 4; row++) {
            result.values[row] +=  A.val[column][row] * B.values[column];
        }
    }
    return result;
}

/*
 * ------------------------------
 * ------- Quaternions ----------
 * ------------------------------
 */

/*
 * OPERATORS
 */

inline Quaternion operator*(Quaternion a, Quaternion b) {
    Quaternion result = {};

    result.w = b.w * a.w - b.v.dot(a.v);
    result.v = a.v * b.w + b.v * a.w + a.v.cross(b.v);

    return result;
}

inline Quaternion operator^(Quaternion a, f32 t) {
    return Quaternion::rotation(a.w * t, a.v);
}

/*
 * Methods
 */

inline Quaternion Quaternion::inverse() {
    Quaternion result;

    result.w = this->w;
    result.v = this->v * -1;

    return result;
}

inline Vec3f Quaternion::rotateVector(Vec3f vec) {
    Quaternion qv = {};
    qv.w == 0;
    qv.v = vec;

    return (*this * qv * this->inverse()).v;
}

inline Quaternion Quaternion::rotation(f32 degrees, Vec3f axis) {
    Quaternion result = {};

    f32 radians = (degrees * PI_32) / 180.0f;

    result.w = cos(radians / 2);

    f32 s = sin(radians / 2);
    result.x = s * axis.x;
    result.y = s * axis.y;
    result.z = s * axis.z;

    return result;
}

inline Quaternion Quaternion::slerp(Quaternion start, Quaternion end, f32 t) {
    // ((end * start.Inverted()) ^ t) * start;
    return ((end * start.inverse()) ^ t) * start;
}


/*
 * ------------------------------
 * ----- General Utilities ------
 * ------------------------------
 */

inline f32 clamp(f32 n, f32 a, f32 b) {
    return (fmin(fmax(n, a), b));
}


#endif //CLIONPROJECT_MATH_UTILS_H
