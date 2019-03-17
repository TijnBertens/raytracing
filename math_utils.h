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

struct Vector2D {
    union {
        struct {
            f32 x;
            f32 y;
        };
        f32 values[2];
    };

public:
    inline f32 dot(Vector2D b);
    inline f32 length();

    inline Vector2D normalized();
};

struct Vector3D {
    union {
        struct {
            f32 x;
            f32 y;
            f32 z;
        };
        f32 values[3];
    };

public:
    inline f32 dot(Vector3D b);
    inline f32 length();
    inline f32 lengthSq();

    inline Vector3D cross(Vector3D b);
    inline Vector3D normalized();
    inline Vector3D reflectIn(Vector3D b);
};

struct Vector4D {
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
    inline f32 dot(Vector4D b);
    inline f32 length();

    inline Vector4D normalized();
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

    inline static Matrix4x4 scale(Vector3D s);
    inline static Matrix4x4 translation(Vector3D t);
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
        Vector3D v;
    };

public:
    inline Quaternion inverse();
    inline Vector3D rotateVector(Vector3D);

    inline static Quaternion rotation(f32 angle, Vector3D axis);
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

inline f32 Vector2D::dot(Vector2D b) {
    return x*b.x + y*b.y;
}


inline f32 Vector2D::length() {
    return sqrt(x*x + y*y);
}

inline Vector2D Vector2D::normalized() {
    Vector2D result = {};

    f32 length = this->length();
    result.x = this->x / length;
    result.y = this->y / length;

    return result;
}

/*
 * OPERATORS
 */

inline Vector2D operator+(Vector2D a, Vector2D b) {
    Vector2D result = {};

    result.x = a.x + b.x;
    result.y = a.y + b.y;

    return result;
}

inline Vector2D operator-(Vector2D a, Vector2D b) {
    Vector2D result = {};

    result.x = a.x - b.x;
    result.y = a.y - b.y;

    return result;
}

inline Vector2D operator-(Vector2D a) {
    Vector2D result = {};

    result.x = - a.x;
    result.y = - a.y;

    return result;
}

inline Vector2D operator*(f32 f, Vector2D v) {
    Vector2D result = {};

    result.x = v.x * f;
    result.y = v.y * f;

    return result;
}

inline Vector2D operator*(Vector2D v, f32 f) {
    Vector2D result = {};

    result.x = v.x * f;
    result.y = v.y * f;

    return result;
}

inline Vector2D operator/(f32 f, Vector2D v) {
    Vector2D result = {};

    result.x = v.x / f;
    result.y = v.y / f;

    return result;
}

inline Vector2D operator/(Vector2D v, f32 f) {
    Vector2D result = {};

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

inline Vector3D operator+(Vector3D a, Vector3D b) {
    Vector3D result = {};

    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;

    return result;
}

inline Vector3D operator-(Vector3D a, Vector3D b) {
    Vector3D result = {};

    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;

    return result;
}

inline Vector3D operator-(Vector3D a) {
    Vector3D result = {};

    result.x = - a.x;
    result.y = - a.y;
    result.z = - a.z;

    return result;
}

inline Vector3D operator*(f32 f, Vector3D v) {
    Vector3D result = {};

    result.x = v.x * f;
    result.y = v.y * f;
    result.z = v.z * f;

    return result;
}

inline Vector3D operator*(Vector3D v, f32 f) {
    Vector3D result = {};

    result.x = v.x * f;
    result.y = v.y * f;
    result.z = v.z * f;

    return result;
}

inline Vector3D operator/(f32 f, Vector3D v) {
    Vector3D result = {};

    result.x = v.x / f;
    result.y = v.y / f;
    result.z = v.z / f;

    return result;
}

inline Vector3D operator/(Vector3D v, f32 f) {
    Vector3D result = {};

    result.x = v.x / f;
    result.y = v.y / f;
    result.z = v.z / f;

    return result;
}

/*
 * Methods
 */

inline f32 Vector3D::dot(Vector3D b) {
    return x*b.x + y*b.y + z*b.z;
}

inline Vector3D Vector3D::cross(Vector3D b) {
    Vector3D result = {};

    result.x = y*b.z - z*b.y;
    result.y = z*b.x - x*b.z;
    result.z = x*b.y - y*b.x;

    return result;
}

inline f32 Vector3D::length() {
    return sqrt(x*x + y*y + z*z);
}

inline f32 Vector3D::lengthSq() {
    return x*x + y*y + z*z;
}

inline Vector3D Vector3D::normalized() {
    Vector3D result = {};

    f32 length = this->length();
    result.x = this->x / length;
    result.y = this->y / length;
    result.z = this->z / length;

    return result;
}

inline Vector3D Vector3D::reflectIn(Vector3D b) {
    Vector3D result = {};

    result = *this - 2.0f*(this->dot(b))*b;

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

inline f32 Vector4D::dot(Vector4D b) {
    return x*b.x + y*b.y + z*b.z + w*b.w;
}

inline f32 Vector4D::length() {
    return sqrt(x*x + y*y + z*z + w*w);
}

inline Vector4D Vector4D::normalized() {
    Vector4D result = {};

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

inline Vector4D operator+(Vector4D a, Vector4D b) {
    Vector4D result = {};

    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;
    result.w = a.w + b.w;

    return result;
}

inline Vector4D operator-(Vector4D a, Vector4D b) {
    Vector4D result = {};

    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;
    result.w = a.w - b.w;

    return result;
}

inline Vector4D operator-(Vector4D a) {
    Vector4D result = {};

    result.x = -a.x;
    result.y = -a.y;
    result.z = -a.z;
    result.w = -a.w;

    return result;
}

inline Vector4D operator*(f32 f, Vector4D v) {
    Vector4D result = {};

    result.x = v.x * f;
    result.y = v.y * f;
    result.z = v.z * f;
    result.w = v.w * f;

    return result;
}

inline Vector4D operator*(Vector4D v, f32 f) {
    Vector4D result = {};

    result.x = v.x * f;
    result.y = v.y * f;
    result.z = v.z * f;
    result.w = v.w * f;

    return result;
}

inline Vector4D operator/(f32 f, Vector4D v) {
    Vector4D result = {};

    result.x = v.x / f;
    result.y = v.y / f;
    result.z = v.z / f;
    result.w = v.w / f;

    return result;
}

inline Vector4D operator/(Vector4D v, f32 f) {
    Vector4D result = {};

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

inline Matrix4x4 Matrix4x4::scale(Vector3D s) {
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

inline Matrix4x4 Matrix4x4::translation(Vector3D t) {
    f32 x = t.x;
    f32 y = t.y;
    f32 z = t.z;
    Matrix4x4 result = {
            {{
                     { 0, 0, 0, 0},
                     { 0, 0, 0, 0},
                     { 0, 0, 0, 0},
                     { x, y, z, 1}}
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

inline Vector3D Quaternion::rotateVector(Vector3D vec) {
    Quaternion qv = {};
    qv.w == 0;
    qv.v = vec;

    return (*this * qv * this->inverse()).v;
}

inline Quaternion Quaternion::rotation(f32 degrees, Vector3D axis) {
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
