#include <malloc.h>
#include <cstdio>
#include "types.h"
#include "math_utils.h"
#include "color.h"
#include "ray.h"
#include "scene.h"
#include "file_io.h"


int main() {
    u32 width = 1920*2;
    u32 height = 1080*2;

    Color *pixels = (Color *) malloc(sizeof(Color) * width *  height);

    Camera camera = {};

    camera.position = {0, 1, 0};
    camera.viewDirection = {0, 0, 1};
    camera.upVector = {0, 1, 0};

    camera.fovy = 90;
    camera.nearClippingDistance = 0.1f;
    camera.aspectRatio = (f32) width / (f32) height;

    // Set up some dummy spheres

#define numTestSpheres 3
    SphereObject spheres[numTestSpheres] = {};
    spheres[0].sphere.position = {0, 2, 10};
    spheres[0].sphere.radius = 2;
    spheres[0].material = PBM_ROUGH_RED;

    spheres[1].sphere.position = {8, 3, 9};
    spheres[1].sphere.radius = 3;
    spheres[1].material = PBM_METALLIC_GREEN;

    spheres[2].sphere.position = {-2, 1, 8};
    spheres[2].sphere.radius = 1;
    spheres[2].material = PBM_SMOOTH_BLUE;

#define numTestTriangles 4
    TriangleObject triangles[numTestTriangles];

    // floor

    triangles[0].triangle.B = {-50, 0, -50};
    triangles[0].triangle.A = { 50, 0, -50};
    triangles[0].triangle.C = {-50, 0,  50};
    triangles[0].material = PBM_GRAY;

    triangles[1].triangle.B = {-50, 0,  50};
    triangles[1].triangle.A = { 50, 0, -50};
    triangles[1].triangle.C = { 50, 0,  50};
    triangles[1].material = PBM_GRAY;

    // back wall

    triangles[2].triangle.B = {-50, 20,  50};
    triangles[2].triangle.A = {-50,  0,  50};
    triangles[2].triangle.C = { 50,  0,  50};
    triangles[2].material = PBM_GRAY;

    triangles[3].triangle.B = {-50, 20,  50};
    triangles[3].triangle.A = { 50,  0,  50};
    triangles[3].triangle.C = { 50, 20,  50};
    triangles[3].material = PBM_GRAY;

#define numTestLights 2
    PointLight lights[numTestLights] = {};
    lights[0].position = {-5,3, 5};
    lights[1].position = {5,5, 4};

    Scene scene = {};
    scene.backgroundColor = {0.05f, 0.05f, 0.05f, 1};
    scene.spheres = spheres;
    scene.numSpheres = numTestSpheres;
    scene.triangles = triangles;
    scene.numTriangles = numTestTriangles;
    scene.lights = lights;
    scene.numLights = numTestLights;

    for(u32 y = 0; y < height; y++) {
        for(u32 x = 0; x < width; x++) {

            Vector3D wPixel = pixelToWorldSpace(camera, x, y, width, height);

            Vector3D rayP = wPixel;
            Vector3D rayD = (wPixel - camera.position).normalized();

            Ray ray = {};
            ray.start = rayP;
            ray.direction = rayD;

            pixels[x + y * width] = traceThroughScene(ray, scene);
        }
    }

    // Make shift anti-aliasing


    u32 halfWidth = width / 2;
    u32 halfHeight = height / 2;

    Color *aaPixels = (Color *) malloc(sizeof(Color) * halfWidth * halfHeight);

    for(u32 xx = 0; xx < halfWidth; xx++) {
        for(u32 yy = 0; yy < halfHeight; yy++) {
            u32 x = 2*xx;
            u32 y = 2*yy;

            f32 r, g, b, a;

            r = pixels[x + y * width].r
                    +pixels[(x+1) + y * width].r
                    +pixels[x + (y+1) * width].r
                    +pixels[(x+1) + (y+1) * width].r;

            g = pixels[x + y * width].g
                +pixels[(x+1) + y * width].g
                +pixels[x + (y+1) * width].g
                +pixels[(x+1) + (y+1) * width].g;

            b = pixels[x + y * width].b
                +pixels[(x+1) + y * width].b
                +pixels[x + (y+1) * width].b
                +pixels[(x+1) + (y+1) * width].b;

            a = pixels[x + y * width].a
                +pixels[(x+1) + y * width].a
                +pixels[x + (y+1) * width].a
                +pixels[(x+1) + (y+1) * width].a;

            aaPixels[xx + yy * halfWidth].r = r /4;
            aaPixels[xx + yy * halfWidth].g = g /4;
            aaPixels[xx + yy * halfWidth].b = b /4;
            aaPixels[xx + yy * halfWidth].a = a /4;
        }
    }

    // Convert color format to that of BMP

    IntColor *bmpPixels = (IntColor *) malloc(sizeof(IntColor) * width * height);
    IntColor *bmpAaPixels = (IntColor *) malloc(sizeof(IntColor) * halfWidth * halfHeight);

    for(u32 y = 0; y < height; y++) {
        for (u32 x = 0; x < width; x++) {
            bmpPixels[x + y * width] = toIntColor(pixels[x + y * width]);
        }
    }

    for (u32 y = 0; y < halfHeight; y++) {
        for(u32 x = 0; x < halfWidth; x++) {
            bmpAaPixels[x + y * halfWidth] = toIntColor(aaPixels[x + y * halfWidth]);
        }
    }


    saveOutBitmap(createBitmap(bmpPixels, width, height), "P:/raytracing/res/reaytracing_test.bmp");
    saveOutBitmap(createBitmap(bmpAaPixels, halfWidth, halfHeight), "P:/raytracing/res/aa_raytracing_test.bmp");

    free(pixels);
    free(aaPixels);

    free(bmpPixels);
    free(bmpAaPixels);


    return 0;
}

