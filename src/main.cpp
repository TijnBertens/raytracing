#include <malloc.h>
#include <random>
#include <thread>
#include "types.h"
#include "math_utils.h"
#include "color.h"
#include "ray.h"
#include "scene.h"
#include "file_io.h"
#include "ray_tracer.h"

#define USE_MULTI_THREADING
#define NUM_THREADS 8
#define TILE_SIZE 32

/**
 * Function called by each thread. Image is divided into tiles, which are divided over the threads.
 */
void threadWork(u32 threadID, u32 imageWidth, u32 imageHeight, Color *pixels, RayTracer *tracer, Camera camera) {
    // Calculate number of tiles in total / horizontal direction / vertical direction
    u32 numTilesHor = (imageWidth + TILE_SIZE - 1) / TILE_SIZE;
    u32 numTilesVert = (imageHeight + TILE_SIZE - 1) / TILE_SIZE;

    u32 totalNumTiles = numTilesHor * numTilesVert;

    // Loop over tiles
    for(u32 tileID = threadID; tileID < totalNumTiles; tileID += NUM_THREADS) {
        u32 startX = (tileID % numTilesHor) * TILE_SIZE;
        u32 startY = (tileID / numTilesHor) * TILE_SIZE;

        printf("Thread %u starting tile %u / %u...\n", threadID, (tileID + 1), totalNumTiles);

        // Process pixels in tile
        for(u32 y = startY; y < startY + TILE_SIZE && y < imageHeight; y++) {

            for(u32 x = startX; x < startX + TILE_SIZE && x < imageWidth; x++) {
                // Number of samples per pixel
                const u32 SAMPLES = 15;

                pixels[x + y * imageWidth] = {0, 0, 0, 0};
                for(u32 i = 0; i < SAMPLES; i++) {
                    Vec3f wPixel = pixelToWorldSpaceRand(camera, x, y, imageWidth, imageHeight);

                    Vec3f rayP = wPixel;
                    Vec3f rayD = (wPixel - camera.position).normalized();

                    Ray ray = {};
                    ray.start = rayP;
                    ray.direction = rayD;

                    TraceReport traceReport = traceRay(tracer, ray);

                    pixels[x + y * imageWidth] = pixels[x + y * imageWidth] + (1.0 / SAMPLES) * traceReport.color;
                }
            }
        }
    }
}

/**
 * Renders the output image using multi-threading. Instantiates NUM_THREADS threads.
 */
void threadedRender(u32 imageWidth, u32 imageHeight, Color *pixels, RayTracer *tracer, Camera camera) {
    // Array to keep track of threads
    std::thread threads[NUM_THREADS];

    // Instantiate all threads
    for(u32 i = 0; i < NUM_THREADS; i++) {
        threads[i] = std::thread(threadWork, i, imageWidth, imageHeight, pixels, tracer, camera);
    }

    // Wait for all threads to finish
    for(u32 i = 0; i < NUM_THREADS; i++) {
        threads[i].join();
    }
}

/**
 * Renders the output image without mulit-threading.
 */
void render(u32 imageWidth, u32 imageHeight, Color *pixels, RayTracer *rayTracer, Camera camera) {
    // Process all pixels
    for(u32 y = 0; y < imageHeight; y++) {

        printf("Reached pixel line %u.\n", y);

        for(u32 x = 0; x < imageWidth; x++) {
            // Number of samples per pixel
            const u32 SAMPLES = 15;

            pixels[x + y * imageWidth] = {0, 0, 0, 0};
            for(u32 i = 0; i < SAMPLES; i++) {
                Vec3f wPixel = pixelToWorldSpaceRand(camera, x, y, imageWidth, imageHeight);

                Vec3f rayP = wPixel;
                Vec3f rayD = (wPixel - camera.position).normalized();

                Ray ray = {};
                ray.start = rayP;
                ray.direction = rayD;

                TraceReport traceReport = traceRay(rayTracer, ray);

                pixels[x + y * imageWidth] = pixels[x + y * imageWidth] + (1.0 / SAMPLES) * traceReport.color;
            }
        }
    }
}

int main() {
    // Seed the RNG
    srand(0);

    // Width and heigh of the screen
    u32 width = 1920;
    u32 height = 1080;

    Color *pixels = (Color *) malloc(sizeof(Color) * width *  height);

    // Set up the camera

    Camera camera = {};

    camera.position = {0, 1, 0};
    camera.viewDirection = {0, 0, 1};
    camera.upVector = {0, 1, 0};

    camera.fovy = 90;
    camera.nearClippingDistance = 0.1f;
    camera.aspectRatio = (f32) width / (f32) height;
    camera.exposure = 1.0f;

    // Set up some dummy spheres

#define numTestSpheres 3
    SphereObject spheres[numTestSpheres] = {};
    spheres[0].sphere.position = {0, 2, 10};
    spheres[0].sphere.radius = 2;
    spheres[0].material = &PBM_ROUGH_RED;

    spheres[1].sphere.position = {8, 3, 9};
    spheres[1].sphere.radius = 3;
    spheres[1].material = &PBM_METALLIC_GREEN;

    spheres[2].sphere.position = {-2, 1, 8};
    spheres[2].sphere.radius = 1;
    spheres[2].material = &PBM_SMOOTH_BLUE;

    // Set up some dummy triangles

#define numTestTriangles 4
    TriangleObject triangles[numTestTriangles];

    // floor

    triangles[0].triangle.B = {-50, 0, -50};
    triangles[0].triangle.A = { 50, 0, -50};
    triangles[0].triangle.C = {-50, 0,  50};
    triangles[0].material = &PBM_GRAY;

    triangles[1].triangle.B = {-50, 0,  50};
    triangles[1].triangle.A = { 50, 0, -50};
    triangles[1].triangle.C = { 50, 0,  50};
    triangles[1].material = &PBM_GRAY;

    // back wall

    triangles[2].triangle.B = {-50, 20,  50};
    triangles[2].triangle.A = {-50,  0,  50};
    triangles[2].triangle.C = { 50,  0,  50};
    triangles[2].material = &PBM_GRAY;

    triangles[3].triangle.B = {-50, 20,  50};
    triangles[3].triangle.A = { 50,  0,  50};
    triangles[3].triangle.C = { 50, 20,  50};
    triangles[3].material = &PBM_GRAY;

    // Set up some dummy lights

#define numTestLights 2
    PointLight lights[numTestLights] = {};
    lights[0].position = {-5,3, 5};
    lights[0].color = {1, 1, 1, 0};

    lights[1].position = {5,5, 4};
    lights[1].color = {1, 1, 1, 0};

    // Load dummy mesh

    Mesh nol;
    printf("loading mesh\n");
    loadMesh(&nol, "../res/objects/nol.obj");
    printf("loaded mesh\n");

    // Set up dummy models

#define numTestModels 1
    Model models[numTestModels] = {};
    models[0].mesh = &nol;
    models[0].transform =
            Matrix4x4::translation({0.4, 2, 4})
            * Matrix4x4::rotationY(-135)
            * Matrix4x4::scale({.5f, .5f, .5f});

    // Build scene

    Scene scene = {};
    scene.backgroundColor = {0.05f, 0.05f, 0.05f, 1};
    scene.spheres = spheres;
    scene.numSpheres = numTestSpheres;
    scene.triangles = triangles;
    scene.numTriangles = numTestTriangles;
    scene.models = models;
    scene.numModels = numTestModels;
    scene.lights = lights;
    scene.numLights = numTestLights;

    // Build ray tracer from scene
    RayTracer rayTracer = createRayTracer(&scene);

#ifdef USE_MULTI_THREADING
    threadedRender(width, height ,pixels, &rayTracer, camera);
#else
    render(width, height ,pixels, &rayTracer, camera);
#endif

    // Tone mapping

    for(u32 y = 0; y < height; y++) {
        for(u32 x = 0; x < width; x++) {
            pixels[x + y * width].r = 1.0 - exp(-pixels[x + y * width].r * camera.exposure);
            pixels[x + y * width].g = 1.0 - exp(-pixels[x + y * width].g * camera.exposure);
            pixels[x + y * width].b = 1.0 - exp(-pixels[x + y * width].b * camera.exposure);
        }
    }

    // Convert color format to that of BMP

    IntColor *bmpPixels = (IntColor *) malloc(sizeof(IntColor) * width * height);

    for(u32 y = 0; y < height; y++) {
        for (u32 x = 0; x < width; x++) {
            bmpPixels[x + y * width] = toIntColor(pixels[x + y * width]);
        }
    }

    // Save bitmap

    saveOutBitmap(createBitmap(bmpPixels, width, height), "../out/output.bmp");

    // Clean up allocated pixel arrays

    free(pixels);
    free(bmpPixels);

    return 0;
}


