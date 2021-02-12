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

/**
 * Sets up a scene with three balls, a floor, and a back-wall.
 */
void setupThreeBallScene(Scene *scene, Camera *camera, u32 width, u32 height) {
    // Setup camera
    camera->position = {0, 1, 0};
    camera->viewDirection = {0, 0, 1};
    camera->upVector = {0, 1, 0};

    camera->fovy = 90;
    camera->nearClippingDistance = 0.1f;
    camera->aspectRatio = (f32) width / (f32) height;
    camera->exposure = 1.0f;

    // Setup spheres
    scene->numSpheres = 3;
    scene->spheres = (SphereObject *) malloc(sizeof(SphereObject) * scene->numSpheres);
    scene->spheres[0].sphere.position = {0, 2, 10};
    scene->spheres[0].sphere.radius = 2;
    scene->spheres[0].material = &PBM_ROUGH_RED;

    scene->spheres[1].sphere.position = {8, 3, 9};
    scene->spheres[1].sphere.radius = 3;
    scene->spheres[1].material = &PBM_METALLIC_GREEN;

    scene->spheres[2].sphere.position = {-2, 1, 8};
    scene->spheres[2].sphere.radius = 1;
    scene->spheres[2].material = &PBM_SMOOTH_BLUE;

    // Setup triangles
    scene->numTriangles = 4;
    scene->triangles = (TriangleObject *) malloc(sizeof(TriangleObject) * scene->numTriangles);

    // floor

    scene->triangles[0].triangle.B = {-50, 0, -50};
    scene->triangles[0].triangle.A = { 50, 0, -50};
    scene->triangles[0].triangle.C = {-50, 0,  50};
    scene->triangles[0].material = &PBM_GRAY;

    scene->triangles[1].triangle.B = {-50, 0,  50};
    scene->triangles[1].triangle.A = { 50, 0, -50};
    scene->triangles[1].triangle.C = { 50, 0,  50};
    scene->triangles[1].material = &PBM_GRAY;

    // back wall

    scene->triangles[2].triangle.B = {-50, 20,  50};
    scene->triangles[2].triangle.A = {-50,  0,  50};
    scene->triangles[2].triangle.C = { 50,  0,  50};
    scene->triangles[2].material = &PBM_GRAY;

    scene->triangles[3].triangle.B = {-50, 20,  50};
    scene->triangles[3].triangle.A = { 50,  0,  50};
    scene->triangles[3].triangle.C = { 50, 20,  50};
    scene->triangles[3].material = &PBM_GRAY;

    // Setup lights
    scene->numLights = 2;
    scene->lights = (PointLight *) malloc(sizeof(PointLight) * scene->numLights);

    scene->lights[0].position = {-5,3, 5};
    scene->lights[0].color = {1, 1, 1, 0};

    scene->lights[1].position = {5,5, 4};
    scene->lights[1].color = {1, 1, 1, 0};

    // Setup models
    scene->numModels = 0;


    // Setup background color
    scene->backgroundColor = {0.05f, 0.05f, 0.05f, 1.0f};
}

/**
 * Sets up the studio scene.
 */
void setupDragonScene(Scene *scene, Camera *camera, u32 width, u32 height, Mesh *studio) {
    // Setup camera
    camera->position = {0, 2, 0.1f};
    camera->viewDirection = {0, 0, 1};
    camera->upVector = {0, 1, 0};

    camera->fovy = 70;
    camera->nearClippingDistance = 0.1f;
    camera->aspectRatio = (f32) width / (f32) height;
    camera->exposure = 0.5f;

    // Setup spheres
    scene->numSpheres = 2;
    scene->spheres = (SphereObject *) malloc(sizeof(SphereObject) * scene->numSpheres);

    scene->spheres[0].sphere.position = {3, 2.5, 5};
    scene->spheres[0].sphere.radius = 0.5;
    scene->spheres[0].material = &PBM_ROUGH_RED;

    scene->spheres[1].sphere.position = {-3, 2, 5};
    scene->spheres[1].sphere.radius = 1;
    scene->spheres[1].material = &PBM_METALLIC_GREEN;

    // Setup triangles
    scene->numTriangles = 0;

    // Setup lights
    scene->numLights = 4;
    scene->lights = (PointLight *) malloc(sizeof(PointLight) * scene->numLights);

    scene->lights[0].position = {-3,3, 1};
    scene->lights[0].color = {3, 3, 3, 0};

    scene->lights[1].position = {3,3, 1};
    scene->lights[1].color = {3, 3, 3, 0};

    scene->lights[2].position = {-3,3, 7};
    scene->lights[2].color = {2, 2, 2, 0};

    scene->lights[3].position = {3,3, 7};
    scene->lights[3].color = {2, 2, 2, 0};

    // Setup models
    scene->numModels = 1;
    scene->models = (Model *) malloc(sizeof(Model) * scene->numModels);

    scene->models[0].mesh = studio;
    scene->models[0].transform = Matrix4x4::identity();

    // Setup background color
    scene->backgroundColor = {0.05f, 0.05f, 0.05f, 1.0f};
}

int main() {
    // Seed the RNG
    srand(0);

    // Width and heigh of the screen
    u32 width = 1920;
    u32 height = 1080;

    Color *pixels = (Color *) malloc(sizeof(Color) * width *  height);

    // Load meshes for scenes
    Mesh dragon;
    loadMesh(&dragon, "../res/objects/dragon.obj");


    // Setup scene and camera
    Scene scene;
    Camera camera;

    //setupThreeBallScene(&scene, &camera, width, height);
    setupDragonScene(&scene, &camera, width, height, &dragon);

    // Build ray tracer from scene
    RayTracer rayTracer = createRayTracer(&scene);

    // Render using either multi- or single-threading
#ifdef USE_MULTI_THREADING
    threadedRender(width, height ,pixels, &rayTracer, camera);
#else
    render(width, height ,pixels, &rayTracer, camera);
#endif

    // Free scene

    freeScene(&scene);

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


