#include <malloc.h>
#include <cstdio>
#include "types.h"
#include "math_utils.h"
#include "color.h"
#include "ray.h"
#include "scene.h"

/**
 * It is important that the bitmap header is not padded!
 */
#pragma pack(push, 1)

/**
 * The header of the .bmp file format.
 */
struct BitmapHeader {
    u16 fileType;
    u32 fileSize;
    u16 reserved0;
    u16 reserved1;
    u32 dataOffset;

    u32 headerSize;
    u32 width;
    u32 height;
    u16 numPlanes;
    u16 bitsPerPixel;
    u32 compression;
    u32 dataSize;
    s32 horzResolution;
    s32 vertResolution;
    u32 numColorsInPalette;
    u32 importantColors;
};
#pragma pack(pop)

/**
 * An internal representation of a bitmap file, containing a .bmp header and color data.
 *
 * @note: There is no extra metadata in this struct, since everything can be found in/constructed from data in the header.
 */
struct Bitmap {
    BitmapHeader header;
    IntColor *data;
};

/**
 * Saves out a bitmap to a specified file path.
 */
bool saveOutBitmap(Bitmap bitmap, const char *filePath) {
    FILE *outFile = fopen(filePath, "wb");

    if(outFile) {
        fwrite(&bitmap, sizeof(bitmap.header), 1, outFile);
        fwrite(bitmap.data, bitmap.header.dataSize, 1, outFile);
        fclose(outFile);
        return true;
    } else {
        printf("Failed to save bitmap file, could not open file.");
        return false;
    }
}

/**
 * Creates a bitmap with the given data, of the specified width and height. The header will automatically be filled in
 * with the standard values used in this program.
 */
Bitmap createBitmap(IntColor *data, u32 width, u32 height) {
    Bitmap result = {};

    result.header.fileType = 0x4D42;                            // 'BM'
    result.header.fileSize = sizeof(BitmapHeader) + (width * height * sizeof(u32));
    result.header.dataOffset = sizeof(BitmapHeader);
    result.header.headerSize = sizeof(BitmapHeader) - 14; // 14 is the size of the file type/size, offset + reserved bytes
    result.header.width = width;
    result.header.height = height;
    result.header.numPlanes = 1;
    result.header.bitsPerPixel = 32;                            // RGBA format
    result.header.compression = 0;                              // no compression
    result.header.dataSize = (width * height * sizeof(u32));    // again, 32 bits per pixels for RGBA
    result.header.horzResolution = 5000;                        // arbitrary....
    result.header.vertResolution = 5000;                        // arbitrary....
    result.header.numColorsInPalette = 0;                       // we don't use a palette
    result.header.importantColors = 0;

    result.data = data;

    return result;
}

/**
 * Representation of a camera object.
 */
struct Camera {
    Vector3D position;          // world space position
    Vector3D viewDirection;     // world space view direction (normalized)
    Vector3D upVector;          // world space up vector (normalized)

    f32 fovy;                   // vertical field of view
    f32 nearClippingDistance;   // distance from camera to screen
    f32 aspectRatio;            // width over height
};

/**
 * Calculates the world space position of a pixel on the screen of a camera.
 */
Vector3D pixelToWorldSpace(Camera camera, u32 x, u32 y, u32 width, u32 height) {
    Vector3D screenCenter = camera.position + (camera.viewDirection * camera.nearClippingDistance);

    Vector3D screenVerticalDirection = camera.upVector;
    Vector3D screenHorizontalDirection = camera.upVector.cross(camera.viewDirection).normalized();

    f32 screenHeight = 2 * (f32) tan((camera.fovy / 2.0) * PI_32 / 180.0) * camera.nearClippingDistance;
    f32 screenWidth =  screenHeight * camera.aspectRatio;

    Vector3D botLeft =
            screenCenter
            - ((screenWidth / 2) * screenHorizontalDirection)
            - ((screenHeight / 2) * screenVerticalDirection);

    Vector3D botRight =
            screenCenter
            + ((screenWidth / 2) * screenHorizontalDirection)
            - ((screenHeight / 2) * screenVerticalDirection);

    Vector3D topLeft =
            screenCenter
            - ((screenWidth / 2) * screenHorizontalDirection)
            + ((screenHeight / 2) * screenVerticalDirection);

    f32 u = (f32) x / (f32) width;
    f32 v = (f32) y / (f32) height;

    Vector3D result =
            botLeft + (u * (botRight - botLeft)) + (v * (topLeft - botLeft));

    return result;
}

//  ------------------------------------------------------------------------------
// The following functions are implementations of Physically Based Shading functions
// all credits to: https://learnopengl.com/PBR/Lighting
//  ------------------------------------------------------------------------------

/**
 * Credits to:
 * https://learnopengl.com/PBR/Lighting
 */
Color fresnelSchlick(float cosTheta, Color F0)
{
    Color one = {1, 1, 1, 1};
    return F0 + (one - F0) * pow(1.0 - cosTheta, 5.0);
}

float DistributionGGX(Vector3D N, Vector3D H, f32 roughness)
{
    f32 a      = roughness*roughness;
    f32 a2     = a*a;
    f32 NdotH  = fmax(N.dot(H), 0.0);
    f32 NdotH2 = NdotH*NdotH;

    f32 num   = a2;
    f32 denom = (NdotH2 * (a2 - 1.0) + 1.0);
    denom = PI_32 * denom * denom;

    return num / denom;
}

f32 GeometrySchlickGGX(f32 NdotV, f32 roughness)
{
    f32 r = (roughness + 1.0);
    f32 k = (r*r) / 8.0;

    f32 num   = NdotV;
    f32 denom = NdotV * (1.0 - k) + k;

    return num / denom;
}
f32 GeometrySmith(Vector3D N, Vector3D V, Vector3D L, f32 roughness)
{
    f32 NdotV = fmax(N.dot(V), 0.0);
    f32 NdotL = fmax(N.dot(L), 0.0);
    f32 ggx2  = GeometrySchlickGGX(NdotV, roughness);
    f32 ggx1  = GeometrySchlickGGX(NdotL, roughness);

    return ggx1 * ggx2;
}

Color calculateSurfaceColorFromLightPBR(Vector3D viewPosition, Vector3D surfacePosition, PBMaterial surface, PointLight light, Vector3D normal) {

    Vector3D V = (viewPosition - surfacePosition).normalized();
    Vector3D N = normal;

    Color F0 = {0.3, 0.3, 0.3, 1};
    F0 = mix(F0, surface.albedo, surface.metallic);

    // calculate per-light radiance
    Vector3D L = (light.position - surfacePosition).normalized();
    Vector3D H = (V + L).normalized();
    f32 distance    = (light.position - surfacePosition).length();
    f32 attenuation = 1.0 / (distance * distance); //todo: hardcoded attenuation

    Color white = {1, 1, 1, 1};
    Color radiance     = white; //white * attenuation; //todo: hardcoded light color

    // cook-torrance brdf
    f32 NDF = DistributionGGX(N, H, surface.roughness);
    f32 G   = GeometrySmith(N, V, L, surface.roughness);
    Color F    = fresnelSchlick(fmax(H.dot(V), 0.0), F0);

    Color kS = F;
    Color kD = {1 - kS.r, 1 - kS.g, 1 - kS.b, 1}; //todo: what about alpha here?  //vec3(1.0) - kS; // here
    kD = kD * (1.0 - surface.metallic);

    Color numerator    = NDF * G * F;
    f32 denominator = 4.0 * fmax(N.dot(V), 0.0) * fmax(N.dot(L), 0.0);
    Color specular     = numerator / fmax(denominator, 0.001);

    // add to outgoing radiance Lo
    f32 NdotL = fmax(N.dot(L), 0.0);
    return (kD * surface.albedo / PI_32 + specular) * radiance * NdotL;
}

//  ------------------------------------------------------------------------------
//  ------------------------------------------------------------------------------

Color calculateSurfaceColorFromLight(Vector3D viewPosition, Vector3D surfacePosition, Material surface, PointLight light, Vector3D normal) {
    Color ambient = light.ambientIntensity * surface.color;

    // diffuse
    Vector3D lightDirection = (light.position - surfacePosition).normalized();
    f32 diff = fmax(normal.dot(lightDirection), 0.0);
    Color diffuse = diff * surface.color;

    // specular
    Vector3D viewDir = (viewPosition - surfacePosition).normalized();
    Vector3D reflectDir = (-lightDirection).reflectIn(normal);
    f32 spec = pow(fmax(viewDir.dot(reflectDir), 0.0), surface.shininess);
    Color specular = surface.specularIntensity * spec * surface.color;

    Color result = ambient + diffuse;// + specular;
    return result;
}

Color traceThroughScene(Ray ray, Scene scene, u32 traceDepth = 5) {
    SceneIntersectReport sceneIntersect = intersectScene(ray, scene);

    // Did not hit anything, so we can return the background color
    if(!sceneIntersect.hit.hit) {
        return {0, 0, 0, 1}; // todo: background color
    }

    Color surfaceColor = {0, 0, 0, 1};
    for(u32 i = 0; i < scene.numLights; i++) {
        Vector3D hitToLight = (scene.lights[i].position - sceneIntersect.hit.hitPosition);

        Ray shadowRay = {};
        shadowRay.start = sceneIntersect.hit.hitPosition;
        shadowRay.direction = hitToLight.normalized();

        SceneIntersectReport shadowTrace = intersectScene(shadowRay, scene);

        // improvised ambient term
        Color ambientTerm = 0.03 *sceneIntersect.hitMaterial.albedo * sceneIntersect.hitMaterial.ao;

        if(shadowTrace.hit.hit && shadowTrace.hit.TOI <= hitToLight.length()) {
            //surfaceColor = surfaceColor + (scene.lights[i].ambientIntensity * sceneIntersect.hitMaterial.color);
            surfaceColor = ambientTerm;
        } else {
            surfaceColor = surfaceColor + calculateSurfaceColorFromLightPBR(ray.start, sceneIntersect.hit.hitPosition, sceneIntersect.hitMaterial, scene.lights[i], sceneIntersect.hit.hitNormal);
            surfaceColor = surfaceColor + ambientTerm;
        }
    }

    //todo: reflection and refraction
    surfaceColor.a = 1; //todo: handle alpha
    return surfaceColor;
}


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

#define numTestSpheres 3
    SphereObject spheres[numTestSpheres] = {};
    spheres[0].sphere.position = {0, 2, 10};
    spheres[0].sphere.radius = 2;
    spheres[0].material = PBM_ROUGH_RED;
    //spheres[0].material.color = {(f32) 0xAA / 255, (f32) 0x00 / 255, (f32) 0x00 / 255, (f32) 0xFF / 255};
    //spheres[0].material.shininess = 0.25f;

    spheres[1].sphere.position = {8, 3, 9};
    spheres[1].sphere.radius = 3;
    spheres[1].material = PBM_METALLIC_GREEN;
    //spheres[1].material.color = {(f32) 0x00 / 255, (f32) 0xAA / 255, (f32) 0x00 / 255, (f32) 0xFF / 255};
    //spheres[1].material.shininess = 0.25f;

    spheres[2].sphere.position = {-2, 1, 8};
    spheres[2].sphere.radius = 1;
    spheres[2].material = PBM_SMOOTH_BLUE;
    //spheres[2].material.color = {(f32) 0x00 / 255, (f32) 0x00 / 255, (f32) 0xAA / 255, (f32) 0xFF / 255};
    //spheres[2].material.shininess = 0.25f;

#define numTestTriangles 4
    TriangleObject triangles[numTestTriangles];

    // floor

    triangles[0].triangle.B = {-50, 0, -50};
    triangles[0].triangle.A = { 50, 0, -50};
    triangles[0].triangle.C = {-50, 0,  50};
    triangles[0].material = PBM_GRAY;
    //triangles[0].material.color = {(f32) 0x20 / 255, (f32) 0x20 / 255, (f32) 0x20 / 255, (f32) 0xFF / 255};
    //triangles[0].material.shininess = 0.078125f;

    triangles[1].triangle.B = {-50, 0,  50};
    triangles[1].triangle.A = { 50, 0, -50};
    triangles[1].triangle.C = { 50, 0,  50};
    triangles[1].material = PBM_GRAY;
    //triangles[1].material.color = {(f32) 0x20 / 255, (f32) 0x20 / 255, (f32) 0x20 / 255, (f32) 0xFF / 255};
    //triangles[1].material.shininess = 0.078125f;

    // back wall

    triangles[2].triangle.B = {-50, 20,  50};
    triangles[2].triangle.A = {-50,  0,  50};
    triangles[2].triangle.C = { 50,  0,  50};
    triangles[2].material = PBM_GRAY;
    //triangles[2].material.color = {(f32) 0x20 / 255, (f32) 0x20 / 255, (f32) 0x20 / 255, (f32) 0xFF / 255};
    //triangles[2].material.shininess = 0.078125f;

    triangles[3].triangle.B = {-50, 20,  50};
    triangles[3].triangle.A = { 50,  0,  50};
    triangles[3].triangle.C = { 50, 20,  50};
    triangles[3].material = PBM_GRAY;
    //triangles[3].material.color = {(f32) 0x20 / 255, (f32) 0x20 / 255, (f32) 0x20 / 255, (f32) 0xFF / 255};
    //triangles[3].material.shininess = 0.078125f;

#define numTestLights 2
    PointLight lights[numTestLights] = {};
    lights[0].position = {-5,3, 5};
    lights[0].ambientIntensity = 0.01f;

    lights[1].position = {5,5, 4};
    lights[1].ambientIntensity = 0.05f;

    Scene scene = {};
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

    #if 0
    for(u32 y = 0; y < height; y++) {
        for(u32 x = 0; x < width; x++) {

            Vector3D wPixel = pixelToWorldSpace(camera, x, y, width, height);

            Vector3D rayP = wPixel;
            Vector3D rayD = (wPixel - camera.position).normalized();

            Ray ray = {};
            ray.start = rayP;
            ray.direction = rayD;

            SceneTraceReport traceReport = intersectScene(ray, scene);

            if(traceReport.hit.hit) {
                Vector3D hitToLight = (light.position - traceReport.hit.hitPosition);

                Ray shadowRay = {};
                shadowRay.start = traceReport.hit.hitPosition;
                shadowRay.direction = hitToLight.normalized();

                SceneTraceReport shadowTrace = intersectScene(shadowRay, scene);

                Color ambient = {0xFF202020};

                if(shadowTrace.hit.hit && shadowTrace.hit.TOI <= hitToLight.length()) {
                    pixels[x + y * width] = ambient;
                } else {
                    pixels[x + y * width] = ambient + fmax(0, shadowRay.direction.dot(traceReport.hit.hitNormal)) * traceReport.color;
                }

            } else {
                pixels[x + y * width] = {0xFF000000};
            }
        }
    }
    #endif

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


    saveOutBitmap(createBitmap(bmpPixels, width, height), "P:/raytracing/res/restracing_test.bmp");
    saveOutBitmap(createBitmap(bmpAaPixels, halfWidth, halfHeight), "P:/raytracing/res/aa_restracing_test.bmp");

    free(pixels);
    free(aaPixels);

    free(bmpPixels);
    free(bmpAaPixels);


    return 0;
}

