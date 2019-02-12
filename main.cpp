#include <malloc.h>
#include <cstdio>
#include "types.h"
#include "math_utils.h"

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
    u32 *data;
};

/**
 * Saves out a bitmap to a specified file path.
 */
bool saveOutBitmap(Bitmap bitmap, char *filePath) {
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
Bitmap createBitmap(u32 *data, u32 width, u32 height) {
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
    f32 screenWidth =  2 * screenHeight * camera.aspectRatio;

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

int main(){
    Camera camera = {};

    camera.position = {0, 0, 0};
    camera.viewDirection = {0, 0, 1};
    camera.upVector = {0, 1, 0};

    camera.fovy = 90;
    camera.nearClippingDistance = 1;
    camera.aspectRatio = 1280.0 / 720.0;


    Vector3D ws = pixelToWorldSpace(camera, 1280, 720, 1280, 720);

    printf("{%f, %f, %f} \n", ws.x, ws.y, ws.z);

    return 0;
}
