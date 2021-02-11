//
// Created by Tijn Bertens on 31-3-2019.
//

#ifndef RAYTRACING_FILE_IO_H
#define RAYTRACING_FILE_IO_H

#include <cstdio>
#include "math_utils.h"
#include "color.h"

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"
// Have to undefine because this is a single header library in a unity build ...
#undef TINYOBJLOADER_IMPLEMENTATION

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
 * Reads in a bitmap from a specified file path.
 */
bool readInBitmap(Bitmap *bitmap, const char *filePath) {
    FILE *inFile = fopen(filePath, "rb");

    if(inFile) {
        // Read bitmap header
        fread(&bitmap->header, 1, sizeof(BitmapHeader), inFile);

        // Allocate space for pixel data
        bitmap->data = (IntColor *) malloc(bitmap->header.dataSize);

        // Move to start of pixel data and read
        fseek(inFile, bitmap->header.dataOffset, SEEK_SET);
        fread(bitmap->data, 1, bitmap->header.dataSize, inFile);

        fclose(inFile);
        return true;
    } else {
        printf("Failed to read bitmap file, could not open file.");
        return false;
    }
}

/**
 * Contents of a .OBJ file.
 * This particular struct is basically a wrapper for the return values of loadObj from the tiny obj loader library.
 */
struct ObjContent {
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
};

/**
 * Reads in an object file from a specified file path.
 */
bool readInObj(ObjContent *objContent, const char *file) {
    std::string warn;
    std::string err;

    bool ret = tinyobj::LoadObj(&objContent->attrib, &objContent->shapes, &objContent->materials, &warn, &err, file, "../res/objects");

    if(warn.length()) {
        printf("%s\n", warn.c_str());
    }

    if(err.length()) {
        printf("%s\n", err.c_str());
    }

    return ret;
}


#endif //RAYTRACING_FILE_IO_H
