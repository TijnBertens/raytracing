//
// Created by Tijn Bertens on 31-3-2019.
//

#ifndef RAYTRACING_FILE_IO_H
#define RAYTRACING_FILE_IO_H

#include <cstdio>
#include "math_utils.h"
#include "color.h"

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
 * Saves out a bitmap to a specified file path.
 */
bool readInBitmap(Bitmap *output, const char *filePath) {
    FILE *inputFile = fopen(filePath, "rb");

    if(!inputFile) {
        printf("Failed to load bitmap file, could not open file.");
        return false;
    }

    fread
}


#endif //RAYTRACING_FILE_IO_H
