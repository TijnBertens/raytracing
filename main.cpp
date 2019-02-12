#include <malloc.h>
#include <cstdio>
#include "types.h"
#include "math_utils.h"

#pragma pack(push, 1)
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

struct Bitmap {
    BitmapHeader header;
    u32 *data;
};
#pragma pack(pop)

void saveOutBitmap(Bitmap bitmap) {
    FILE *outFile = fopen("P:/raytracing/res/test.bmp", "wb");

    if(outFile) {
        fwrite(&bitmap, sizeof(bitmap.header), 1, outFile);
        fwrite(bitmap.data, bitmap.header.dataSize, 1, outFile);
        fclose(outFile);
    } else {
        printf("Failed to save bitmap file, could not open file.");
    }
}

Bitmap createBitmap(u32 *data, u32 width, u32 height) {
    Bitmap result = {};

    result.header.fileType = 0x4D42; // 'BM'
    result.header.fileSize = sizeof(BitmapHeader) + (width * height * sizeof(u32));
    result.header.dataOffset = sizeof(BitmapHeader);
    result.header.headerSize = sizeof(BitmapHeader) - 14;
    result.header.width = width;
    result.header.height = height;
    result.header.numPlanes = 1;
    result.header.bitsPerPixel = 32;
    result.header.compression = 0;
    result.header.dataSize = (width * height * sizeof(u32));
    result.header.horzResolution = 5000;
    result.header.vertResolution = 5000;
    result.header.numColorsInPalette = 0;
    result.header.importantColors = 0;

    result.data = data;

    return result;
}

int main(){
    u32 width = 1280;
    u32 height = 720;

    u32 *data = (u32 *) malloc(sizeof(u32) * width * height);

    for(u32 x = 0; x < width; x++) {
        for(u32 y = 0; y < height; y++) {

            s32 xx = (s32) x - width/2;
            s32 yy = (s32) y - height/2;

            if(xx * xx + yy * yy < 200 * 200) {
                data[x + y * width] = 0xFFFFFFFF;
            } else {
                data[x + y * width] = 0xFF000000;
            }

        }
    }

    Bitmap bitmap = createBitmap(data, width, height);
    saveOutBitmap(bitmap);

    return 0;
}
