//
// Created by Tijn Bertens on 26-5-2019.
//

#ifndef RAYTRACING_TEXTURE_H
#define RAYTRACING_TEXTURE_H

#include "types.h"
#include "color.h"
#include "file_io.h"

/**
 * Raw pixel data that can be used for sampling.
 */
struct Texture {
    u32 width;              // Width of texture in pixels
    u32 height;             // Height of texture in pixels

    Color *pixelData;       // Pixel data
};

/**
 * Creates a texture from the contents of a bitmap.
 */
Texture textureFromBitmap(Bitmap bitmap) {
    Texture result = {};

    u32 width = bitmap.header.width;
    u32 height = bitmap.header.height;

    result.width = width;
    result.width = height;

    result.pixelData = (Color *) malloc(width * height);

    for(u32 x = 0; x < width; x++) {
        for(u32 y = 0; y < height; y++) {
            u32 i = y*width + x;
            result.pixelData[i] = toColor(bitmap.data[i]);
        }
    }

    return result;
}

/**
 * Samples color data from a texture, given normalized coordinates.
 * @param texture The texture to sample from.
 * @param u [0, 1] horizontal coordinate
 * @param v [0, 1] vertical coordinate
 */
Color sample(Texture texture, f32 u, f32 v) {
    u32 x = (u32) u * texture.width;
    u32 y = (u32) v * texture.height;

    return(texture.pixelData[y*texture.width + x]);
}

#endif //RAYTRACING_TEXTURE_H
