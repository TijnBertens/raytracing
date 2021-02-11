//
// Created by Tijn Bertens on 31-7-2019.
//

#ifndef RAYTRACING_MATERIAL_H
#define RAYTRACING_MATERIAL_H

//  ----------------------
// Physically based material
//  ----------------------

/**
 * Physically based material.
 */
struct PBMaterial {
    Color albedo;
    f32 metallic;
    f32 roughness;
    f32 ao;
};

static const PBMaterial PBM_ROUGH_RED = {
        {1, 0, 0, 1},
        0,
        1,
        1
};

static const PBMaterial PBM_METALLIC_GREEN = {
        {0, 1, 0, 1},
        1,
        0.2,
        1
};

static const PBMaterial PBM_SMOOTH_BLUE = {
        {0, 0, 1, 1},
        0,
        0.1,
        1
};

static const PBMaterial PBM_GRAY = {
        {0.8, 0.8, 0.8, 8},
        0,
        1,
        1
};

static const PBMaterial PBM_SMOOTH_OFF_WHITE = {
        {0.9, 0.9, 0.9, 1},
        0,
        0.2,
        1
};

static const PBMaterial PBM_ROUGH_OFF_WHITE = {
        {0.95, 0.95, 0.95, 1},
        0,
        0.95,
        1
};


#endif //RAYTRACING_MATERIAL_H
