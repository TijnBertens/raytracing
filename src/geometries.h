//
// Created by s150818 on 14-2-2019.
//

#ifndef RAYTRACING_GEOMETRY_H
#define RAYTRACING_GEOMETRY_H

#include "math_utils.h"
#include "file_io.h"
#include "material.h"

/**
 * 3D triangle, defined by three points.
 */
struct Triangle {
    Vec3f A;
    Vec3f B;
    Vec3f C;
};

/**
 * 3D sphere, defined by center and radius.
 */
struct Sphere {
    Vec3f position;      // center
    f32 radius;
};

/**
 * 3D axix aligned bounding box.
 */
struct AABB {
    Vec3f min;
    Vec3f max;
};

/**
 * 3D mesh data. List of vertex coordinates and a list of indices that form triangles.
 * Also stores materials and per-face material indices.
 */
struct Mesh {
    u32 numTriangles;               // Number of triangles (faces) in the mesh

    u32 numVertices;                // Number of vertices
    Vec3f *vertices;                // Array of vertex positions

    u32 numIndices;                 // Number of vertex indices
    u32 *indices;                   // Array of vertex indices, creating triangles

    u32 *materialIndices;           // Per-triangle material indices

    u32 numMaterials;               // Number of materials in the material library
    PBMaterial *materials;          // Library of materials
};

/**
 * Load a mesh from a given file.
 */
bool loadMesh(Mesh *mesh, const char *file) {
    ObjContent objContent;
    if(!readInObj(&objContent, file)) {
        printf("Could not load mesh from file %s\n", file);
        return false;
    }

    // Copy vertex data

    mesh->numVertices = objContent.attrib.vertices.size() / 3;
    mesh->vertices = (Vec3f *) malloc(sizeof(Vec3f) * mesh->numVertices);
    // NOTE: kind of banking on the fact that real_t is a 32 bit float in the tiny obj library
    memcpy(mesh->vertices, objContent.attrib.vertices.data(), objContent.attrib.vertices.size() * sizeof(f32));

    // Count number of indices over all shapes

    mesh->numIndices = 0;
    for(const auto &shape : objContent.shapes) {
        mesh->numIndices += shape.mesh.indices.size();
    }
    mesh->numTriangles = mesh->numIndices / 3;

    // Copy vertex index data

    mesh->indices = (u32 *) malloc(sizeof(u32) * mesh->numIndices);
    u32 i = 0;
    for(const auto &shape : objContent.shapes) {
        for(const auto &index : shape.mesh.indices) {
            mesh->indices[i++] = index.vertex_index;
        }
    }

    // Copy material indices

    mesh->materialIndices = (u32 *) malloc(sizeof(u32) * mesh->numTriangles);
    u32 faceIndex = 0;
    for(const auto &shape : objContent.shapes) {
        for(const auto &matIndex : shape.mesh.material_ids) {
            mesh->materialIndices[faceIndex++] = matIndex;
        }
    }

    // Copy materials

    mesh->numMaterials = objContent.materials.size();
    mesh->materials = (PBMaterial *) malloc(sizeof(PBMaterial) * mesh->numMaterials);
    for(u32 i = 0; i < mesh->numMaterials; i++) {
        tinyobj::material_t &mat = objContent.materials[i];

        // Note: These are some gruesome hacks that get PBR material values from an MTL material.
        // We do this because Blender does not support the PBR extension for MTL, which tinyOBJ uses.
        // We could adapt the Blender exporter to support it.
        // https://developer.blender.org/diffusion/BA/browse/master/io_scene_obj/export_obj.py

        mesh->materials[i].albedo.r = mat.diffuse[0];
        mesh->materials[i].albedo.g = mat.diffuse[1];
        mesh->materials[i].albedo.b = mat.diffuse[2];
        mesh->materials[i].albedo.a = 1;

        mesh->materials[i].roughness = -(sqrtf(mat.shininess) / 30 - 1);

        if(mat.ambient[0] == 1.0 && mat.ambient[1] == 1.0 && mat.ambient[2] == 1.0) {
            mesh->materials[i].metallic = 0;
        } else {
            mesh->materials[i].metallic = mat.ambient[0];
        }

        mesh->materials[i].ao = 1;
    }

    return true;
}

#endif //RAYTRACING_GEOMETRY_H

