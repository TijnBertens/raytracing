//
// Created by s150818 on 14-2-2019.
//

#ifndef RAYTRACING_GEOMETRY_H
#define RAYTRACING_GEOMETRY_H

#include "math_utils.h"
#include "file_io.h"

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
 */
struct Mesh {
    u32 numTriangles;

    u32 numVertices;
    Vec3f *vertices;

    u32 numIndices;
    u32 *indices;
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

    // Copy index data

    mesh->indices = (u32 *) malloc(sizeof(u32) * mesh->numIndices);
    u32 i = 0;
    for(const auto &shape : objContent.shapes) {
        for(const auto &index : shape.mesh.indices) {
            mesh->indices[i++] = index.vertex_index;
        }
    }

    return true;
}

#endif //RAYTRACING_GEOMETRY_H

