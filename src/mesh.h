#ifndef MESH_H
#define MESH_H

#include "rtweekend.h"
#include "triangle.h"
#include "tiny_obj_loader.h"
#include <string>
#include <vector>
#include <iostream>
#include <limits>

// Simple struct so we can transform mesh vertices
struct mesh_transform
{
    vec3 scale;
    vec3 translate;

    mesh_transform(vec3 s = vec3(1, 1, 1), vec3 t = vec3(0, 0, 0))
        : scale(s), translate(t) {}

    point3 apply(const point3 &p) const
    {
        return point3(
            p.x() * scale.x() + translate.x(),
            p.y() * scale.y() + translate.y(),
            p.z() * scale.z() + translate.z());
    }

    vec3 apply_normal(const vec3 &n) const
    {
        // For uniform scaling, normal just needs renormalization.
        return unit_vector(vec3(
            n.x() / scale.x(),
            n.y() / scale.y(),
            n.z() / scale.z()));
    }
};

inline std::vector<shared_ptr<hittable>> load_obj_mesh_as_triangles(
    const std::string &filename,
    shared_ptr<material> default_mat,
    const mesh_transform &xf = mesh_transform())
{
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string warn, err;

    bool ok = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, filename.c_str());

    if (!warn.empty())
        std::cerr << "tinyobj warning: " << warn << "\n";
    if (!err.empty())
        std::cerr << "tinyobj error: " << err << "\n";
    if (!ok)
    {
        std::cerr << "ERROR: Failed to load OBJ: " << filename << "\n";
        return {};
    }

    if (attrib.vertices.empty())
    {
        std::cerr << "ERROR: OBJ has no vertices: " << filename << "\n";
        return {};
    }

    // ----------------------------
    // AUTO-CENTER + AUTO-SCALE
    // ----------------------------
    double min_x = std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double min_z = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();
    double max_z = -std::numeric_limits<double>::infinity();

    for (size_t i = 0; i < attrib.vertices.size(); i += 3)
    {
        double x = attrib.vertices[i + 0];
        double y = attrib.vertices[i + 1];
        double z = attrib.vertices[i + 2];

        min_x = fmin(min_x, x);
        max_x = fmax(max_x, x);
        min_y = fmin(min_y, y);
        max_y = fmax(max_y, y);
        min_z = fmin(min_z, z);
        max_z = fmax(max_z, z);
    }

    point3 center_raw(
        0.5 * (min_x + max_x),
        0.5 * (min_y + max_y),
        0.5 * (min_z + max_z));

    double extent_x = max_x - min_x;
    double extent_y = max_y - min_y;
    double extent_z = max_z - min_z;
    double max_extent = fmax(extent_x, fmax(extent_y, extent_z));

    // scale mesh so its largest dimension becomes ~1.0
    double auto_scale = (max_extent > 0) ? (1.0 / max_extent) : 1.0;

    // ----------------------------
    // Helpers to read OBJ arrays
    // ----------------------------
    auto get_position = [&](int idx) -> point3
    {
        return point3(
            attrib.vertices[3 * idx + 0],
            attrib.vertices[3 * idx + 1],
            attrib.vertices[3 * idx + 2]);
    };

    auto get_normal = [&](int idx) -> vec3
    {
        return vec3(
            attrib.normals[3 * idx + 0],
            attrib.normals[3 * idx + 1],
            attrib.normals[3 * idx + 2]);
    };

    auto get_texcoord = [&](int idx) -> vec3
    {
        // texcoords are (u,v)
        return vec3(
            attrib.texcoords[2 * idx + 0],
            attrib.texcoords[2 * idx + 1],
            0.0);
    };

    std::vector<shared_ptr<hittable>> tris;

    for (const auto &shape : shapes)
    {
        size_t index_offset = 0;

        for (size_t f = 0; f < shape.mesh.num_face_vertices.size(); f++)
        {
            int fv = shape.mesh.num_face_vertices[f];

            if (fv < 3)
            { // degenerate face
                index_offset += fv;
                continue;
            }

            // Load all indices of this face
            std::vector<tinyobj::index_t> face_indices;
            face_indices.reserve(fv);
            for (int k = 0; k < fv; k++)
            {
                face_indices.push_back(shape.mesh.indices[index_offset + k]);
            }

            // Fan triangulate: (v0, v1, v2), (v0, v2, v3), ...
            for (int k = 1; k < fv - 1; k++)
            {
                tinyobj::index_t i0 = face_indices[0];
                tinyobj::index_t i1 = face_indices[k];
                tinyobj::index_t i2 = face_indices[k + 1];

                // raw positions
                point3 p0_raw = get_position(i0.vertex_index);
                point3 p1_raw = get_position(i1.vertex_index);
                point3 p2_raw = get_position(i2.vertex_index);

                // auto normalize: recenter then scale
                point3 p0 = (p0_raw - center_raw) * auto_scale;
                point3 p1 = (p1_raw - center_raw) * auto_scale;
                point3 p2 = (p2_raw - center_raw) * auto_scale;

                // apply user transform after normalization
                p0 = xf.apply(p0);
                p1 = xf.apply(p1);
                p2 = xf.apply(p2);

                bool has_normals = !attrib.normals.empty() && i0.normal_index >= 0 && i1.normal_index >= 0 && i2.normal_index >= 0;

                bool has_uvs = !attrib.texcoords.empty() && i0.texcoord_index >= 0 && i1.texcoord_index >= 0 && i2.texcoord_index >= 0;

                // TEMP FIX: ignore OBJ normals for teapot visibility.
                // This forces geometric (flat) normals via triangle.h
                has_normals = false;

                if (has_normals && has_uvs)
                {
                    vec3 n0 = xf.apply_normal(get_normal(i0.normal_index));
                    vec3 n1 = xf.apply_normal(get_normal(i1.normal_index));
                    vec3 n2 = xf.apply_normal(get_normal(i2.normal_index));

                    vec3 uv0 = get_texcoord(i0.texcoord_index);
                    vec3 uv1 = get_texcoord(i1.texcoord_index);
                    vec3 uv2 = get_texcoord(i2.texcoord_index);

                    tris.push_back(make_shared<triangle>(
                        p0, p1, p2,
                        n0, n1, n2,
                        uv0, uv1, uv2,
                        default_mat));
                }
                else
                {
                    tris.push_back(make_shared<triangle>(p0, p1, p2, default_mat));
                }
            }

            index_offset += fv;
        }
    }

    std::cerr << "Loaded OBJ '" << filename << "' -> " << tris.size() << " triangles.\n";
    return tris;
}

#endif
