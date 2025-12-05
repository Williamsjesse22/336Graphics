#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "rtweekend.h"
#include "hittable.h"
#include "aabb.h"

class triangle : public hittable
{
public:
    point3 v0, v1, v2;
    vec3 n0, n1, n2;
    vec3 uv0, uv1, uv2;
    bool has_vertex_normals = false;
    bool has_uvs = false;
    shared_ptr<material> mat;

    triangle(const point3 &_v0, const point3 &_v1, const point3 &_v2, shared_ptr<material> m)
        : v0(_v0), v1(_v1), v2(_v2), mat(m) {}

    triangle(const point3 &_v0, const point3 &_v1, const point3 &_v2,
             const vec3 &_n0, const vec3 &_n1, const vec3 &_n2,
             const vec3 &_uv0, const vec3 &_uv1, const vec3 &_uv2,
             shared_ptr<material> m)
        : v0(_v0), v1(_v1), v2(_v2),
          n0(_n0), n1(_n1), n2(_n2),
          uv0(_uv0), uv1(_uv1), uv2(_uv2),
          has_vertex_normals(true), has_uvs(true),
          mat(m) {}

    bool hit(const ray &r, interval ray_t, hit_record &rec) const override
    {
        vec3 edge1 = v1 - v0;
        vec3 edge2 = v2 - v0;

        vec3 h = cross(r.direction(), edge2);
        double a = dot(edge1, h);

        if (fabs(a) < 1e-8)
            return false;
        double f = 1.0 / a;

        vec3 s = r.origin() - v0;
        double u = f * dot(s, h);
        if (u < 0.0 || u > 1.0)
            return false;

        vec3 q = cross(s, edge1);
        double v = f * dot(r.direction(), q);
        if (v < 0.0 || u + v > 1.0)
            return false;

        double t = f * dot(edge2, q);
        if (!ray_t.surrounds(t))
            return false;

        rec.t = t;
        rec.p = r.at(t);
        rec.mat = mat;

        vec3 outward_normal;

        if (has_vertex_normals)
        {
            // barycentric normal interpolation
            outward_normal = unit_vector((1 - u - v) * n0 + u * n1 + v * n2);
        }
        else
        {
            outward_normal = unit_vector(cross(edge1, edge2));
        }

        rec.set_face_normal(r, outward_normal);

        if (has_uvs)
        {
            vec3 uv = (1 - u - v) * uv0 + u * uv1 + v * uv2;
            rec.u = uv.x();
            rec.v = uv.y();
        }

        return true;
    }

    bool bounding_box(aabb &output_box) const override
    {
        auto min_x = fmin(v0.x(), fmin(v1.x(), v2.x()));
        auto min_y = fmin(v0.y(), fmin(v1.y(), v2.y()));
        auto min_z = fmin(v0.z(), fmin(v1.z(), v2.z()));

        auto max_x = fmax(v0.x(), fmax(v1.x(), v2.x()));
        auto max_y = fmax(v0.y(), fmax(v1.y(), v2.y()));
        auto max_z = fmax(v0.z(), fmax(v1.z(), v2.z()));

        // pad a bit so thin triangles still have volume
        const double eps = 1e-4;
        output_box = aabb(
            point3(min_x - eps, min_y - eps, min_z - eps),
            point3(max_x + eps, max_y + eps, max_z + eps));
        return true;
    }
};

#endif
