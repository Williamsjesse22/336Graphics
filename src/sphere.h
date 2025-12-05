#ifndef SPHERE_H
#define SPHERE_H

#include "rtweekend.h"
#include "hittable.h"
#include "aabb.h"
#include "onb.h"

class sphere : public hittable
{
public:
    point3 center;
    double radius;
    shared_ptr<material> mat;

    sphere() {}

    sphere(point3 cen, double r, shared_ptr<material> m)
        : center(cen), radius(r), mat(m) {}

    bool hit(const ray &r, interval ray_t, hit_record &rec) const override
    {
        vec3 oc = r.origin() - center;
        auto a = r.direction().length_squared();
        auto half_b = dot(oc, r.direction());
        auto c = oc.length_squared() - radius * radius;

        auto discriminant = half_b * half_b - a * c;
        if (discriminant < 0)
            return false;
        auto sqrtd = sqrt(discriminant);

        auto root = (-half_b - sqrtd) / a;
        if (!ray_t.surrounds(root))
        {
            root = (-half_b + sqrtd) / a;
            if (!ray_t.surrounds(root))
                return false;
        }

        rec.t = root;
        rec.p = r.at(rec.t);
        vec3 outward_normal = (rec.p - center) / radius;
        rec.set_face_normal(r, outward_normal);
        rec.mat = mat;

        // sphere UVs
        get_sphere_uv(outward_normal, rec.u, rec.v);

        return true;
    }

    bool bounding_box(aabb &output_box) const override
    {
        output_box = aabb(
            center - vec3(radius, radius, radius),
            center + vec3(radius, radius, radius));
        return true;
    }

    double pdf_value(const point3 &o, const vec3 &v) const override
    {
        hit_record rec;
        if (!this->hit(ray(o, v), interval(0.001, infinity), rec))
            return 0;

        auto cos_theta_max = sqrt(1 - radius * radius / (center - o).length_squared());
        auto solid_angle = 2 * pi * (1 - cos_theta_max);

        return 1 / solid_angle;
    }

    vec3 random(const point3 &o) const override
    {
        vec3 direction = center - o;
        auto distance_squared = direction.length_squared();
        onb uvw;
        uvw.build_from_w(direction);
        return uvw.local(random_to_sphere(radius, distance_squared));
    }

private:
    static void get_sphere_uv(const point3 &p, double &u, double &v)
    {
        auto theta = acos(-p.y());
        auto phi = atan2(-p.z(), p.x()) + pi;
        u = phi / (2 * pi);
        v = theta / pi;
    }
};

#endif
