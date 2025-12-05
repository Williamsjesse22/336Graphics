#ifndef MOVING_SPHERE_H
#define MOVING_SPHERE_H

#include "rtweekend.h"
#include "hittable.h"
#include "aabb.h"

class moving_sphere : public hittable
{
public:
    point3 center0, center1;
    double time0, time1;
    double radius;
    shared_ptr<material> mat;

    moving_sphere() {}

    moving_sphere(
        point3 _center0, point3 _center1,
        double _time0, double _time1,
        double _radius,
        shared_ptr<material> _mat)
        : center0(_center0), center1(_center1),
          time0(_time0), time1(_time1),
          radius(_radius), mat(_mat) {}

    point3 center(double time) const
    {
        return center0 + ((time - time0) / (time1 - time0)) * (center1 - center0);
    }

    virtual bool hit(const ray &r, interval ray_t, hit_record &rec) const override
    {
        vec3 oc = r.origin() - center(r.time());
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
        vec3 outward_normal = (rec.p - center(r.time())) / radius;
        rec.set_face_normal(r, outward_normal);
        rec.mat = mat;

        return true;
    }

    // IMPORTANT: matches your hittable base class exactly
    virtual bool bounding_box(aabb &output_box) const override
    {
        aabb box0(
            center0 - vec3(radius, radius, radius),
            center0 + vec3(radius, radius, radius));

        aabb box1(
            center1 - vec3(radius, radius, radius),
            center1 + vec3(radius, radius, radius));

        output_box = surrounding_box(box0, box1);
        return true;
    }
};

#endif
