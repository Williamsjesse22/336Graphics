#ifndef HITTABLE_H
#define HITTABLE_H

#include "rtweekend.h"
#include "ray.h"
#include "interval.h"
#include "aabb.h"

class material;

struct hit_record
{
    point3 p;
    vec3 normal;
    shared_ptr<material> mat;
    double t;
    double u;
    double v;
    bool front_face;

    void set_face_normal(const ray &r, const vec3 &outward_normal)
    {
        front_face = dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal : -outward_normal;
    }
};

class hittable
{
public:
    virtual ~hittable() = default;

    virtual bool hit(const ray &r, interval ray_t, hit_record &rec) const = 0;

    // BVH needs this
    virtual bool bounding_box(aabb &output_box) const = 0;

    virtual double pdf_value(const point3 &o, const vec3 &v) const
    {
        return 0.0;
    }

    virtual vec3 random(const point3 &o) const
    {
        return vec3(1, 0, 0);
    }
};

#endif
