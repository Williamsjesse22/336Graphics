#ifndef QUAD_H
#define QUAD_H

#include "rtweekend.h"
#include "hittable.h"
#include "aabb.h"

class quad : public hittable
{
public:
    point3 Q; // one corner of the quad
    vec3 u;   // edge vector along "u"
    vec3 v;   // edge vector along "v"
    shared_ptr<material> mat;

    aabb bbox;
    vec3 normal;
    vec3 w; // helper for barycentric-style test
    double D;
    double area;

    quad() {}

    quad(const point3 &_Q, const vec3 &_u, const vec3 &_v, shared_ptr<material> m)
        : Q(_Q), u(_u), v(_v), mat(m)
    {
        vec3 n = cross(u, v);
        area = n.length();
        normal = unit_vector(n);
        w = n / dot(n, n);
        D = dot(normal, Q);

        // Build bounding box from the four corners
        point3 p0 = Q;
        point3 p1 = Q + u;
        point3 p2 = Q + v;
        point3 p3 = Q + u + v;

        auto xmin = fmin(fmin(p0.x(), p1.x()), fmin(p2.x(), p3.x()));
        auto ymin = fmin(fmin(p0.y(), p1.y()), fmin(p2.y(), p3.y()));
        auto zmin = fmin(fmin(p0.z(), p1.z()), fmin(p2.z(), p3.z()));

        auto xmax = fmax(fmax(p0.x(), p1.x()), fmax(p2.x(), p3.x()));
        auto ymax = fmax(fmax(p0.y(), p1.y()), fmax(p2.y(), p3.y()));
        auto zmax = fmax(fmax(p0.z(), p1.z()), fmax(p2.z(), p3.z()));

        bbox = aabb(point3(xmin, ymin, zmin), point3(xmax, ymax, zmax));
    }

    bool hit(const ray &r, interval ray_t, hit_record &rec) const override
    {
        // Intersect ray with the plane of the quad
        double denom = dot(normal, r.direction());
        if (fabs(denom) < 1e-8)
            return false;

        double t = (D - dot(normal, r.origin())) / denom;
        if (!ray_t.contains(t))
            return false;

        point3 p = r.at(t);

        // Express p in (u, v) coordinates:
        // planar_p = alpha * u + beta * v
        vec3 planar_p = p - Q;

        double alpha = dot(w, cross(planar_p, v));
        double beta = dot(w, cross(u, planar_p));

        if (alpha < 0.0 || alpha > 1.0 || beta < 0.0 || beta > 1.0)
            return false;

        rec.t = t;
        rec.p = p;
        rec.mat = mat;

        // Normal and face orientation
        vec3 outward_normal = normal;
        rec.set_face_normal(r, outward_normal);

        // Simple UV mapping across the quad
        rec.u = alpha;
        rec.v = beta;

        return true;
    }

    bool bounding_box(aabb &output_box) const override
    {
        output_box = bbox;
        return true;
    }

    // Optional: make quads usable as lights for importance sampling later
    double pdf_value(const point3 &origin, const vec3 &direction) const override
    {
        hit_record rec;
        if (!this->hit(ray(origin, direction), interval(0.001, infinity), rec))
            return 0;

        auto distance_squared = rec.t * rec.t * direction.length_squared();
        auto cosine = fabs(dot(direction, rec.normal)) / direction.length();
        return distance_squared / (cosine * area);
    }

    vec3 random(const point3 &origin) const override
    {
        auto r1 = random_double();
        auto r2 = random_double();
        point3 p = Q + r1 * u + r2 * v;
        return p - origin;
    }
};

#endif
