#ifndef MATERIAL_H
#define MATERIAL_H

#include "rtweekend.h"
#include "ray.h"
#include "vec3.h"
#include "hittable.h"
#include "texture.h"

struct hit_record;

class material
{
public:
    virtual ~material() = default;

    virtual bool scatter(
        const ray &r_in, const hit_record &rec,
        color &attenuation, ray &scattered) const = 0;

    virtual color emitted(double u, double v, const point3 &p) const
    {
        return color(0, 0, 0);
    }

    // NEW: for importance sampling
    virtual double scattering_pdf(
        const ray &r_in, const hit_record &rec, const ray &scattered) const
    {
        return 0;
    }

    // NEW: lets ray_color know whether to use pdf math
    virtual bool is_specular() const { return false; }
};

// -------------------------
// Random helpers
// -------------------------
inline vec3 random_in_unit_sphere()
{
    while (true)
    {
        auto p = vec3::random(-1, 1);
        if (p.length_squared() >= 1)
            continue;
        return p;
    }
}

inline vec3 random_unit_vector()
{
    return unit_vector(random_in_unit_sphere());
}

inline vec3 reflect(const vec3 &v, const vec3 &n)
{
    return v - 2 * dot(v, n) * n;
}

inline vec3 refract(const vec3 &uv, const vec3 &n, double etai_over_etat)
{
    auto cos_theta = fmin(dot(-uv, n), 1.0);
    vec3 r_out_perp = etai_over_etat * (uv + cos_theta * n);
    vec3 r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.length_squared())) * n;
    return r_out_perp + r_out_parallel;
}

inline double reflectance(double cosine, double ref_idx)
{
    auto r0 = (1 - ref_idx) / (1 + ref_idx);
    r0 = r0 * r0;
    return r0 + (1 - r0) * pow((1 - cosine), 5);
}

// -------------------------
// Lambertian diffuse (TEXTURED)
// -------------------------
#include "onb.h" // add near top of file

class lambertian : public material
{
public:
    shared_ptr<texture> albedo;

    lambertian(const color &a)
        : albedo(make_shared<solid_color>(a)) {}

    lambertian(shared_ptr<texture> a)
        : albedo(a) {}

    bool scatter(const ray &r_in, const hit_record &rec,
                 color &attenuation, ray &scattered) const override
    {
        onb uvw;
        uvw.build_from_w(rec.normal);

        auto direction = uvw.local(random_cosine_direction());
        scattered = ray(rec.p, unit_vector(direction));

        attenuation = albedo->value(rec.u, rec.v, rec.p);
        return true;
    }

    double scattering_pdf(const ray &, const hit_record &rec,
                          const ray &scattered) const override
    {
        auto cosine = dot(rec.normal, unit_vector(scattered.direction()));
        return (cosine < 0) ? 0 : cosine / pi;
    }
};

// -------------------------
// Metal
// -------------------------
class metal : public material
{
public:
    color albedo;
    double fuzz;

    metal(const color &a, double f) : albedo(a), fuzz(f < 1 ? f : 1) {}

    bool is_specular() const override { return true; }

    bool scatter(const ray &r_in, const hit_record &rec,
                 color &attenuation, ray &scattered) const override
    {
        vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
        scattered = ray(rec.p, reflected + fuzz * random_in_unit_sphere());
        attenuation = albedo;
        return (dot(scattered.direction(), rec.normal) > 0);
    }
};

// -------------------------
// Dielectric (glass)
// -------------------------
class dielectric : public material
{
public:
    double ir;

    dielectric(double index_of_refraction) : ir(index_of_refraction) {}

    bool is_specular() const override { return true; }

    bool scatter(const ray &r_in, const hit_record &rec,
                 color &attenuation, ray &scattered) const override
    {
        attenuation = color(1.0, 1.0, 1.0);
        double refraction_ratio = rec.front_face ? (1.0 / ir) : ir;

        vec3 unit_direction = unit_vector(r_in.direction());
        double cos_theta = fmin(dot(-unit_direction, rec.normal), 1.0);
        double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

        bool cannot_refract = refraction_ratio * sin_theta > 1.0;
        vec3 direction;

        if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double())
            direction = reflect(unit_direction, rec.normal);
        else
            direction = refract(unit_direction, rec.normal, refraction_ratio);

        scattered = ray(rec.p, direction);
        return true;
    }
};

class diffuse_light : public material
{
public:
    shared_ptr<texture> emit;

    diffuse_light(shared_ptr<texture> a) : emit(a) {}
    diffuse_light(color c) : emit(make_shared<solid_color>(c)) {}

    bool scatter(const ray &, const hit_record &, color &, ray &) const override
    {
        return false; // light doesn't scatter
    }

    color emitted(double u, double v, const point3 &p) const override
    {
        return emit->value(u, v, p);
    }
};

#endif
