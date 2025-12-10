#ifndef TEXTURE_H
#define TEXTURE_H

#include "rtweekend.h"
#include "vec3.h"
#include "perlin.h"
#include <memory>

class texture
{
public:
    virtual ~texture() = default;
    virtual color value(double u, double v, const point3 &p) const = 0;
};

class solid_color : public texture
{
public:
    color albedo;

    solid_color() {}
    solid_color(const color &c) : albedo(c) {}
    solid_color(double r, double g, double b) : albedo(r, g, b) {}

    color value(double u, double v, const point3 &p) const override
    {
        return albedo;
    }
};

class checker_texture : public texture
{
public:
    std::shared_ptr<texture> even;
    std::shared_ptr<texture> odd;

    checker_texture() {}

    checker_texture(std::shared_ptr<texture> _even, std::shared_ptr<texture> _odd)
        : even(_even), odd(_odd) {}

    checker_texture(const color &c1, const color &c2)
        : even(std::make_shared<solid_color>(c1)),
          odd(std::make_shared<solid_color>(c2)) {}

    color value(double u, double v, const point3 &p) const override
    {
        auto sines = sin(10 * p.x()) * sin(10 * p.y()) * sin(10 * p.z());
        if (sines < 0)
            return odd->value(u, v, p);
        else
            return even->value(u, v, p);
    }
};

// ------------------------------------------------------------
// Perlin noise texture (for clouds / foggy materials)
// ------------------------------------------------------------
class noise_texture : public texture
{
public:
    perlin noise;
    double scale;

    noise_texture() : scale(1.0) {}
    noise_texture(double sc) : scale(sc) {}

    color value(double u, double v, const point3 &p) const override
    {
        // Classic RTW-style "marble / cloud" turbulence
        auto t = 0.5 * (1.0 + sin(scale * p.z() + 10.0 * noise.turb(p)));
        return color(1.0, 1.0, 1.0) * t; // grayscale foggy look
    }
};

#endif
