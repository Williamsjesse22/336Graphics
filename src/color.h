#ifndef COLOR_H
#define COLOR_H

#include "vec3.h"
#include "rtweekend.h"
#include <iostream>

inline double linear_to_gamma(double linear_component)
{
    if (linear_component <= 0)
        return 0;
    return sqrt(linear_component);
}

inline void write_color(std::ostream &out, color pixel_color, int samples_per_pixel)
{
    auto r = pixel_color.x();
    auto g = pixel_color.y();
    auto b = pixel_color.z();

    // Divide the color by the number of samples and gamma-correct for gamma=2.0.
    auto scale = 1.0 / samples_per_pixel;
    r = sqrt(scale * r);
    g = sqrt(scale * g);
    b = sqrt(scale * b);

    // Clamp and convert to [0,255]
    unsigned char pixel[3];
    pixel[0] = static_cast<unsigned char>(256 * clamp(r, 0.0, 0.999));
    pixel[1] = static_cast<unsigned char>(256 * clamp(g, 0.0, 0.999));
    pixel[2] = static_cast<unsigned char>(256 * clamp(b, 0.0, 0.999));

    out.write(reinterpret_cast<char *>(pixel), 3);
}

inline void write_color_hdr(std::ostream &out, color pixel_color, int samples_per_pixel)
{
    // Average the samples
    auto scale = 1.0 / samples_per_pixel;
    double r = pixel_color.x() * scale;
    double g = pixel_color.y() * scale;
    double b = pixel_color.z() * scale;

    // Clamp negatives to 0, but DO NOT clamp the upper bound.
    if (r < 0.0)
        r = 0.0;
    if (g < 0.0)
        g = 0.0;
    if (b < 0.0)
        b = 0.0;

    // Store as 32-bit floats, linear space (no gamma).
    float fr = static_cast<float>(r);
    float fg = static_cast<float>(g);
    float fb = static_cast<float>(b);

    out.write(reinterpret_cast<char *>(&fr), sizeof(float));
    out.write(reinterpret_cast<char *>(&fg), sizeof(float));
    out.write(reinterpret_cast<char *>(&fb), sizeof(float));
}

#endif
