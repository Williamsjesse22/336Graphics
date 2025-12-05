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

void write_color(std::ostream &out, color pixel_color, int samples_per_pixel)
{
    auto r = pixel_color.x();
    auto g = pixel_color.y();
    auto b = pixel_color.z();

    // Average over samples
    auto scale = 1.0 / samples_per_pixel;
    r *= scale;
    g *= scale;
    b *= scale;

    // Gamma correction (gamma=2.0)
    r = linear_to_gamma(r);
    g = linear_to_gamma(g);
    b = linear_to_gamma(b);

    out << static_cast<int>(256 * clamp(r, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(g, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(b, 0.0, 0.999)) << '\n';
}

#endif
