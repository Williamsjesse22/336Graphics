#ifndef BRICK_TEXTURE_H
#define BRICK_TEXTURE_H

#include "texture.h"
#include <cmath>

// Procedural brick texture:
// - brick_color: fill color of the bricks
// - mortar_color: color of the joints between bricks
// - bricks_u: how many bricks across (horizontal)
// - bricks_v: how many brick rows (vertical)
// - mortar_frac: thickness of mortar as a fraction of brick size
class brick_texture : public texture
{
public:
    color brick_color;
    color mortar_color;
    int bricks_u;
    int bricks_v;
    double mortar_frac; // [0, 0.5) reasonable

    brick_texture(
        const color &brick,
        const color &mortar,
        int bricks_u_ = 12, // more columns ⇒ thinner bricks
        int bricks_v_ = 6,  // fewer rows ⇒ taller bricks
        double mortar_frac_ = 0.10)
        : brick_color(brick),
          mortar_color(mortar),
          bricks_u(bricks_u_),
          bricks_v(bricks_v_),
          mortar_frac(mortar_frac_)
    {
    }

    virtual color value(double u, double v, const point3 &p) const override
    {
        // Make sure UVs are positive so repetition works cleanly
        u = std::fabs(u);
        v = std::fabs(v);

        // Scale UVs into "brick space"
        double u_scaled = u * bricks_u;
        double v_scaled = v * bricks_v;

        // Which row are we in?
        int row = static_cast<int>(std::floor(v_scaled));

        // Stagger every other row by half a brick horizontally
        double u_offset = u_scaled + 0.5 * (row & 1); // row&1 is 0 or 1

        // Local coordinates inside the current brick cell
        double fu = u_offset - std::floor(u_offset); // [0,1) in brick X
        double fv = v_scaled - row;                  // [0,1) in brick Y

        // Determine if we are in the mortar region.
        // mortar_frac is the thickness of mortar on each side of the brick.
        if (fu < mortar_frac || fu > 1.0 - mortar_frac ||
            fv < mortar_frac || fv > 1.0 - mortar_frac)
        {
            return mortar_color;
        }

        return brick_color;
    }
};

#endif // BRICK_TEXTURE_H
