#ifndef IMAGE_TEXTURE_H
#define IMAGE_TEXTURE_H

#include "texture.h"
#include "stb_image.h"
#include <string>

class image_texture : public texture
{
public:
    image_texture() : data(nullptr), width(0), height(0), bytes_per_scanline(0) {}

    image_texture(const std::string &filename)
    {
        auto components_per_pixel = bytes_per_pixel;

        data = stbi_load(filename.c_str(), &width, &height, &components_per_pixel, components_per_pixel);

        if (!data)
        {
            width = height = 0;
            bytes_per_scanline = 0;
            std::cerr << "ERROR: Could not load texture image file '" << filename << "'.\n";
        }
        else
        {
            bytes_per_scanline = bytes_per_pixel * width;
        }
    }

    ~image_texture() override
    {
        if (data)
            stbi_image_free(data);
    }

    color value(double u, double v, const point3 &p) const override
    {
        if (!data)
            return color(0, 1, 1); // cyan debug fallback

        // Clamp input texture coords to [0,1]
        u = clamp(u, 0.0, 1.0);
        v = 1.0 - clamp(v, 0.0, 1.0); // flip V to image coords

        auto i = static_cast<int>(u * width);
        auto j = static_cast<int>(v * height);

        if (i >= width)
            i = width - 1;
        if (j >= height)
            j = height - 1;

        const auto color_scale = 1.0 / 255.0;
        auto pixel = data + j * bytes_per_scanline + i * bytes_per_pixel;

        return color(color_scale * pixel[0],
                     color_scale * pixel[1],
                     color_scale * pixel[2]);
    }

private:
    unsigned char *data;
    int width, height;
    int bytes_per_scanline;
    static const int bytes_per_pixel = 3;
};

#endif
