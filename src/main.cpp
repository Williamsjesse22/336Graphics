#include "rtweekend.h"

#include "color.h"
#include "camera.h"
#include "hittable_list.h"
#include "sphere.h"
#include "moving_sphere.h"
#include "triangle.h"
#include "material.h"
#include "texture.h"
#include "image_texture.h"
#include "mesh.h"
#include "bvh_node.h"
#include "quad.h"
#include "brick_texture.h"
#include <thread>
#include <atomic>
#include <vector>

#include <iostream>
#include <chrono>
#include <fstream>

// Axis-aligned box made from 6 quads.
// p0 = minimum corner (x0,y0,z0), p1 = maximum corner (x1,y1,z1).
void add_box(hittable_list &world, const point3 &p0, const point3 &p1,
             shared_ptr<material> mat)
{
    const auto x0 = p0.x();
    const auto y0 = p0.y();
    const auto z0 = p0.z();
    const auto x1 = p1.x();
    const auto y1 = p1.y();
    const auto z1 = p1.z();

    // FRONT (+z)
    world.add(make_shared<quad>(
        point3(x0, y0, z1),
        vec3(x1 - x0, 0, 0),
        vec3(0, y1 - y0, 0),
        mat));

    // BACK (−z)
    world.add(make_shared<quad>(
        point3(x1, y0, z0),
        vec3(x0 - x1, 0, 0),
        vec3(0, y1 - y0, 0),
        mat));

    // LEFT (−x)
    world.add(make_shared<quad>(
        point3(x0, y0, z0),
        vec3(0, 0, z1 - z0),
        vec3(0, y1 - y0, 0),
        mat));

    // RIGHT (+x)
    world.add(make_shared<quad>(
        point3(x1, y0, z1),
        vec3(0, 0, z0 - z1),
        vec3(0, y1 - y0, 0),
        mat));

    // TOP (+y)
    world.add(make_shared<quad>(
        point3(x0, y1, z1),
        vec3(x1 - x0, 0, 0),
        vec3(0, 0, z0 - z1),
        mat));

    // BOTTOM (−y)
    world.add(make_shared<quad>(
        point3(x0, y0, z0),
        vec3(x1 - x0, 0, 0),
        vec3(0, 0, z1 - z0),
        mat));
}

// -----------------------------------------------------------------------------
// Instanced chair builder
// -----------------------------------------------------------------------------

enum class chair_facing
{
    PosX, // sitter faces +x
    NegX, // sitter faces -x
    PosZ, // sitter faces +z
    NegZ  // sitter faces -z
};

// Build a single chair given its seat footprint and a facing direction.
// The Y-dimensions, leg sizes, and backrest height match the tuned chair
// you already liked.
void add_chair(hittable_list &world,
               double seat_x_min, double seat_x_max,
               double seat_z_min, double seat_z_max,
               chair_facing facing,
               shared_ptr<material> material_chair)
{
    // -----------------------------
    // Shared dimensions (same as before)
    // -----------------------------
    const double seat_y_bottom = 0.25;
    const double seat_y_top = 0.35;

    // Seat box
    add_box(world,
            point3(seat_x_min, seat_y_bottom, seat_z_min),
            point3(seat_x_max, seat_y_top, seat_z_max),
            material_chair);

    // -----------------------------
    // Backrest: which side depends on facing
    // -----------------------------
    const double backrest_thickness = 0.10;

    switch (facing)
    {
    case chair_facing::PosX:
    {
        // Sitter faces +x, so backrest is on the low-x side
        double back_x_min = seat_x_min - backrest_thickness;
        double back_x_max = seat_x_min;
        add_box(world,
                point3(back_x_min, seat_y_top, seat_z_min),
                point3(back_x_max, 1.30, seat_z_max),
                material_chair);
        break;
    }
    case chair_facing::NegX:
    {
        // Sitter faces -x, backrest on high-x side
        double back_x_min = seat_x_max;
        double back_x_max = seat_x_max + backrest_thickness;
        add_box(world,
                point3(back_x_min, seat_y_top, seat_z_min),
                point3(back_x_max, 1.30, seat_z_max),
                material_chair);
        break;
    }
    case chair_facing::PosZ:
    {
        // Sitter faces +z, backrest on low-z side
        double back_z_min = seat_z_min - backrest_thickness;
        double back_z_max = seat_z_min;
        add_box(world,
                point3(seat_x_min, seat_y_top, back_z_min),
                point3(seat_x_max, 1.30, back_z_max),
                material_chair);
        break;
    }
    case chair_facing::NegZ:
    {
        // Sitter faces -z, backrest on high-z side
        double back_z_min = seat_z_max;
        double back_z_max = seat_z_max + backrest_thickness;
        add_box(world,
                point3(seat_x_min, seat_y_top, back_z_min),
                point3(seat_x_max, 1.30, back_z_max),
                material_chair);
        break;
    }
    }

    // -----------------------------
    // Legs (same logic for all facings)
    // -----------------------------
    const double leg_y_bottom = -0.50;
    const double leg_y_top = seat_y_bottom;

    const double leg_width = 0.12;
    const double leg_depth = 0.12;
    const double inset_x = 0.10;
    const double inset_z = 0.00; // you liked these values as-is

    // X positions
    const double left_leg_x_min = seat_x_min + inset_x;
    const double left_leg_x_max = left_leg_x_min + leg_width;
    const double right_leg_x_max = seat_x_max - inset_x;
    const double right_leg_x_min = right_leg_x_max - leg_width;

    // Z positions
    const double front_leg_z_max = seat_z_max - inset_z;
    const double front_leg_z_min = front_leg_z_max - leg_depth;
    const double back_leg_z_min = seat_z_min + inset_z;
    const double back_leg_z_max = back_leg_z_min + leg_depth;

    // Front-left
    add_box(world,
            point3(left_leg_x_min, leg_y_bottom, front_leg_z_min),
            point3(left_leg_x_max, leg_y_top, front_leg_z_max),
            material_chair);

    // Front-right
    add_box(world,
            point3(right_leg_x_min, leg_y_bottom, front_leg_z_min),
            point3(right_leg_x_max, leg_y_top, front_leg_z_max),
            material_chair);

    // Back-left
    add_box(world,
            point3(left_leg_x_min, leg_y_bottom, back_leg_z_min),
            point3(left_leg_x_max, leg_y_top, back_leg_z_max),
            material_chair);

    // Back-right
    add_box(world,
            point3(right_leg_x_min, leg_y_bottom, back_leg_z_min),
            point3(right_leg_x_max, leg_y_top, back_leg_z_max),
            material_chair);
}

color ray_color(const ray &r, const hittable &world, const hittable &lights, int depth)
{
    if (depth <= 0)
        return color(0, 0, 0);

    hit_record rec;
    if (!world.hit(r, interval(0.001, infinity), rec))
        return color(0, 0, 0);

    color emitted = rec.mat->emitted(rec.u, rec.v, rec.p);

    ray scattered;
    color attenuation;

    if (!rec.mat->scatter(r, rec, attenuation, scattered))
        return emitted;

    // Specular: one continuation ray, no PDFs
    if (rec.mat->is_specular())
        return emitted + attenuation * ray_color(scattered, world, lights, depth - 1);

    // --- MIS for diffuse (ONE recursive ray) ---
    ray new_ray;
    double pdf;

    if (random_double() < 0.5)
    {
        // Strategy A: sample a direction toward the lights
        vec3 light_dir = lights.random(rec.p);
        new_ray = ray(rec.p, unit_vector(light_dir), r.time());

        double light_pdf = lights.pdf_value(rec.p, new_ray.direction());
        double scatter_pdf = rec.mat->scattering_pdf(r, rec, new_ray);

        pdf = 0.5 * light_pdf + 0.5 * scatter_pdf;
    }
    else
    {
        // Strategy B: use the cosine-weighted scattered ray you already generated
        new_ray = scattered;

        double light_pdf = lights.pdf_value(rec.p, new_ray.direction());
        double scatter_pdf = rec.mat->scattering_pdf(r, rec, new_ray);

        pdf = 0.5 * light_pdf + 0.5 * scatter_pdf;
    }

    if (pdf <= 1e-8)
        pdf = 1e-8;

    double scattering_pdf = rec.mat->scattering_pdf(r, rec, new_ray);

    return emitted + attenuation * scattering_pdf *
                         ray_color(new_ray, world, lights, depth - 1) / pdf;
}

int main()
{
    // ------------------------------------------------------------
    // IMAGE SETTINGS
    // ------------------------------------------------------------
    const auto aspect_ratio = 16.0 / 9.0;

    const int image_width = 1280; // <- bump to 4000 for final
    const int image_height = static_cast<int>(image_width / aspect_ratio);

    // Adaptive sampling settings
    const int max_samples_per_pixel = 200; // hard cap
    const int min_samples_per_pixel = 20;  // minimum before we even check convergence
    const int sample_batch_size = 10;      // check every N samples
    const double target_noise = 0.02;

    const int max_depth = 20;

    // ------------------------------------------------------------
    // TEXTURES
    // ------------------------------------------------------------
    auto checker_tex = make_shared<checker_texture>(
        color(0.4, 0.4, 0.4),
        color(0.98, 0.98, 0.98));

    auto earth_tex = make_shared<image_texture>("earthmap.jpg");

    auto brick_tex = make_shared<brick_texture>(
        color(0.60, 0.10, 0.08),
        color(0.96, 0.93, 0.85),
        12,
        6,
        0.10);

    // ------------------------------------------------------------
    // MATERIALS
    // ------------------------------------------------------------
    auto material_ground = make_shared<lambertian>(checker_tex);
    auto material_teapot = make_shared<lambertian>(color(0.9, 0.2, 0.2));
    auto material_earth = make_shared<lambertian>(earth_tex);
    auto material_fore = make_shared<lambertian>(color(0.2, 0.8, 0.9));
    auto material_back = make_shared<lambertian>(color(0.9, 0.9, 0.2));
    auto material_moving = make_shared<lambertian>(color(0.7, 0.2, 0.9));
    auto light = make_shared<diffuse_light>(color(7, 7, 7));
    auto material_wall = make_shared<lambertian>(brick_tex);
    auto material_chair = make_shared<lambertian>(color(0.45, 0.25, 0.08));

    auto perlin_tex = make_shared<noise_texture>(2.0);
    auto material_fog = make_shared<lambertian>(perlin_tex);

    auto material_table = make_shared<lambertian>(color(0.55, 0.27, 0.07));

    // ------------------------------------------------------------
    // WORLD CONTENT
    // ------------------------------------------------------------
    hittable_list world;

    // Floor
    world.add(make_shared<sphere>(
        point3(0.0, -100.5, -1.0),
        100.0,
        material_ground));

    // Back wall
    {
        point3 wall_origin(-6.0, -0.5, -8.0);
        vec3 wall_u(12.0, 0.0, 0.0);
        vec3 wall_v(0.0, 4.5, 0.0);

        world.add(make_shared<quad>(wall_origin, wall_u, wall_v, material_wall));
    }

    // Left wall
    {
        point3 left_origin(-6.0, -0.5, 6.0);
        vec3 left_u(0.0, 0.0, -14.0);
        vec3 left_v(0.0, 4.5, 0.0);

        world.add(make_shared<quad>(left_origin, left_u, left_v, material_wall));
    }

    // Tabletop
    {
        point3 table_origin(-1.5, 0.5, -3.0);
        vec3 table_u(3.0, 0.0, 0.0);
        vec3 table_v(0.0, 0.0, 2.4);

        world.add(make_shared<quad>(table_origin, table_u, table_v, material_table));
    }

    // Table pedestal
    {
        const double pedestal_bottom_y = -0.5;
        const double pedestal_top_y = 0.48;
        const double pedestal_half_x = 0.25;
        const double pedestal_half_z = 0.5;

        const double z_front = -1.3;
        const double z_back = -2.3;

        // FRONT
        {
            point3 origin(-pedestal_half_x, pedestal_bottom_y, z_front);
            vec3 u(2.0 * pedestal_half_x, 0.0, 0.0);
            vec3 v(0.0, pedestal_top_y - pedestal_bottom_y, 0.0);
            world.add(make_shared<quad>(origin, u, v, material_table));
        }
        // BACK
        {
            point3 origin(-pedestal_half_x, pedestal_bottom_y, z_back);
            vec3 u(2.0 * pedestal_half_x, 0.0, 0.0);
            vec3 v(0.0, pedestal_top_y - pedestal_bottom_y, 0.0);
            world.add(make_shared<quad>(origin, u, v, material_table));
        }
        // LEFT
        {
            point3 origin(-pedestal_half_x, pedestal_bottom_y, z_front);
            vec3 u(0.0, 0.0, z_back - z_front);
            vec3 v(0.0, pedestal_top_y - pedestal_bottom_y, 0.0);
            world.add(make_shared<quad>(origin, u, v, material_table));
        }
        // RIGHT
        {
            point3 origin(pedestal_half_x, pedestal_bottom_y, z_front);
            vec3 u(0.0, 0.0, z_back - z_front);
            vec3 v(0.0, pedestal_top_y - pedestal_bottom_y, 0.0);
            world.add(make_shared<quad>(origin, u, v, material_table));
        }
    }

    // -----------------------------------------------------------------
    // CHAIRS: four instances around the table using add_chair()
    // -----------------------------------------------------------------

    // Chair 1: left of table, facing +x (what you had tuned)
    {
        const double seat_width_x = 1.0;
        const double seat_depth_z = 0.8;
        const double seat_x_max = -1.7;
        const double seat_x_min = seat_x_max - seat_width_x;
        const double seat_z_center = -1.8;
        const double seat_z_min = seat_z_center - 0.5 * seat_depth_z;
        const double seat_z_max = seat_z_center + 0.5 * seat_depth_z;

        add_chair(world,
                  seat_x_min, seat_x_max,
                  seat_z_min, seat_z_max,
                  chair_facing::PosX,
                  material_chair);
    }

    // Chair 2: right of table, facing -x
    {
        const double seat_width_x = 1.0;
        const double seat_depth_z = 0.8;
        const double seat_x_min = 1.7;
        const double seat_x_max = seat_x_min + seat_width_x;
        const double seat_z_center = -1.8;
        const double seat_z_min = seat_z_center - 0.5 * seat_depth_z;
        const double seat_z_max = seat_z_center + 0.5 * seat_depth_z;

        add_chair(world,
                  seat_x_min, seat_x_max,
                  seat_z_min, seat_z_max,
                  chair_facing::NegX,
                  material_chair);
    }

    // Chair 3: in front of table, facing -z (toward table)
    {
        const double seat_width_z = 1.2;
        const double seat_depth_x = 0.8;
        const double seat_x_center = 0.0;
        const double seat_z_max = -0.35;
        const double seat_z_min = seat_z_max - seat_width_z;
        const double seat_x_min = seat_x_center - 0.5 * seat_depth_x;
        const double seat_x_max = seat_x_center + 0.5 * seat_depth_x;

        add_chair(world,
                  seat_x_min, seat_x_max,
                  seat_z_min, seat_z_max,
                  chair_facing::NegZ,
                  material_chair);
    }

    // Chair 4: behind table, facing +z
    {
        const double seat_width_z = 1.2;
        const double seat_depth_x = 0.8;
        const double seat_x_center = 0.0;
        const double seat_z_min = -3.25;
        const double seat_z_max = seat_z_min + seat_width_z;
        const double seat_x_min = seat_x_center - 0.5 * seat_depth_x;
        const double seat_x_max = seat_x_center + 0.5 * seat_depth_x;

        add_chair(world,
                  seat_x_min, seat_x_max,
                  seat_z_min, seat_z_max,
                  chair_facing::PosZ,
                  material_chair);
    }

    // Perlin fog sphere
    world.add(make_shared<sphere>(
        point3(0.0, 1.2, -6.0),
        1.6,
        material_fog));

    // Earth sphere
    world.add(make_shared<sphere>(
        point3(2.2, 1.0, -3.5),
        1.0,
        material_earth));

    // Teapot
    mesh_transform teapot_xf(
        vec3(1.4, 1.4, 1.4),
        vec3(0.0, 0.9, -1.8));

    auto teapot_tris = load_obj_mesh_as_triangles(
        "teapot.obj",
        material_teapot,
        teapot_xf);

    for (auto &t : teapot_tris)
        world.add(t);

    // Overhead light sphere
    auto light_center = point3(0.0, 3.0, -1.8);
    auto light_radius = 1.5;
    world.add(make_shared<sphere>(
        light_center,
        light_radius,
        light));

    // Foreground / background spheres
    world.add(make_shared<sphere>(
        point3(-1.0, 0.6, 1.5),
        0.6,
        material_fore));

    world.add(make_shared<sphere>(
        point3(2.5, 1.0, -8.0),
        1.2,
        material_back));

    // Moving sphere (motion blur)
    double shutter_open = 0.0;
    double shutter_close = 1.0;

    world.add(make_shared<moving_sphere>(
        point3(-1.2, 0.5, -1.2),
        point3(-1.2, 0.5, -2.6),
        shutter_open, shutter_close,
        0.5,
        material_moving));

    // BVH
    bvh_node world_bvh(world.objects);

    // ------------------------------------------------------------
    // LIGHTS LIST
    // ------------------------------------------------------------
    hittable_list lights;
    lights.add(make_shared<sphere>(light_center, light_radius, light));

    // ------------------------------------------------------------
    // CAMERA
    // ------------------------------------------------------------
    point3 lookfrom(6.0, 2.8, 6.5);
    point3 lookat(0.0, 0.4, -1.8);
    vec3 vup(0, 1, 0);
    double vfov = 30.0;

    double focus_dist = (lookfrom - lookat).length();
    double aperture = 0.6;

    camera cam(lookfrom, lookat, vup, vfov, aspect_ratio,
               aperture, focus_dist,
               shutter_open, shutter_close);

    // ------------------------------------------------------------
    // RENDER (parallel with adaptive sampling + HDR)
    // ------------------------------------------------------------
    std::ofstream out("image.ppm", std::ios::binary);
    if (!out)
    {
        std::cerr << "ERROR: Could not open image.ppm for writing.\n";
        return 1;
    }

    out << "P6\n"
        << image_width << " " << image_height << "\n255\n";

    std::ofstream out_hdr("image_hdr.pfm", std::ios::binary);
    if (!out_hdr)
    {
        std::cerr << "ERROR: Could not open image_hdr.pfm for writing.\n";
        return 1;
    }

    out_hdr << "PF\n"
            << image_width << " " << image_height << "\n-1.0\n";

    struct PixelData
    {
        color sum;
        int samples;
    };

    std::vector<PixelData> framebuffer(image_width * image_height);

    std::atomic<int> next_row(image_height - 1);
    std::atomic<long long> total_samples_used(0);

    auto render_worker = [&](int thread_id)
    {
        while (true)
        {
            int j = next_row.fetch_sub(1);
            if (j < 0)
                break;

            auto row_start = std::chrono::high_resolution_clock::now();

            for (int i = 0; i < image_width; ++i)
            {
                color pixel_sum(0, 0, 0);
                color pixel_sum_sq(0, 0, 0);
                int samples_used = 0;

                for (int s = 0; s < max_samples_per_pixel; ++s)
                {
                    auto u = (i + random_double()) / (image_width - 1);
                    auto v = (j + random_double()) / (image_height - 1);
                    ray r = cam.get_ray(u, v);

                    color sample = ray_color(r, world_bvh, lights, max_depth);

                    pixel_sum += sample;
                    pixel_sum_sq += sample * sample;
                    samples_used = s + 1;

                    if (samples_used >= min_samples_per_pixel &&
                        (samples_used % sample_batch_size) == 0)
                    {
                        double inv_n = 1.0 / samples_used;

                        color mean = pixel_sum * inv_n;
                        color mean_sq = pixel_sum_sq * inv_n;
                        color var = mean_sq - mean * mean;

                        double lum_var = 0.2126 * var.x() + 0.7152 * var.y() + 0.0722 * var.z();
                        if (lum_var < 0.0)
                            lum_var = 0.0;

                        double lum_std = std::sqrt(lum_var);

                        if (lum_std < target_noise)
                            break;
                    }
                }

                total_samples_used += samples_used;

                framebuffer[j * image_width + i].sum = pixel_sum;
                framebuffer[j * image_width + i].samples = samples_used;
            }

            auto row_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = row_end - row_start;

            std::cerr << "\rScanlines remaining: " << j
                      << " | last row time (thread " << thread_id << "): "
                      << elapsed.count() << "s " << std::flush;
        }
    };

    int num_threads = static_cast<int>(std::thread::hardware_concurrency());
    if (num_threads <= 0)
        num_threads = 4;

    std::cerr << "Launching " << num_threads << " render threads...\n";

    auto t_start = std::chrono::high_resolution_clock::now();

    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    for (int t = 0; t < num_threads; ++t)
    {
        threads.emplace_back(render_worker, t);
    }

    for (auto &t : threads)
    {
        t.join();
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_elapsed = t_end - t_start;
    std::cerr << "\nRendering finished in " << total_elapsed.count() << " seconds.\n";

    // Write out the final images from the framebuffer (single-threaded, in scanline order)
    for (int j = image_height - 1; j >= 0; --j)
    {
        for (int i = 0; i < image_width; ++i)
        {
            const PixelData &px = framebuffer[j * image_width + i];
            write_color(out, px.sum, px.samples);
            write_color_hdr(out_hdr, px.sum, px.samples);
        }
    }

    long long total_samples_used_final = total_samples_used.load();
    std::cerr << "Total samples used: " << total_samples_used_final << "\n";
    std::cerr << "Average samples per pixel: "
              << static_cast<double>(total_samples_used_final) / (image_width * image_height)
              << "\n";

    std::cerr << "Done.\n";
    return 0;
}