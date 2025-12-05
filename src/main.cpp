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

#include <iostream>
#include <chrono>

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

    // You said you’ll run width ~4000 for submission.
    // Keep a smaller one for quick tests, then crank later.
    const int image_width = 4000; // <- set to 4000 for your final
    const int image_height = static_cast<int>(image_width / aspect_ratio);

    const int samples_per_pixel = 200;
    const int max_depth = 20;

    // ------------------------------------------------------------
    // TEXTURES
    // ------------------------------------------------------------
    auto checker_tex = make_shared<checker_texture>(
        color(0.2, 0.3, 0.1),
        color(0.9, 0.9, 0.9));

    auto earth_tex = make_shared<image_texture>("earthmap.jpg");

    // ------------------------------------------------------------
    // MATERIALS
    // ------------------------------------------------------------
    auto material_ground = make_shared<lambertian>(checker_tex);

    auto material_teapot = make_shared<lambertian>(color(0.9, 0.2, 0.2));

    auto material_earth = make_shared<lambertian>(earth_tex);

    auto material_fore = make_shared<lambertian>(color(0.2, 0.8, 0.9)); // cyan
    auto material_back = make_shared<lambertian>(color(0.9, 0.9, 0.2)); // yellow

    auto material_moving = make_shared<lambertian>(color(0.7, 0.2, 0.9)); // purple moving blur

    auto light = make_shared<diffuse_light>(color(7, 7, 7));

    // ------------------------------------------------------------
    // WORLD CONTENT
    // ------------------------------------------------------------
    hittable_list world;

    // Ground (big sphere)
    world.add(make_shared<sphere>(
        point3(0.0, -100.5, -1.0),
        100.0,
        material_ground));

    // Earth-textured sphere (shows texture support)
    world.add(make_shared<sphere>(
        point3(2.2, 1.0, -3.5), // pushed back ~1 unit
        1.0,
        material_earth));

    // Teapot OBJ mesh (triangles + mesh loading)
    mesh_transform teapot_xf(
        vec3(1.4, 1.4, 1.4),
        vec3(0.0, 0.2, -1.8));

    auto teapot_tris = load_obj_mesh_as_triangles(
        "teapot.obj",
        material_teapot,
        teapot_xf);

    for (auto &t : teapot_tris)
        world.add(t);

    // Emissive light sphere overhead (area light)
    auto light_center = point3(0.0, 3.0, -1.8);
    auto light_radius = 1.5;
    world.add(make_shared<sphere>(
        light_center,
        light_radius,
        light));

    // Foreground sphere (very close) for DOF proof
    world.add(make_shared<sphere>(
        point3(-1.0, 0.6, 1.5),
        0.6,
        material_fore));

    // Background sphere (very far) for DOF proof
    world.add(make_shared<sphere>(
        point3(2.5, 1.0, -8.0),
        1.2,
        material_back));

    // Motion-blurred moving sphere (motion blur proof)
    double shutter_open = 0.0;
    double shutter_close = 1.0;

    world.add(make_shared<moving_sphere>(
        point3(-1.2, 0.5, -1.2), // start
        point3(-1.2, 0.5, -2.6), // end (longer move = clearer blur)
        shutter_open, shutter_close,
        0.5,
        material_moving));

    // Build BVH AFTER all objects are added
    bvh_node world_bvh(world.objects);

    // ------------------------------------------------------------
    // LIGHTS LIST (for importance sampling)
    // ------------------------------------------------------------
    hittable_list lights;
    lights.add(make_shared<sphere>(light_center, light_radius, light));

    // ------------------------------------------------------------
    // CAMERA (DOF + motion blur)
    // ------------------------------------------------------------
    point3 lookfrom(6.0, 2.8, 6.5);
    point3 lookat(0.0, 0.4, -1.8);
    vec3 vup(0, 1, 0);
    double vfov = 30.0;

    // Focus on teapot
    double focus_dist = (lookfrom - lookat).length();
    double aperture = 0.6; // visible DOF but not insane

    camera cam(lookfrom, lookat, vup, vfov, aspect_ratio,
               aperture, focus_dist,
               shutter_open, shutter_close);

    // ------------------------------------------------------------
    // RENDER
    // ------------------------------------------------------------
    std::cout << "P3\n"
              << image_width << " " << image_height << "\n255\n";

    for (int j = image_height - 1; j >= 0; --j)
    {
        auto start = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < image_width; ++i)
        {
            color pixel_color(0, 0, 0);

            for (int s = 0; s < samples_per_pixel; ++s)
            {
                auto u = (i + random_double()) / (image_width - 1);
                auto v = (j + random_double()) / (image_height - 1);
                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, world_bvh, lights, max_depth);
            }

            write_color(std::cout, pixel_color, samples_per_pixel);
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;

        std::cerr << "\rScanlines remaining: " << j
                  << " | last line: " << elapsed.count() << "s " << std::flush;
    }

    std::cerr << "\nDone.\n";
    return 0;
}
