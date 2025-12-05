#ifndef HITTABLE_LIST_H
#define HITTABLE_LIST_H

#include "hittable.h"
#include "aabb.h"

#include <vector>

class hittable_list : public hittable
{
public:
    std::vector<shared_ptr<hittable>> objects;

    hittable_list() {}
    hittable_list(shared_ptr<hittable> object) { add(object); }

    void clear() { objects.clear(); }

    void add(shared_ptr<hittable> object)
    {
        objects.push_back(object);
    }

    bool hit(const ray &r, interval ray_t, hit_record &rec) const override
    {
        hit_record temp_rec;
        bool hit_anything = false;
        auto closest_so_far = ray_t.max;

        for (const auto &object : objects)
        {
            if (object->hit(r, interval(ray_t.min, closest_so_far), temp_rec))
            {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                rec = temp_rec;
            }
        }

        return hit_anything;
    }

    bool bounding_box(aabb &output_box) const override
    {
        if (objects.empty())
            return false;

        aabb temp_box;
        bool first_box = true;

        for (const auto &object : objects)
        {
            if (!object->bounding_box(temp_box))
                return false;
            output_box = first_box ? temp_box : surrounding_box(output_box, temp_box);
            first_box = false;
        }

        return true;
    }
    double pdf_value(const point3 &o, const vec3 &v) const override
    {
        auto weight = 1.0 / objects.size();
        auto sum = 0.0;

        for (const auto &object : objects)
            sum += weight * object->pdf_value(o, v);

        return sum;
    }

    vec3 random(const point3 &o) const override
    {
        auto int_size = static_cast<int>(objects.size());
        return objects[random_int(0, int_size - 1)]->random(o);
    }
};

#endif
