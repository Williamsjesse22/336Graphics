#ifndef BVH_NODE_H
#define BVH_NODE_H

#include "rtweekend.h"
#include "hittable.h"
#include "aabb.h"

#include <algorithm>
#include <vector>

class bvh_node : public hittable
{
public:
    shared_ptr<hittable> left;
    shared_ptr<hittable> right;
    aabb box;

    bvh_node() {}

    bvh_node(const std::vector<shared_ptr<hittable>> &src_objects,
             size_t start, size_t end)
    {
        auto objects = src_objects; // make a modifiable copy

        int axis = static_cast<int>(3 * random_double()); // 0,1,2

        auto comparator = (axis == 0)   ? box_x_compare
                          : (axis == 1) ? box_y_compare
                                        : box_z_compare;

        size_t object_span = end - start;

        if (object_span == 1)
        {
            left = right = objects[start];
        }
        else if (object_span == 2)
        {
            if (comparator(objects[start], objects[start + 1]))
            {
                left = objects[start];
                right = objects[start + 1];
            }
            else
            {
                left = objects[start + 1];
                right = objects[start];
            }
        }
        else
        {
            std::sort(objects.begin() + start, objects.begin() + end, comparator);

            auto mid = start + object_span / 2;
            left = make_shared<bvh_node>(objects, start, mid);
            right = make_shared<bvh_node>(objects, mid, end);
        }

        aabb box_left, box_right;

        if (!left->bounding_box(box_left) || !right->bounding_box(box_right))
            std::cerr << "ERROR: No bounding box in bvh_node constructor.\n";

        box = surrounding_box(box_left, box_right);
    }

    // convenience ctor
    bvh_node(const std::vector<shared_ptr<hittable>> &objects)
        : bvh_node(objects, 0, objects.size()) {}

    bool hit(const ray &r, interval ray_t, hit_record &rec) const override
    {
        if (!box.hit(r, ray_t))
            return false;

        bool hit_left = left->hit(r, ray_t, rec);
        bool hit_right = right->hit(r, interval(ray_t.min, hit_left ? rec.t : ray_t.max), rec);

        return hit_left || hit_right;
    }

    bool bounding_box(aabb &output_box) const override
    {
        output_box = box;
        return true;
    }

private:
    static bool box_compare(const shared_ptr<hittable> a,
                            const shared_ptr<hittable> b, int axis)
    {
        aabb box_a;
        aabb box_b;

        if (!a->bounding_box(box_a) || !b->bounding_box(box_b))
            std::cerr << "ERROR: No bounding box in bvh_node comparison.\n";

        return box_a.axis_interval(axis).min < box_b.axis_interval(axis).min;
    }

    static bool box_x_compare(const shared_ptr<hittable> a,
                              const shared_ptr<hittable> b)
    {
        return box_compare(a, b, 0);
    }

    static bool box_y_compare(const shared_ptr<hittable> a,
                              const shared_ptr<hittable> b)
    {
        return box_compare(a, b, 1);
    }

    static bool box_z_compare(const shared_ptr<hittable> a,
                              const shared_ptr<hittable> b)
    {
        return box_compare(a, b, 2);
    }
};

#endif
