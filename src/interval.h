#ifndef INTERVAL_H
#define INTERVAL_H

#include <limits>

class interval
{
public:
    double min, max;

    interval()
        : min(std::numeric_limits<double>::infinity()),
          max(-std::numeric_limits<double>::infinity()) {}

    interval(double _min, double _max)
        : min(_min), max(_max) {}

    bool contains(double x) const
    {
        return min <= x && x <= max;
    }

    bool surrounds(double x) const
    {
        return min < x && x < max;
    }
};

#endif
