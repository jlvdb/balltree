#include <math.h>

#include "point.h"

Point point_create(double x, double y, double z) {
    return point_create_weighted(x, y, z, 1.0);
}

Point point_create_weighted(double x, double y, double z, double weight) {
    return (Point){
        .x = x,
        .y = y,
        .z = z,
        .weight = weight,
    };
}

double point_dist(const Point *p1, const Point *p2) {
    return sqrt(EUCLIDEAN_DIST_SQ(p1, p2));
}

double point_dist_sq(const Point *p1, const Point *p2) {
    return EUCLIDEAN_DIST_SQ(p1, p2);
}
