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
    return sqrt(point_dist_sq(p1, p2));
}

double point_dist_sq(const Point *p1, const Point *p2) {
    double dx = p1->x - p2->x;
    double dy = p1->y - p2->y;
    double dz = p1->z - p2->z;
    return dx*dx + dy*dy + dz*dz;
}
