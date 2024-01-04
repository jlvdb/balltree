#include <math.h>

#include "point.h"

inline void point_swap(Point *p1, Point *p2);
inline double point_get_coord(const Point *point, enum Axis axis);
int pointslice_partition(PointSlice *slice, int pivot_idx, enum Axis axis);
int pointslice_quickselect(PointSlice *slice, int partition_idx, enum Axis axis);

inline void point_swap(Point *p1, Point *p2) {
    Point temp = *p1;
    *p1 = *p2;
    *p2 = temp;
}

inline double point_get_coord(const Point *point, enum Axis axis) {
    return *((double*)point + axis);
}

Point pointslice_get_mean(const PointSlice *slice) {
    double center_x = 0.0;
    double center_y = 0.0;
    double center_z = 0.0;
    int total = 0;

    Point *points = slice->points;
    for (int i = slice->start; i < slice->end; ++i) {
        ++total;
        double scale = (double)total;
        Point point = points[i];
        center_x += (point.x - center_x) / scale;
        center_y += (point.y - center_y) / scale;
        center_z += (point.z - center_z) / scale;
    }
    return point_create(center_x, center_y, center_z);
}

double pointslice_get_maxdist(const PointSlice *slice, Point *center) {
    Point *points = slice->points;
    double dist_squared_max = 0.0;
    for (int i = slice->start; i < slice->end; ++i) {
        double dist_squared = point_dist_squared(points + i, center);
        if (dist_squared > dist_squared_max) {
            dist_squared_max = dist_squared;
        }
    }
    return sqrt(dist_squared_max);
}

enum Axis pointslice_get_maxspread_axis(const PointSlice *slice) {
    double x_min = INFINITY;
    double y_min = x_min;
    double z_min = x_min;

    double x_max = -INFINITY;
    double y_max = x_max;
    double z_max = x_max;

    Point *points = slice->points;
    double xi, yi, zi;
    for (int i = slice->start; i < slice->end; ++i) {
        Point point = points[i];

        xi = point.x;
        if (xi < x_min) {
            x_min = xi;
        } else if (xi > x_max) {
            x_max = xi;
        }
        yi = point.y;
        if (yi < y_min) {
            y_min = yi;
        } else if (yi > y_max) {
            y_max = yi;
        }
        zi = point.z;
        if (zi < z_min) {
            z_min = zi;
        } else if (zi > z_max) {
            z_max = zi;
        }
    }

    double x_spread = x_max - x_min;
    double y_spread = y_max - y_min;
    double z_spread = z_max - z_min;
    if (x_spread > y_spread && x_spread > z_spread) {
        return (enum Axis)X;
    } else if (y_spread > z_spread) {
        return (enum Axis)Y;
    } else {
        return (enum Axis)Z;
    }
}

int pointslice_partition(PointSlice *slice, int pivot_idx, enum Axis axis) {
    Point *points = slice->points;
    int last_idx = slice->end - 1;

    double pivot = point_get_coord(points + pivot_idx, axis);
    point_swap(points + pivot_idx, points + last_idx);

    int partition_idx = slice->start;
    for (int i = slice->start; i < last_idx; ++i) {
        if (point_get_coord(points + i, axis) < pivot) {
            if (partition_idx != i) {
                point_swap(points + i, points + partition_idx);
            }
            ++partition_idx;
        }
    }

    point_swap(points + last_idx, points + partition_idx);
    return partition_idx;
}

int pointslice_quickselect(PointSlice *slice, int partition_idx, enum Axis axis) {
    if (slice->start < slice->end) {
        int pivot_idx = (slice->start + slice->end) / 2;
        pivot_idx = pointslice_partition(slice, pivot_idx, axis);

        if (pivot_idx < partition_idx) {
            PointSlice subslice = {
                .start = pivot_idx + 1,
                .end = slice->end,
                .points = slice->points,
            };
            pivot_idx = pointslice_quickselect(&subslice, partition_idx, axis);
        } else if (pivot_idx > partition_idx) {
            PointSlice subslice = {
                .start = slice->start,
                .end = pivot_idx,
                .points = slice->points,
            };
            pivot_idx = pointslice_quickselect(&subslice, partition_idx, axis);
        }
        return pivot_idx;
    } else {
        return -1;
    }
}

int pointslice_median_partition(PointSlice *slice, enum Axis axis) {
    int median_idx = (slice->end + slice->start) / 2;
    return pointslice_quickselect(slice, median_idx, axis);
}
