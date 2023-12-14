#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "point.h"

#define POINT_ACCESS_BY_INDEX(ptr, index) (*((double*)((char*)(ptr) + (index) * sizeof(double))))
#define SWAP(temp, a, b) do { (temp) = (a); (a) = (b); (b) = (temp); } while (0)

void print_point(const struct Point *point)
{
    printf("{%lf, %lf, %lf}\n", point->x, point->y, point->z);
}

void swap_points(struct Point *p1, struct Point *p2)
{
    double temp;
    SWAP(temp, p1->x, p2->x);
    SWAP(temp, p1->y, p2->y);
    SWAP(temp, p1->z, p2->z);
}

double points_distance2(const struct Point *p1, const struct Point *p2)
{
    double dx = p1->x - p2->x;
    double dy = p1->y - p2->y;
    double dz = p1->z - p2->z;
    return dx*dx + dy*dy + dz*dz;
}

double points_distance(const struct Point *p1, const struct Point *p2)
{
    return sqrt(points_distance2(p1, p2));
}

void print_pointslice(const struct PointSlice *slice)
{
    printf("{\n");
    for (size_t i = slice->start; i < slice->end; i++) {
        printf("    ");
        print_point(&slice->points[i]);
    }
    printf("}\n");
}

size_t get_pointslice_size(const struct PointSlice *slice)
{
    return slice->end - slice->start;
}

struct Point get_center_point(const struct PointSlice *slice)
{
    double center_x = 0.0;
    double center_y = 0.0;
    double center_z = 0.0;
    int total = 1;

    struct Point *points = slice->points;
    for (size_t i = slice->start; i < slice->end; i++, total++) {
        struct Point point = points[i];
        double scale = (double)total;
        center_x += (point.x - center_x) / scale;
        center_y += (point.y - center_y) / scale;
        center_z += (point.z - center_z) / scale;
    }

    struct Point center = {
        .x = center_x,
        .y = center_y,
        .z = center_z
    };
    return center;
}

double get_maxdist_from_center(const struct PointSlice *slice, struct Point center)
{
    double maxdist = 0.0;
    struct Point *points = slice->points;
    for (size_t i = slice->start; i < slice->end; i++) {
        double points_distance = points_distance2(&points[i], &center);
        if (points_distance > maxdist) {
            maxdist = points_distance;
        }
    }
    return sqrt(maxdist);
}

enum Axis get_max_spread_axis(const struct PointSlice *slice)
{
    double x_min = INFINITY;
    double y_min = INFINITY;
    double z_min = INFINITY;
    double x_max = -INFINITY;
    double y_max = -INFINITY;
    double z_max = -INFINITY;

    struct Point *points = slice->points;
    double xi, yi, zi;
    for (size_t i = slice->start; i < slice->end; i++) {
        struct Point point = points[i];

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
    enum Axis axis;
    if (x_spread > y_spread && x_spread > z_spread) {
        axis = X;
    } else if (y_spread > z_spread) {
        axis = Y;
    } else {
        axis = Z;
    }
    return axis;
}

size_t partition_points(struct PointSlice *slice, size_t i_pivot, enum Axis axis)
{
    struct Point *points = slice->points;
    size_t i_last = slice->end - 1;

    double pivot = POINT_ACCESS_BY_INDEX(&points[i_pivot], axis);
    swap_points(&points[i_pivot], &points[i_last]);

    size_t i_partition = slice->start;
    for (size_t i = slice->start; i < i_last; i++) {
        if (POINT_ACCESS_BY_INDEX(&points[i], axis) < pivot) {
            if (i_partition != i) {
                swap_points(&points[i], &points[i_partition]);
            }
            i_partition++;
        }
    }

    swap_points(&points[i_last], &points[i_partition]);
    return i_partition;
}

int quickselect(struct PointSlice *slice, size_t k, enum Axis axis)
{
    if (slice->start < slice->end) {
        int i_pivot = (slice->start + slice->end) / 2;
        i_pivot = partition_points(slice, i_pivot, axis);

        if (i_pivot < k) {
            struct PointSlice subslice = {
                .start = i_pivot + 1,
                .end = slice->end,
                .points = slice->points,
            };
            i_pivot = quickselect(&subslice, k, axis);
        } else if (i_pivot > k) {
            struct PointSlice subslice = {
                .start = slice->start,
                .end = i_pivot,
                .points = slice->points,
            };
            i_pivot = quickselect(&subslice, k, axis);
        }
        return i_pivot;
    } else {
        return -1;
    }
}

int partial_median_sort(struct PointSlice *slice, enum Axis axis)
{
    size_t i_median = (slice->end + slice->start) / 2;
    return quickselect(slice, i_median, axis);
}
