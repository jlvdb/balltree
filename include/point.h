#ifndef POINT_H
#define POINT_H

#include <stdlib.h>

struct Point {
    double x;
    double y;
    double z;
};

enum Axis {X, Y, Z};  // field to index mapping for struct Point

struct PointBuffer {
    size_t size;
    struct Point *points;
};

struct PointSlice {
    size_t start;
    size_t end;
    struct Point *points;
};

void print_point(const struct Point *point);
double points_distance(const struct Point *p1, const struct Point *p2);

struct PointBuffer pointbuffer_create(int size);
int pointbuffer_resize(struct PointBuffer *buffer, int newsize);

struct PointSlice pointslice_from_buffer(const struct PointBuffer buffer);
void print_pointslice(const struct PointSlice *slice);
size_t get_pointslice_size(const struct PointSlice *slice);
struct Point get_center_point(const struct PointSlice *slice);
double get_maxdist_from_center(const struct PointSlice *slice, struct Point center);
enum Axis get_max_spread_axis(const struct PointSlice *slice);
int partial_median_sort(struct PointSlice *slice, enum Axis axis);


#endif /* POINT_H */
