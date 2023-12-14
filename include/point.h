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

void print_point(const struct Point*);
double points_distance(const struct Point*, const struct Point*);
double points_distance2(const struct Point*, const struct Point*);

struct PointBuffer pointbuffer_create(int);
int pointbuffer_resize(struct PointBuffer*, int);

struct PointSlice pointslice_from_buffer(const struct PointBuffer);
void print_pointslice(const struct PointSlice*);
size_t get_pointslice_size(const struct PointSlice*);
struct Point get_center_point(const struct PointSlice*);
double get_maxdist_from_center(const struct PointSlice*, struct Point);
enum Axis get_max_spread_axis(const struct PointSlice*);
int partial_median_sort(struct PointSlice*, enum Axis);


#endif /* POINT_H */
