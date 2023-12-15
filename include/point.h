#ifndef POINT_H
#define POINT_H

#include <stdlib.h>

struct Point {
    double x;
    double y;
    double z;
    double weight;
};

enum Axis {X, Y, Z};  // field to index mapping for struct Point

struct PointBuffer {
    int size;
    struct Point *points;
};

struct PointSlice {
    int start;
    int end;
    struct Point *points;
};

struct Point create_point_weighted(double, double, double, double);
struct Point create_point_unweighted(double, double, double);
void print_point(const struct Point*);
double points_distance(const struct Point*, const struct Point*);
double points_distance2(const struct Point*, const struct Point*);

struct PointBuffer* pointbuffer_create(int);
void pointbuffer_free(struct PointBuffer*);
int pointbuffer_resize(struct PointBuffer*, int);
void print_pointbuffer(const struct PointBuffer*);
double count_within_radius(struct PointBuffer*, struct Point*, double);
double count_within_range(struct PointBuffer*, struct Point*, double, double);

struct PointSlice* pointslice_from_buffer(const struct PointBuffer*);
void pointslice_free(struct PointSlice*);
void print_pointslice(const struct PointSlice*);
int get_pointslice_size(const struct PointSlice*);
struct Point get_center_point(const struct PointSlice*);
double get_maxdist_from_center(const struct PointSlice*, struct Point*);
enum Axis get_max_spread_axis(const struct PointSlice*);
int partial_median_sort(struct PointSlice *, enum Axis);

#endif /* POINT_H */
