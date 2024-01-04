#ifndef POINT_H
#define POINT_H

#include <stdlib.h>
#include <stdio.h>

typedef struct {
    double x;
    double y;
    double z;
    double weight;
} Point;

enum Axis {X, Y, Z};

Point point_create_weighted(double, double, double, double);
Point point_create(double, double, double);
double point_dist_squared(const Point *, const Point *);
double point_dist(const Point *, const Point *);

typedef struct {
    int size;
    Point *points;
} PointBuffer;

PointBuffer *pointbuffer_new(int);
void pointbuffer_free(PointBuffer *);
int pointbuffer_resize(PointBuffer *, int);

typedef struct {
    int start;
    int end;
    Point *points;
} PointSlice;

PointSlice *pointslice_from_buffer(const PointBuffer *);
int pointslice_get_size(const PointSlice *);
int pointslice_update_start(PointSlice *, int);
int pointslice_update_end(PointSlice *, int);

Point pointslice_get_mean(const PointSlice *);
double pointslice_get_maxdist(const PointSlice *, Point *);
enum Axis pointslice_get_maxspread_axis(const PointSlice *);
int pointslice_median_partition(PointSlice *, enum Axis);

#endif /* POINT_H */
