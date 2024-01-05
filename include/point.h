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
double point_dist_sq(const Point *, const Point *);
double point_dist(const Point *, const Point *);

typedef struct {
    int size;
    Point *points;
} PointBuffer;

PointBuffer *ptbuf_new(int);
void ptbuf_free(PointBuffer *);
int ptbuf_resize(PointBuffer *, int);
PointBuffer *ptbuf_copy(const PointBuffer *);

typedef struct {
    int start;
    int end;
    Point *points;
} PointSlice;

PointSlice *ptslc_from_buffer(const PointBuffer *);
int ptslc_get_size(const PointSlice *);

#endif /* POINT_H */
