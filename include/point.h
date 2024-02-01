#ifndef POINT_H
#define POINT_H

#include <stdint.h>

enum Axis {X, Y, Z};

typedef struct {
    double x;
    double y;
    double z;
    double weight;
    int64_t index;
} Point;

typedef struct {
    Point *points;
    int64_t size;
} PointBuffer;

typedef struct {
    Point *start;
    Point *end;
} PointSlice;

#define EUCLIDEAN_DIST_SQ(p1, p2) \
    (((p1)->x - (p2)->x) * ((p1)->x - (p2)->x) + \
     ((p1)->y - (p2)->y) * ((p1)->y - (p2)->y) + \
     ((p1)->z - (p2)->z) * ((p1)->z - (p2)->z))

// from point.c
Point point_create(double, double, double);
Point point_create_weighted(double, double, double, double);
double point_dist(const Point *, const Point *);
double point_dist_sq(const Point *, const Point *);

// from pointbuffers.c
PointBuffer *ptbuf_new(int64_t);
void ptbuf_free(PointBuffer *);
int ptbuf_resize(PointBuffer *, int64_t);
PointBuffer *ptbuf_copy(const PointBuffer *);
PointBuffer *ptbuf_gen_random(double, double, int64_t);

PointSlice *ptslc_from_buffer(const PointBuffer *);
int64_t ptslc_get_size(const PointSlice *);
double ptslc_sumw_in_radius_sq(const PointSlice *, const Point *, double);
double ptslc_sumw_in_range_sq(const PointSlice *, const Point *, double, double);
double ptslc_dualsumw_in_radius_sq(const PointSlice *, const PointSlice *, double);
double ptslc_dualsumw_in_range_sq(const PointSlice *, const PointSlice *, double, double);

#endif /* POINT_H */
