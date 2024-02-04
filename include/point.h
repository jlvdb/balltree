#ifndef POINT_H
#define POINT_H

#include <math.h>
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

#define POINT_CREATE(x, y, z) (Point){(x), (y), (z), .weight = 1.0, .index = 0}
#define POINT_CREATE_WEIGHTED(x, y, z, weight) (Point){(x), (y), (z), .weight = (weight), .index = 0}

inline double point_dist(const Point *p1, const Point *p2) {
    return sqrt(EUCLIDEAN_DIST_SQ(p1, p2));
}

inline double point_dist_sq(const Point *p1, const Point *p2) {
    return EUCLIDEAN_DIST_SQ(p1, p2);
}

// from pointbuffers.c
PointBuffer *ptbuf_new(int64_t);
void ptbuf_free(PointBuffer *);
int ptbuf_resize(PointBuffer *, int64_t);
PointBuffer *ptbuf_copy(const PointBuffer *);
PointBuffer *ptbuf_gen_random(double, double, int64_t);

PointSlice *ptslc_from_buffer(const PointBuffer *);
double ptslc_sumw_in_radius_sq(const PointSlice *, const Point *, double);
double ptslc_sumw_in_range_sq(const PointSlice *, const Point *, double, double);
double ptslc_dualsumw_in_radius_sq(const PointSlice *, const PointSlice *, double);
double ptslc_dualsumw_in_range_sq(const PointSlice *, const PointSlice *, double, double);

inline int64_t ptslc_get_size(const PointSlice *slice) {
    return slice->end - slice->start;
}

#endif /* POINT_H */
