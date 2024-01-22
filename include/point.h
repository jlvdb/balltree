#ifndef POINT_H
#define POINT_H

#define EUCLIDEAN_DIST_SQ(p1, p2) \
    ({ \
        double __dx = (p1)->x - (p2)->x; \
        double __dy = (p1)->y - (p2)->y; \
        double __dz = (p1)->z - (p2)->z; \
        __dx * __dx + __dy * __dy + __dz * __dz; \
    })

enum Axis {X, Y, Z};

typedef struct {
    double x;
    double y;
    double z;
    double weight;
} Point;

typedef struct {
    Point *points;
    int size;
} PointBuffer;

typedef struct {
    Point *start;
    Point *end;
} PointSlice;

// from point.c
Point point_create(double, double, double);
Point point_create_weighted(double, double, double, double);
double point_dist(const Point *, const Point *);
double point_dist_sq(const Point *, const Point *);

// from pointbuffers.c
PointBuffer *ptbuf_new(int);
PointBuffer *ptbuf_from_buffers(int, double *, double *, double *);
PointBuffer *ptbuf_from_buffers_weighted(int, double *, double *, double *, double *);
void ptbuf_free(PointBuffer *);
int ptbuf_resize(PointBuffer *, int);
PointBuffer *ptbuf_copy(const PointBuffer *);
PointBuffer *ptbuf_gen_random(double, double, int);

PointSlice *ptslc_from_buffer(const PointBuffer *);
int ptslc_get_size(const PointSlice *);

#endif /* POINT_H */
