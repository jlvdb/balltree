#ifndef POINT_H
#define POINT_H

enum Axis {X, Y, Z};

typedef struct {
    double x;
    double y;
    double z;
    double weight;
} Point;

typedef struct {
    int size;
    Point *points;
} PointBuffer;

typedef struct {
    int start;
    int end;
    Point *points;
} PointSlice;

// from point.c
Point point_create(double, double, double);
Point point_create_weighted(double, double, double, double);
double point_dist(const Point *, const Point *);
double point_dist_sq(const Point *, const Point *);

// from pointbuffers.c
PointBuffer *ptbuf_new(int);
void ptbuf_free(PointBuffer *);
int ptbuf_resize(PointBuffer *, int);
PointBuffer *ptbuf_copy(const PointBuffer *);
PointSlice *ptslc_from_buffer(const PointBuffer *);
int ptslc_get_size(const PointSlice *);

#endif /* POINT_H */
