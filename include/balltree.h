#ifndef BALLTREE_H
#define BALLTREE_H

#include "point.h"
#include "ballnode.h"

typedef struct {
    BallNode *root;
    PointBuffer data;
    int leafsize;
} BallTree;

// from balltree.c
BallTree *balltree_build(const PointBuffer *);
BallTree *balltree_build_leafsize(const PointBuffer *, int);
void balltree_free(BallTree *);

int balltree_count_nodes(const BallTree *);
double balltree_count_radius(const BallTree *, const Point *, double);
double balltree_count_range(const BallTree *, const Point *, double, double);
double balltree_dualcount_radius(const BallTree *, const BallTree *, double);
double balltree_dualcount_range(const BallTree *, const BallTree *, double, double);

// from balltree_serialize.c
int balltree_to_file(const BallTree *, const char *);
BallTree *balltree_from_file(const char *);

#endif /* BALLTREE_H */
