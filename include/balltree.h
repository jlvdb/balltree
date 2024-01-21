#ifndef BALLTREE_H
#define BALLTREE_H

#include "point.h"
#include "ballnode.h"

#define DEFAULT_LEAFSIZE 20

typedef struct {
    BallNode *root;
    PointBuffer data;
    int leafsize;
    int data_owned;
} BallTree;

// from balltree.c
BallTree *balltree_build(const PointBuffer *);
BallTree *balltree_build_leafsize(const PointBuffer *, int);
BallTree* balltree_build_nocopy(PointBuffer *, int);
void balltree_free(BallTree *);

double balltree_count_radius(const BallTree *, const Point *, double);
double balltree_count_range(const BallTree *, const Point *, double, double);
double balltree_dualcount_radius(const BallTree *, const BallTree *, double);
double balltree_dualcount_range(const BallTree *, const BallTree *, double, double);

int balltree_count_nodes(const BallTree *);
StatsVector *balltree_collect_stats(const BallTree *);

// from balltree_serialize.c
int balltree_to_file(const BallTree *, const char *);
BallTree *balltree_from_file(const char *);

#endif /* BALLTREE_H */
