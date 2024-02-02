#ifndef BALLTREE_H
#define BALLTREE_H

#include <stdint.h>

#include "point.h"
#include "containers.h"
#include "ballnode.h"

#define DEFAULT_LEAFSIZE 20

typedef struct {
    BallNode *root;
    PointBuffer *data;
    int leafsize;
    int data_owned;
} BallTree;

// from balltree.c
BallTree *balltree_build(const PointBuffer *);
BallTree *balltree_build_leafsize(const PointBuffer *, int);
BallTree* balltree_build_nocopy(PointBuffer *, int);
void balltree_free(BallTree *);

KnnQueue *balltree_nearest_neighbours(const BallTree *, const Point *, int64_t, double);
double balltree_brute_radius(const BallTree *, const Point *, double);
void balltree_brute_range(const BallTree *, const Point *, DistHistogram *);
double balltree_count_radius(const BallTree *, const Point *, double);
void balltree_count_range(const BallTree *, const Point *, DistHistogram *);
double balltree_dualcount_radius(const BallTree *, const BallTree *, double);
void balltree_dualcount_range(const BallTree *, const BallTree *, DistHistogram *);

int64_t balltree_count_nodes(const BallTree *);
StatsVector *balltree_collect_stats(const BallTree *);

// from balltree_serialize.c
int balltree_to_file(const BallTree *, const char *);
BallTree *balltree_from_file(const char *);

#endif /* BALLTREE_H */
