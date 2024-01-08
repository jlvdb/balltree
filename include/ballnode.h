#ifndef BALLNODE_H
#define BALLNODE_H

#include "point.h"

struct BallNode {
    Point center;
    double radius;
    double sum_weight;
    PointSlice data;
    struct BallNode *left;
    struct BallNode *right;
};
typedef struct BallNode BallNode;

typedef struct {
    int depth;
    int n_points;
    double sum_weight;
    double x, y, z;
    double radius;
} NodeStats;

typedef struct {
    int size;
    int end;
    NodeStats *stats;
} StatsVector;

// from ballnode.c
BallNode *bnode_build(PointBuffer *, int, int, int);
void bnode_free(BallNode *);
int bnode_is_leaf(const BallNode *);

// from ballnode_query.c
double bnode_count_radius(const BallNode *, const Point *, double);
double bnode_count_range(const BallNode *, const Point *, double, double);
double bnode_dualcount_radius(const BallNode *, const BallNode *, double);
double bnode_dualcount_range(const BallNode *, const BallNode *, double, double);

// from ballnode_stats.c
StatsVector *statvec_new_reserve(int);
void statvec_free(StatsVector *);
int statvec_resize(StatsVector *, int);
int statvec_append(StatsVector *, const NodeStats *);
int bnode_collect_stats(const BallNode *, StatsVector *, int);
int bnode_count_nodes(const BallNode *);

#endif /* BALLNODE_H */
