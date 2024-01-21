#ifndef BALLNODE_H
#define BALLNODE_H

#include "point.h"

typedef struct {
    double x;
    double y;
    double z;
    double radius;
} Ball;

struct BallNode;

typedef struct {
    struct BallNode *left;
    struct BallNode *right;
} Childs;

struct BallNode {
    Ball ball;
    union {
        Childs childs;    // is_leaf == 0
        PointSlice data;  // is_leaf == 1
    };
    struct {
        unsigned long is_leaf : 1;
        unsigned long num_points : 63;
    };
    double sum_weight;
};
typedef struct BallNode BallNode;

#define BALLNODE_IS_LEAF(node) (node)->is_leaf

typedef struct {
    long depth;
    long num_points;
    double sum_weight;
    double x, y, z;
    double radius;
} NodeStats;

typedef struct {
    NodeStats *stats;
    int size;
    int end;
} StatsVector;

// from ballnode.c
BallNode *bnode_build(PointSlice *, int);
void bnode_free(BallNode *);
int bnode_is_leaf(const BallNode *);
PointSlice bnode_get_ptslc(const BallNode *);

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
