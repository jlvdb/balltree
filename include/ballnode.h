#ifndef BALLNODE_H
#define BALLNODE_H

#include <stdint.h>

#include "point.h"
#include "containers.h"

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
        uint64_t is_leaf : 1;
        uint64_t num_points : 63;
    };
    double sum_weight;
};
typedef struct BallNode BallNode;

typedef struct {
    int64_t depth;
    int64_t num_points;
    double sum_weight;
    double x, y, z;
    double radius;
} NodeStats;

typedef struct {
    NodeStats *stats;
    int64_t capacity;
    int64_t size;
} StatsVector;

// from ballnode.c
BallNode *bnode_build(PointSlice *, int);
void bnode_free(BallNode *);
PointSlice bnode_get_ptslc(const BallNode *);

// from ballnode_query.c
void bnode_nearest_neighbours(const BallNode *, const Point *, KnnQueue *);
double bnode_count_radius(const BallNode *, const Point *, double);
void bnode_count_range(const BallNode *, const Point *, DistHistogram *);
double bnode_dualcount_radius(const BallNode *, const BallNode *, double);
void bnode_dualcount_range(const BallNode *, const BallNode *, DistHistogram *);

// from ballnode_stats.c
StatsVector *statvec_new(void);
StatsVector *statvec_new_reserve(int64_t);
void statvec_free(StatsVector *);
int statvec_resize(StatsVector *, int64_t);
int statvec_append(StatsVector *, const NodeStats *);
int bnode_collect_stats(const BallNode *, StatsVector *, int);
int64_t bnode_count_nodes(const BallNode *);

#endif /* BALLNODE_H */
