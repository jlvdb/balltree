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

// from ballnode.c
BallNode *bnode_build(PointBuffer *, int, int, int);
void bnode_free(BallNode *);
int bnode_is_leaf(const BallNode *);

// from ballnode_query.c
int bnode_count_nodes(const BallNode *);
double bnode_count_radius(const BallNode *, const Point *, double);
double bnode_count_range(const BallNode *, const Point *, double, double);
double bnode_dualcount_radius(const BallNode *, const BallNode *, double);
double bnode_dualcount_range(const BallNode *, const BallNode *, double, double);

#endif /* BALLNODE_H */
