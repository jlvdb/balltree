#ifndef BALLNODE_H
#define BALLNODE_H

#include "point.h"

struct BallNode {
    Point center;
    double radius;
    struct BallNode *left;
    struct BallNode *right;
    PointSlice data;
    double sum_weight;
};
typedef struct BallNode BallNode;

BallNode *bnode_build_recursive(PointBuffer *, int, int, int);
void bnode_free(BallNode *);
int bnode_is_leaf(const BallNode *);

int bnode_count_nodes(const BallNode *);
double bnode_count_radius(const BallNode *, const Point *, double);
double bnode_count_range(const BallNode *, const Point *, double, double);
double bnode_dualcount_radius(const BallNode *, const BallNode *, double);
double bnode_dualcount_range(const BallNode *, const BallNode *, double, double);

#endif /* BALLNODE_H */
