#ifndef BALLTREE_H
#define BALLTREE_H

#include "point.h"

struct BallNode {
    struct Point center;
    double radius;
    struct BallNode *left;
    struct BallNode *right;
    struct PointSlice data;
};

struct BallTree {
    struct BallNode *root;
    struct PointBuffer data;
    int leafsize;
};

void balltree_free(struct BallTree*);
struct BallTree* balltree_build(const struct PointBuffer*, int);
double balltree_count_radius(const struct BallTree*, const struct Point*, double);
double balltree_count_range(const struct BallTree*, const struct Point*, double, double);
double balltree_dualcount_radius(const struct BallTree*, const struct BallTree*, double);

#endif /* BALLTREE_H */
