#ifndef BALLTREE_H
#define BALLTREE_H

#include "point.h"

struct BallTree {
    struct Point center;
    double radius;
    struct BallTree *left;
    struct BallTree *right;
    struct PointBuffer data;
};

int balltree_is_leaf(const struct BallTree*);
void balltree_free(struct BallTree*);
void balltree_print(const struct BallTree*);
struct BallTree* balltree_build_recursive(struct PointSlice*, int);
struct BallTree* balltree_build(struct PointBuffer*, int);
double balltree_count_radius(struct BallTree*, struct Point*, double);
double balltree_count_range(struct BallTree*, struct Point*, double, double);

#endif /* BALLTREE_H */
