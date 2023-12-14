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

int balltree_is_leaf(const struct BallTree *node);
void balltree_print(const struct BallTree *node);
struct BallTree* balltree_build(struct PointSlice *points, int leafsize);
void balltree_free(struct BallTree *node);

#endif /* BALLTREE_H */
