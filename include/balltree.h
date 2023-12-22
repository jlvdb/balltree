#ifndef BALLTREE_H
#define BALLTREE_H

#include "point.h"
#include "ballnode.h"

struct BallTree {
    struct BallNode *root;
    struct PointBuffer data;
    int leafsize;
};

struct BallTree* balltree_build(const struct PointBuffer*);
struct BallTree* balltree_build_leafsize(const struct PointBuffer*, int);
void balltree_free(struct BallTree*);
int balltree_count_nodes(const struct BallTree*);
int balltree_to_file(const struct BallTree*, const char*);
struct BallTree* balltree_from_file(const char*);

double balltree_count_radius(const struct BallTree*, const struct Point*, double);
double balltree_count_range(const struct BallTree*, const struct Point*, double, double);
double balltree_dualcount_radius(const struct BallTree*, const struct BallTree*, double);

#endif /* BALLTREE_H */
