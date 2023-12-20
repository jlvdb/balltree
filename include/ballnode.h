#ifndef BALLNODE_H
#define BALLNODE_H

#include "point.h"

struct BallNode {
    struct Point center;
    double radius;
    struct BallNode *left;
    struct BallNode *right;
    struct PointSlice data;
};

struct BallNodeSerialized {
    double center_x, center_y, center_z;
    double radius;
    int left, right;  // struct BallNode*
    int data_start, data_end;  // struct PointSlice
};

struct BallNodeBuffer {
    int size;
    int *next_free;
    struct BallNodeSerialized *buffer;
};

struct BallNode* ballnode_build_recursive(struct PointSlice*, int);
int ballnode_count_nodes(const struct BallNode*);
double ballnode_count_radius(const struct BallNode*, const struct Point*, double);
double ballnode_count_range(const struct BallNode*, const struct Point*, double, double);
double ballnode_dualcount_radius(const struct BallNode*, const struct BallNode*, double);
int ballnode_serialise_recursive(struct BallNodeBuffer, struct BallNode*, int);
struct BallNode* ballnode_deserialise_recursive(struct BallNodeSerialized*, int, const struct PointBuffer*, int);

#endif /* BALLNODE_H */
