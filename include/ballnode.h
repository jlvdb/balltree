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


typedef struct {
    double center_x, center_y, center_z;
    double radius;
    double sum_weight;
    int left, right;  // BallNode* -> index in BallNodeBuffer
    int data_start, data_end;  // PointSlice
} BallNodeSerialized;

typedef struct {
    int size;
    BallNodeSerialized *buffer;
    int *next_free;
} BallNodeBuffer;

BallNode *ballnode_build_recursive(PointSlice *, int);
int ballnode_count_nodes(const BallNode *);
double ballnode_count_radius(const BallNode *, const Point *, double);
double ballnode_count_range(const BallNode *, const Point *, double, double);
double ballnode_dualcount_radius(const BallNode *, const BallNode *, double);
int ballnode_serialise_recursive(BallNodeBuffer, BallNode *, int);
BallNode *ballnode_deserialise_recursive(BallNodeSerialized *, int, const PointBuffer *, int);

#endif /* BALLNODE_H */
