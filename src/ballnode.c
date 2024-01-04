#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "point.h"
#include "ballnode.h"

#define SUCCESS 1
#define FAILED  0

#define TRUE  1
#define FALSE 0

BallNode *ballnode_create_node(Point center, double radius, const PointSlice *slice) {
    BallNode *node = (BallNode *)malloc(sizeof(BallNode));
    if (!node) {
        fprintf(stderr, "ERROR: BallNode node allocation failed\n");
        return NULL;
    }
    node->center = center;
    node->radius = radius;
    node->sum_weight = 0.0;
    node->left = NULL;
    node->right = NULL;
    node->data = (PointSlice){
        .start = slice->start,
        .end = slice->end,
        .points = slice->points,
    };
    return node;
}

inline int ballnode_is_leaf(const BallNode *node) {
    if (node->left && node->right) {
        return FALSE;
    }
    return TRUE;
}

void ballnode_free(BallNode *node) {
    if (!ballnode_is_leaf(node)) {
        ballnode_free(node->left);
        ballnode_free(node->right);
    }
    free(node);
}

void attach_childs(BallNode *node, PointSlice *points, int leafsize) {
    enum Axis split_axis = pointslice_get_maxspread_axis(points);
    int split_idx = pointslice_median_partition(points, split_axis);
    if (split_idx == -1) {
        fprintf(stderr, "ERROR: could not determine the median element for the next split\n");
        return;
    }

    PointSlice child = {.points = points->points};
    child.start = node->data.start;
    child.end = split_idx;
    node->left = ballnode_build_recursive(&child, leafsize);
    if (!node->left) {
        return;
    }

    child.start = split_idx;
    child.end = node->data.end;
    node->right = ballnode_build_recursive(&child, leafsize);
    if (!node->right) {
        ballnode_free(node->left);
        node->left = NULL;
    }
}

double sum_weights(const PointSlice *slice) {
    double sumw = 0.0;
    Point *points = slice->points;
    for (size_t i = slice->start; i < slice->end; ++i) {
        Point *point_i = points + i;
        sumw += point_i->weight;
    }
    return sumw;
}

BallNode *ballnode_build_recursive(PointSlice *slice, int leafsize) {
    int size = pointslice_get_size(slice);
    if (size < 1) {
        return NULL;
    }
    Point center = pointslice_get_mean(slice);
    double radius = pointslice_get_maxdist(slice, &center);

    BallNode *node = ballnode_create_node(center, radius, slice);
    // TODO: check null pointer
    if (size <= leafsize) {
        node->sum_weight = sum_weights(&node->data);
        return node;
    }

    if (node) {
        attach_childs(node, slice, leafsize);
    }
    // TODO: use is_leaf()?
    if (!node->left || !node->right) {
        return NULL;
    }
    node->sum_weight = node->left->sum_weight + node->right->sum_weight;
    return node;
}

int ballnode_count_nodes(const BallNode *node) {
    int count = 1;
    if (node->left) {
        count += ballnode_count_nodes(node->left);
    }
    if (node->right) {
        count += ballnode_count_nodes(node->right);
    }
    return count;
}

double sum_weights_within_radius2(const PointSlice *slice, const Point *point, double radius2) {
    double sumw = 0.0;
    Point *points = slice->points;
    for (size_t i = slice->start; i < slice->end; ++i) {
        Point *point_i = points + i;
        double distance2 = point_dist_squared(point_i, point);
        int distance_mask = distance2 <= radius2;
        sumw += point_i->weight * (double)distance_mask;  // branchless
    }
    return sumw;
}

double ballnode_count_radius(const BallNode *node, const Point *point, double radius) {
    double distance = point_dist(&node->center, point);
    if (distance <= radius - node->radius) {
        // all points are pairs, stop distance evaluation
        return point->weight * node->sum_weight;
    } else if (distance <= radius + node->radius) {
        // some points may be pairs
        if (ballnode_is_leaf(node)) {
            // check each point individually
            return point->weight * sum_weights_within_radius2(&node->data, point, radius * radius);
        }
        // traverse the tree to narrow down the node size
        double counts = 0.0;
        counts += ballnode_count_radius(node->left, point, radius);
        counts += ballnode_count_radius(node->right, point, radius);
        return counts;
    }
    return 0.0;
}

double sum_weights_within_range2(const PointSlice *slice, const Point *point, double rmin2, double rmax2) {
    double sumw = 0.0;
    Point *points = slice->points;
    for (size_t i = slice->start; i < slice->end; ++i) {
        Point *point_i = points + i;
        double distance2 = point_dist_squared(point_i, point);
        int distance_mask = rmin2 < distance2 || distance2 <= rmax2;
        sumw += point_i->weight * (double)distance_mask; // branchless
    }
    return sumw;
}

double ballnode_count_range(const BallNode *node, const Point *point, double rmin, double rmax) {
    double distance = point_dist(&node->center, point);
    if (rmin + node->radius < distance && distance <= rmax - node->radius) {
        // all points are pairs, stop distance evaluation
        return point->weight * node->sum_weight;
    } else if (rmin - node->radius < distance || distance <= rmax + node->radius) {
        // some points may be pairs
        if (ballnode_is_leaf(node)) {
            // check each point individually
            return point->weight * sum_weights_within_range2(&node->data, point, rmin * rmin, rmax * rmax);
        }
        // traverse the tree to narrow down the node size
        double counts = 0.0;
        counts += ballnode_count_range(node->left, point, rmin, rmax);
        counts += ballnode_count_range(node->right, point, rmin, rmax);
        return counts;
    }
    return 0.0;
}

double dualsum_weights_within_radius2(const PointSlice *slice1, const PointSlice *slice2, double radius2) {
    double sumw = 0.0;
    Point *points1 = slice1->points;
    for (size_t i = slice1->start; i < slice1->end; ++i) {
        Point *point1_i = points1 + i;
        sumw += point1_i->weight * sum_weights_within_radius2(slice2, point1_i, radius2);
    }
    return sumw;
}

double ballnode_dualcount_radius(const BallNode *node1, const BallNode *node2, double radius) {
    double distance = point_dist(&node1->center, &node2->center);
    if (distance <= radius - node1->radius - node2->radius) {
        // all points are pairs, stop distance evaluation
        return node1->sum_weight * node2->sum_weight;
    } else if (distance <= radius + node1->radius + node2->radius) {
        // some points may be pairs
        int node1_is_leaf = ballnode_is_leaf(node1);
        int node2_is_leaf = ballnode_is_leaf(node2);
        double counts = 0.0;
        if (node1_is_leaf && node2_is_leaf) {
            // check all pairs of points individually
            return dualsum_weights_within_radius2(&node1->data, &node2->data, radius * radius);
        }
        // traverse the tree to narrow down the node size
        if (node1_is_leaf) {
            counts += ballnode_dualcount_radius(node1, node2->left, radius);
            counts += ballnode_dualcount_radius(node1, node2->right, radius);
        } else if (node2_is_leaf) {
            counts += ballnode_dualcount_radius(node1->left, node2, radius);
            counts += ballnode_dualcount_radius(node1->right, node2, radius);
        } else {
            counts += ballnode_dualcount_radius(node1->left, node2->left, radius);
            counts += ballnode_dualcount_radius(node1->left, node2->right, radius);
            counts += ballnode_dualcount_radius(node1->right, node2->left, radius);
            counts += ballnode_dualcount_radius(node1->right, node2->right, radius);
        }
        return counts;
    }
    return 0.0;
}

int ballnode_serialise_recursive(BallNodeBuffer buffer, BallNode *node, int insertion_index) {
    if (*buffer.next_free > buffer.size) {
        return FAILED;
    }
    int is_leaf = ballnode_is_leaf(node);

    int index_left, index_right;
    if (is_leaf) {
        index_left = -1;
        index_right = -1;
    } else {
        index_left = (*buffer.next_free)++;
        index_right = (*buffer.next_free)++;
    }
    BallNodeSerialized serialized = {
        .center_x = node->center.x,
        .center_y = node->center.y,
        .center_z = node->center.z,
        .radius = node->radius,
        .sum_weight = node->sum_weight,
        .left = index_left,
        .right = index_right,
        .data_start = node->data.start,
        .data_end = node->data.end,
    };
    buffer.buffer[insertion_index] = serialized;

    if (!is_leaf) {
        if (!ballnode_serialise_recursive(buffer, node->left, serialized.left)) {
            return FAILED;
        }
        if (!ballnode_serialise_recursive(buffer, node->right, serialized.right)) {
            return FAILED;
        }
    }
    return SUCCESS;
}

BallNode *ballnode_deserialise_recursive(BallNodeSerialized *buffer, int buffer_size, const PointBuffer *points, int index) {
    if (index >= buffer_size) {
        return NULL;
    }
    BallNode *node = (BallNode *)calloc(sizeof(BallNode), 1);
    if (!node) {
        return NULL;
    }
    BallNodeSerialized *serialized = buffer + index;

    if (serialized->data_end > points->size) {
        fprintf(stderr, "ERROR: point buffer does not match data slice expected by node\n");
        free(node);
        return NULL;
    }
    node->center = point_create(serialized->center_x, serialized->center_y, serialized->center_z);
    node->radius = serialized->radius;
    node->sum_weight = serialized->sum_weight;
    node->data = (PointSlice){
        .start = serialized->data_start,
        .end = serialized->data_end,
        .points = points->points,
    };

    // recurse into potential leaf nodes
    if (serialized->left != -1) {
        node->left = ballnode_deserialise_recursive(buffer, buffer_size, points, serialized->left);
        if (!node->left) {
            free(node);
            return NULL;
        }
    }
    if (serialized->right != -1) {
        node->right = ballnode_deserialise_recursive(buffer, buffer_size, points, serialized->right);
        if (!node->right) {
            free(node);
            return NULL;
        }
    }
    return node;
}
