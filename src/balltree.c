#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "point.h"
#include "balltree.h"

#define TRUE  1
#define FALSE 0

struct BallNode* ballnode_build_recursive(struct PointSlice*, int);

struct BallNode* ballnode_create_node(struct Point center, double radius, const struct PointSlice *slice)
{
    struct BallNode *node = (struct BallNode*)calloc(1, sizeof(struct BallNode));
    if (!node) {
        fprintf(stderr, "ERROR: BallNode node allocation failed");
        return NULL;
    }
    node->center = center;
    node->radius = radius;
    node->data = (struct PointSlice){
        .start = slice->start,
        .end = slice->end,
        .points = slice->points,
    };
    return node;
}

int ballnode_is_leaf(const struct BallNode *node)
{
    if (node->left && node->right) {
        return FALSE;
    }
    return TRUE;
}

void ballnode_free(struct BallNode *node)
{
    if (!ballnode_is_leaf(node)) {
        if (node->left) {
            ballnode_free(node->left);
        }
        if (node->right) {
            ballnode_free(node->right);
        }
    }
    free(node);
}

struct BallNode* create_child(struct PointSlice *parent, int leafsize, int split_index, int left)
{
    struct PointSlice child;
    child.points = parent->points;
    if (left) {
        child.start = parent->start;
        child.end = split_index;
    } else {
        child.start = split_index;
        child.end = parent->end;
    }
    return ballnode_build_recursive(&child, leafsize);
}

void attach_childs(struct BallNode *node, struct PointSlice *points, int leafsize)
{
    enum Axis split_axis = get_max_spread_axis(points);
    int i_split = partial_median_sort(points, split_axis);
    if (i_split == -1) {
        fprintf(stderr, "ERROR: could not determine the median element for the next split");
        print_pointslice(points);
        return;
    }

    node->left = create_child(points, leafsize, i_split, TRUE);
    if (!node->left) {
        return;
    }
    node->right = create_child(points, leafsize, i_split, FALSE);
    if (!node->right) {
        ballnode_free(node->left);
        node->left = NULL;
    }
}

struct BallNode* ballnode_build_recursive(struct PointSlice *slice, int leafsize)
{
    int size = get_pointslice_size(slice);
    if (size < 1) {
        return NULL;
    }
    struct Point center = get_center_point(slice);
    double radius = get_maxdist_from_center(slice, &center);

    struct BallNode *node = ballnode_create_node(center, radius, slice);
    if (size <= leafsize) {
        return node;
    }

    if (node) {
        attach_childs(node, slice, leafsize);
    }
    if (!node->left || !node->right) {
        return NULL;
    }
    return node;
}

double sum_weights(struct PointSlice *slice)
{
    double sumw = 0.0;
    struct Point *points = slice->points;
    for (size_t i = slice->start; i < slice->end; ++i) {
        struct Point *point_i = points + i;
        sumw += point_i->weight;
    }
    return sumw;
}

double sum_weights_within_radius2(struct PointSlice *slice, const struct Point *point, double radius2)
{
    double sumw = 0.0;
    struct Point *points = slice->points;
    for (size_t i = slice->start; i < slice->end; ++i) {
        struct Point *point_i = points + i;
        double distance2 = points_distance2(point_i, point);
        if (distance2 <= radius2) {
            sumw += point_i->weight;
        }
    }
    return sumw;
}

double ballnode_count_radius(struct BallNode *node, struct Point *point, double radius)
{
    double distance = points_distance(&node->center, point);
    double node_radius = node->radius;
    if (distance <= radius - node_radius) {
        // all points are pairs, stop distance evaluation
        return point->weight * sum_weights(&node->data);;
    } else if (distance <= radius + node_radius) {
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

double sum_weights_within_range2(struct PointSlice *slice, const struct Point *point, double rmin2, double rmax2)
{
    double sumw = 0.0;
    struct Point *points = slice->points;
    for (size_t i = slice->start; i < slice->end; ++i) {
        struct Point *point_i = points + i;
        double distance2 = points_distance2(point_i, point);
        if (rmin2 < distance2 || distance2 <= rmax2) {
            sumw += point_i->weight;
        }
    }
    return sumw;
}

double ballnode_count_range(struct BallNode *node, struct Point *point, double rmin, double rmax)
{
    double distance = points_distance(&node->center, point);
    double node_radius = node->radius;
    if (rmin + node_radius < distance && distance <= rmax - node_radius) {
        // all points are pairs, stop distance evaluation
        return point->weight * sum_weights(&node->data);;
    } else if (rmin - node_radius < distance || distance <= rmax + node_radius) {
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

double dualsum_weights_within_radius2(struct PointSlice *slice1, struct PointSlice *slice2, double radius2)
{
    double sumw = 0.0;
    struct Point *points1 = slice1->points;
    for (size_t i = slice1->start; i < slice1->end; ++i) {
        struct Point *point1_i = points1 + i;
        sumw += point1_i->weight * sum_weights_within_radius2(slice2, point1_i, radius2);
    }
    return sumw;
}

double ballnode_dualcount_radius(struct BallNode *node1, struct BallNode *node2, double radius)
{
    double distance = points_distance(&node1->center, &node2->center);
    double node1_radius = node1->radius;
    double node2_radius = node2->radius;
    if (distance <= radius - node1_radius - node2_radius) {
        // all points are pairs, stop distance evaluation
        return sum_weights(&node1->data) * sum_weights(&node2->data);;
    } else if (distance <= radius + node1_radius + node2_radius) {
        // some points may be pairs
        int node1_is_leaf = ballnode_is_leaf(node1);
        int node2_is_leaf = ballnode_is_leaf(node2);
        double counts = 0.0;
        if (node1_is_leaf && node2_is_leaf) {
            // check all pairs of points individually
            return dualsum_weights_within_radius2(&node1->data, &node2->data, radius * radius);
        }
        // traverse the tree to narrow down the node size
        if (!node1_is_leaf) {  // node2 is leaf
            counts += ballnode_dualcount_radius(node1->left, node2, radius);
            counts += ballnode_dualcount_radius(node1->right, node2, radius);
        } else if (!node2_is_leaf) {  // node1 is leaf
            counts += ballnode_dualcount_radius(node1, node2->left, radius);
            counts += ballnode_dualcount_radius(node1, node2->right, radius);
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

// public interface

struct BallTree* balltree_build(struct PointBuffer *buffer, int leafsize)
{
    struct BallTree *tree = (struct BallTree*)malloc(sizeof(struct BallTree));
    if (!tree) {
        return NULL;
    }
    tree->leafsize = leafsize;

    // copy the input data buffer, since it will be partitioned while build the tree
    struct PointBuffer data = {buffer->size, NULL};
    size_t n_bytes = buffer->size * sizeof(struct Point);
    data.points = (struct Point*)malloc(n_bytes);
    if (!data.points) {
        free(tree);
        return NULL;
    }
    memcpy(data.points, buffer->points, n_bytes);
    tree->data = data;

    struct PointSlice *slice = pointslice_from_buffer(&data);
    tree->root = ballnode_build_recursive(slice, leafsize);
    if (!tree->root) {
        free(data.points);
        free(tree);
        return NULL;
    }
    return tree;
}

void balltree_free(struct BallTree *tree)
{
    if (tree->data.points) {
        free(tree->data.points);
    }
    if (tree->root) {
        free(tree->root);
    }
    free(tree);
    tree = NULL;
}

double balltree_count_radius(struct BallTree *tree, struct Point *point, double radius)
{
    return ballnode_count_radius(tree->root, point, radius);
}

double balltree_count_range(struct BallTree *tree, struct Point *point, double rmin, double rmax)
{
    return ballnode_count_range(tree->root, point, rmin, rmax);
}

double balltree_dualcount_radius(struct BallTree *tree1, struct BallTree *tree2, double radius)
{
    return ballnode_dualcount_radius(tree1->root, tree2->root, radius);
}
