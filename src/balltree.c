#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "point.h"
#include "balltree.h"

#define TRUE  1
#define FALSE 0

struct BallTree* balltree_create_node(struct Point center, double radius, const struct PointSlice *slice)
{
    struct BallTree *node = (struct BallTree*)calloc(1, sizeof(struct BallTree));
    if (!node) {
        fprintf(stderr, "ERROR: BallTree node allocation failed");
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

int balltree_is_leaf(const struct BallTree *node)
{
    if (node->left && node->right) {
        return FALSE;
    }
    return TRUE;
}

void balltree_free(struct BallTree *node)
{
    if (!balltree_is_leaf(node)) {
        if (node->left) {
            balltree_free(node->left);
        }
        if (node->right) {
            balltree_free(node->right);
        }
    }
    free(node);
}

void balltree_print_indented(const struct BallTree *node, int indent)
{
    int padding = 2;
    if (balltree_is_leaf(node)) {
        printf(
            "%*sLeaf(size=%d)\n",
            indent * padding,
            "",
            get_pointslice_size(&node->data)
        );
    } else {
        printf(
            "%*sNode(radius=%.3f)\n",
            indent * padding,
            "",
            node->radius
        );
        ++indent;
        balltree_print_indented(node->left, indent);
        balltree_print_indented(node->right, indent);
    }
}

void balltree_print(const struct BallTree *node)
{
    balltree_print_indented(node, 0);
}

struct BallTree* create_child(struct PointSlice *parent, int leafsize, int split_index, int left)
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
    return balltree_build_recursive(&child, leafsize);
}

void attach_childs(struct BallTree *node, struct PointSlice *points, int leafsize)
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
        balltree_free(node->left);
        node->left = NULL;
    }
}

struct BallTree* balltree_build_recursive(struct PointSlice *slice, int leafsize)
{
    int size = get_pointslice_size(slice);
    if (size < 1) {
        return NULL;
    }
    struct Point center = get_center_point(slice);
    double radius = get_maxdist_from_center(slice, &center);

    struct BallTree *node = balltree_create_node(center, radius, slice);
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

struct BallTree* balltree_build(struct PointBuffer *buffer, int leafsize)
{
    struct PointSlice *slice = pointslice_from_buffer(buffer);
    struct BallTree *tree = balltree_build_recursive(slice, leafsize);
    free(slice);
    return tree;
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

double balltree_count_radius(struct BallTree *node, struct Point *point, double radius)
{
    double distance = points_distance(&node->center, point);
    double node_radius = node->radius;
    if (distance <= radius - node_radius) {
        // all points are pairs, stop distance evaluation
        return point->weight * sum_weights(&node->data);;
    } else if (distance <= radius + node_radius) {
        // some points may be pairs
        if (balltree_is_leaf(node)) {
            // check each point individually
            return point->weight * sum_weights_within_radius2(&node->data, point, radius * radius);
        }
        // traverse the tree to narrow down the node size
        double counts = 0.0;
        counts += balltree_count_radius(node->left, point, radius);
        counts += balltree_count_radius(node->right, point, radius);
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

double balltree_count_range(struct BallTree *node, struct Point *point, double rmin, double rmax)
{
    double distance = points_distance(&node->center, point);
    double node_radius = node->radius;
    if (rmin + node_radius < distance && distance <= rmax - node_radius) {
        // all points are pairs, stop distance evaluation
        return point->weight * sum_weights(&node->data);;
    } else if (rmin - node_radius < distance || distance <= rmax + node_radius) {
        // some points may be pairs
        if (balltree_is_leaf(node)) {
            // check each point individually
            return point->weight * sum_weights_within_range2(&node->data, point, rmin * rmin, rmax * rmax);
        }
        // traverse the tree to narrow down the node size
        double counts = 0.0;
        counts += balltree_count_range(node->left, point, rmin, rmax);
        counts += balltree_count_range(node->right, point, rmin, rmax);
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

double balltree_dualcount_radius(struct BallTree *node1, struct BallTree *node2, double radius)
{
    double distance = points_distance(&node1->center, &node2->center);
    double node1_radius = node1->radius;
    double node2_radius = node2->radius;
    if (distance <= radius - node1_radius - node2_radius) {
        // all points are pairs, stop distance evaluation
        return sum_weights(&node1->data) * sum_weights(&node2->data);;
    } else if (distance <= radius + node1_radius + node2_radius) {
        // some points may be pairs
        int node1_is_leaf = balltree_is_leaf(node1);
        int node2_is_leaf = balltree_is_leaf(node2);
        double counts = 0.0;
        if (node1_is_leaf && node2_is_leaf) {
            // check all pairs of points individually
            return dualsum_weights_within_radius2(&node1->data, &node2->data, radius * radius);
        }
        // traverse the tree to narrow down the node size
        if (!node1_is_leaf) {  // node2 is leaf
            counts += balltree_dualcount_radius(node1->left, node2, radius);
            counts += balltree_dualcount_radius(node1->right, node2, radius);
        } else if (!node2_is_leaf) {  // node1 is leaf
            counts += balltree_dualcount_radius(node1, node2->left, radius);
            counts += balltree_dualcount_radius(node1, node2->right, radius);
        } else {
            counts += balltree_dualcount_radius(node1->left, node2->left, radius);
            counts += balltree_dualcount_radius(node1->left, node2->right, radius);
            counts += balltree_dualcount_radius(node1->right, node2->left, radius);
            counts += balltree_dualcount_radius(node1->right, node2->right, radius);
        }
        return counts;
    }
    return 0.0;
}
