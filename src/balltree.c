#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "point.h"
#include "balltree.h"

#define TRUE  1
#define FALSE 0

struct BallTree* balltree_create_node(struct Point center, double radius)
{
    struct BallTree *node = (struct BallTree*)calloc(1, sizeof(struct BallTree));
    if (!node) {
        fprintf(stderr, "ERROR: BallTree node allocation failed");
        return NULL;
    }
    node->center = center;
    node->radius = radius;
    return node;
}

struct BallTree* balltree_create_leaf(struct Point center, double radius, const struct PointSlice *slice)
{
    int size = get_pointslice_size(slice);
    size_t n_bytes = size * sizeof(struct Point);

    struct BallTree *node = balltree_create_node(center, radius);
    if (!node) {
        return NULL;
    }
    node->data.size = size;
    node->data.points = (struct Point*)malloc(n_bytes);
    if (!node->data.points) {
        fprintf(stderr, "ERROR: BallTree leaf data allocation failed");
        free(node);
        return NULL;
    }

    const struct Point* slice_offset = slice->points + slice->start;
    memcpy(node->data.points, slice_offset, n_bytes);
    return node;
}

int balltree_is_leaf(const struct BallTree *node)
{
    if (node->data.points == NULL) {
        return FALSE;
    } else {
        return TRUE;
    }
}

void balltree_free(struct BallTree *node)
{
    if (balltree_is_leaf(node)) {
        free(node->data.points);
    } else {
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
            node->data.size
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

    if (size <= leafsize) {
        return balltree_create_leaf(center, radius, slice);
    }

    struct BallTree *node = balltree_create_node(center, radius);
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

double bulk_count(struct BallTree *node) {
    double counts = 0.0;
    if (balltree_is_leaf(node)) {
        counts += (double)node->data.size;
    } else {
        counts += bulk_count(node->left);
        counts += bulk_count(node->right);
    }
    return counts;
}

double balltree_count_radius(struct BallTree *node, struct Point *point, double radius) {
    double counts = 0.0;
    double distance = points_distance(&node->center, point);
    double node_radius = node->radius;

    if (balltree_is_leaf(node)) {
        if (distance <= radius - node_radius) {
            counts += bulk_count(node);  // we want weights later
        } else if (distance <= radius + node_radius) {
            counts += count_within_radius(&node->data, point, radius);
        }
    } else {
        counts += balltree_count_radius(node->left, point, radius);
        counts += balltree_count_radius(node->right, point, radius);
    }
    return counts;
}

double balltree_count_range(struct BallTree *node, struct Point *point, double rmin, double rmax) {
    double counts = 0.0;
    double distance = points_distance(&node->center, point);
    double node_radius = node->radius;

    if (balltree_is_leaf(node)) {
        if (rmin + node_radius < distance && distance <= rmax - node_radius) {
            counts += bulk_count(node);  // we want weights later
        } else if (rmin - node_radius < distance || distance <= rmax + node_radius) {
            counts += count_within_range(&node->data, point, rmin, rmax);
        }
    } else {
        // traverse
        counts += balltree_count_range(node->left, point, rmin, rmax);
        counts += balltree_count_range(node->right, point, rmin, rmax);
    }
    return counts;
}
