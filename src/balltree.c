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
        perror("BallTree node allocation failed");
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
        perror("BallTree leaf data allocation failed");
        return NULL;
    }

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
    return balltree_build(&child, leafsize);
}

void attach_childs(struct BallTree *node, struct PointSlice *points, int leafsize)
{
    enum Axis split_axis = get_max_spread_axis(points);
    int i_split = partial_median_sort(points, split_axis);
    if (i_split == -1) {
        perror("could not determine the median element for the next split");
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
    }
}

struct BallTree* balltree_build(struct PointSlice *points, int leafsize)
{
    int size = get_pointslice_size(points);
    if (size < 1) {
        return NULL;
    }
    struct Point center = get_center_point(points);
    double radius = get_maxdist_from_center(points, center);

    if (size <= leafsize) {
        return balltree_create_leaf(center, radius, points);
    }

    struct BallTree *node = balltree_create_node(center, radius);
    if (node) {
        attach_childs(node, points, leafsize);
    }
    if (!node->left || !node->right) {
        return NULL;
    }
    return node;
}

void balltree_free(struct BallTree *node)
{
    if (balltree_is_leaf(node)) {
        free(node->data.points);
        node->data.points = NULL;
    } else {
        if (node->left) {
            balltree_free(node->left);
        }
        if (node->right) {
            balltree_free(node->right);
        }
    }
    free(node);
    node = NULL;
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

double count_within_radius(struct PointBuffer *buffer, struct Point *point, double radius) {
    double radius2 = radius * radius;
    double counts = 0.0;

    struct Point *points = buffer->points;
    for (size_t i = 0; i < buffer->size; ++i) {
        double distance2 = points_distance2(points + i, point);
        if (distance2 <= radius2) {
            counts += 1.0;  // we want weights later
        }
    }
    return counts;
}

double count_within_range(struct PointBuffer *buffer, struct Point *point, double rmin, double rmax) {
    double rmin2 = rmin * rmin;
    double rmax2 = rmax * rmax;
    double counts = 0.0;

    struct Point *points = buffer->points;
    for (size_t i = 0; i < buffer->size; ++i) {
        double distance2 = points_distance2(points + i, point);
        if (rmin2 < distance2 && distance2 <= rmax2) {
            counts += 1.0;  // we want weights later
        }
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
