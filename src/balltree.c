#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "point.h"
#include "balltree.h"

#define TRUE  1
#define FALSE 0

#define eprintf(str) fprintf (stderr, "ERROR: %s\n", str)

struct BallTree *balltree_create_node(struct Point center, double radius)
{
    struct BallTree *node = (struct BallTree*)calloc(1, sizeof(struct BallTree));
    if (!node) {
        eprintf("BallTree node allocation failed");
        return NULL;
    }
    node->center = center;
    node->radius = radius;
    return node;
}

struct BallTree *balltree_create_leaf(struct Point center, double radius, const struct PointSlice *slice)
{
    size_t size = get_pointslice_size(slice);
    size_t n_bytes = size * sizeof(struct Point);

    struct BallTree *node = balltree_create_node(center, radius);
    if (!node) {
        return NULL;
    }
    node->data.size = size;
    node->data.points = (struct Point*)malloc(n_bytes);
    if (!node->data.points) {
        eprintf("BallTree leaf data allocation failed");
        return NULL;
    }
    memcpy(node->data.points, slice->points, n_bytes);

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
            "%*sLeaf(size=%zu)\n",
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
        indent++;
        balltree_print_indented(node->left, indent);
        balltree_print_indented(node->right, indent);
    }
}

void balltree_print(const struct BallTree *node)
{
    balltree_print_indented(node, 0);
}

struct BallTree *create_child(struct PointSlice *parent, int leafsize, size_t split_index, int left)
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
        eprintf("could not determine the median element for the next split");
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

struct BallTree *balltree_build(struct PointSlice *points, int leafsize)
{
    size_t size = get_pointslice_size(points);
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
