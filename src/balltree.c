#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "point.h"
#include "balltree.h"

#define SUCCESS 1
#define FAILED  0

#define TRUE  1
#define FALSE 0

struct BallTree *balltree_create_node(struct Point center, double radius)
{
    struct BallTree *node = (struct BallTree*)calloc(1, sizeof(struct BallTree));
    if (node == NULL)
        return NULL;
    node->center = center;
    node->radius = radius;
    return node;
}

struct BallTree *balltree_create_leaf(struct Point center, double radius, const struct PointSlice *slice)
{
    size_t size = get_pointslice_size(slice);
    size_t n_bytes = size * sizeof(struct Point);

    struct BallTree *node = balltree_create_node(center, radius);
    node->data.size = size;
    node->data.points = (struct Point*)malloc(n_bytes);
    if (node->data.points == NULL)
        return NULL;
    memcpy(node->data.points, slice->points, n_bytes);

    return node;
}

int balltree_is_leaf(const struct BallTree *node)
{
    if (node->data.points == NULL)
        return FALSE;
    else
        return TRUE;
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

int node_add_childs(struct BallTree *node, struct PointSlice *points, int leafsize)
{
    enum Axis split_axis = get_max_spread_axis(points);
    printf("splitting on %u\n", split_axis);
    int i_split = partial_median_sort(points, split_axis);
    printf("split done: start=%zu, split=%d, end=%zu\n", points->start, i_split, points->end);
    if (i_split == -1)
        return FAILED;

    struct PointSlice left_points = {
        .start = points->start,
        .end = i_split,
        .points = points->points,
    };
    node->left = balltree_build(&left_points, leafsize);
    if (!node->left)
        return FAILED;

    struct PointSlice right_points = {
        .start = i_split,
        .end = points->end,
        .points = points->points,
    };
    node->right = balltree_build(&right_points, leafsize);
    if (!node->right)
        return FAILED;
    return SUCCESS;
}

struct BallTree *balltree_build(struct PointSlice *points, int leafsize)
{
    size_t size = get_pointslice_size(points);
    if (size < 1)
        return NULL;
    struct Point center = get_center_point(points);
    double radius = get_maxdist_from_center(points, center);

    printf("center: "); print_point(&center);
    printf("radius: %.3f\n", radius);

    struct BallTree *node;
    if (size <= leafsize) {
        printf("building leaf [%zu,%zu)\n", points->start, points->end);
        node = balltree_create_leaf(center, radius, points);
    } else {
        printf("building node [%zu,%zu)\n", points->start, points->end);
        node = balltree_create_node(center, radius);
        if (node && !node_add_childs(node, points, leafsize)){
            node = NULL;
        }
    }
    return node;
}

void balltree_free(struct BallTree *node)
{
    if (node->data.points) {
        free(node->data.points);
        node->data.points = NULL;
    }
    if (node->left)
        balltree_free(node->left);
    if (node->right)
        balltree_free(node->right);
    free(node);
    node = NULL;
}
