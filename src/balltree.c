#include <stdio.h>
#include <stdlib.h>

#include "point.h"
#include "ballnode.h"
#include "balltree.h"
#include "balltree_macros.h"

#define DEFAULT_LEAFSIZE 20

BallTree* balltree_build(const PointBuffer *buffer) {
    return balltree_build_leafsize(buffer, DEFAULT_LEAFSIZE);
}

BallTree* balltree_build_leafsize(const PointBuffer *buffer, int leafsize) {
    int size = buffer->size;
    if (size < 1) {
        PRINT_ERR_MSG("need at least one input data point to build a tree\n");
        return NULL;
    }

    BallTree *tree = (BallTree*)malloc(sizeof(BallTree));
    if (tree == NULL) {
        PRINT_ERR_MSG("BallTree root allocation failed\n");
        return NULL;
    }
    tree->leafsize = leafsize;

    PointBuffer *data = ptbuf_copy(buffer);
    if (data == NULL) {
        balltree_free(tree);
        return NULL;
    }
    tree->data = *data;

    BallNode *root = bnode_build(data, 0, size, leafsize);
    if (root == NULL) {
        balltree_free(tree);
        return NULL;
    }
    tree->root = root;

    return tree;
}

void balltree_free(BallTree *tree) {
    if (tree->data.points != NULL) {
        free(tree->data.points);
    }
    if (tree->root != NULL) {
        bnode_free(tree->root);
    }
    free(tree);
}

int balltree_count_nodes(const BallTree *tree) {
    return bnode_count_nodes(tree->root);
}

double balltree_count_radius(
    const BallTree *tree,
    const Point *point,
    double radius
) {
    return bnode_count_radius(tree->root, point, radius);
}

double balltree_count_range(
    const BallTree *tree,
    const Point *point,
    double rmin,
    double rmax
) {
    return bnode_count_range(tree->root, point, rmin, rmax);
}

double balltree_dualcount_radius(
    const BallTree *tree1,
    const BallTree *tree2,
    double radius
) {
    return bnode_dualcount_radius(tree1->root, tree2->root, radius);
}

double balltree_dualcount_range(
    const BallTree *tree1,
    const BallTree *tree2,
    double rmin,
    double rmax
) {
    return bnode_dualcount_range(tree1->root, tree2->root, rmin, rmax);
}
