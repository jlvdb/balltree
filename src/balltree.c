#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "point.h"
#include "histogram.h"
#include "ballnode.h"
#include "balltree.h"
#include "balltree_macros.h"

BallTree *balltree_build(const PointBuffer *buffer) {
    return balltree_build_leafsize(buffer, DEFAULT_LEAFSIZE);
}

BallTree *balltree_build_leafsize(const PointBuffer *buffer, int leafsize) {
    PointBuffer *data = ptbuf_copy(buffer);
    if (data == NULL) {
        return NULL;
    }

    BallTree *tree = balltree_build_nocopy(data, leafsize);
    if (tree == NULL) {
        ptbuf_free(data);
        return NULL;
    }
    tree->data_owned = 1;
    return tree;
}

BallTree *balltree_build_nocopy(PointBuffer *buffer, int leafsize) {
    int64_t size = buffer->size;
    if (size < 1) {
        EMIT_ERR_MSG(ValueError, "need at least one input data point to build a tree");
        return NULL;
    }

    PointSlice *slice = ptslc_from_buffer(buffer);
    if (slice == NULL) {
        return NULL;
    }
    BallNode *root = bnode_build(slice, leafsize);
    free(slice);
    if (root == NULL) {
        return NULL;
    }

    BallTree *tree = malloc(sizeof(BallTree));
    if (tree == NULL) {
        EMIT_ERR_MSG(MemoryError, "BallTree root allocation failed");
        bnode_free(root);
        return NULL;
    }
    tree->leafsize = leafsize;
    tree->data_owned = 0;
    tree->data = buffer;
    tree->root = root;
    return tree;
}

void balltree_free(BallTree *tree) {
    if (tree->data_owned && tree->data != NULL) {
        ptbuf_free(tree->data);
    }
    if (tree->root != NULL) {
        bnode_free(tree->root);
    }
    free(tree);
}

int64_t balltree_count_nodes(const BallTree *tree) {
    return bnode_count_nodes(tree->root);
}

double balltree_brute_radius(
    const BallTree *tree,
    const Point *point,
    double radius
) {
    // avoid using *ptslc_from_buffer as the struct allocation could fail
    PointSlice slice = {
        .start = tree->data->points,
        .end = tree->data->points + tree->data->size,
    };
    return point->weight * ptslc_sumw_in_radius_sq(&slice, point, radius * radius);
}

void balltree_brute_range(
    const BallTree *tree,
    const Point *point,
    DistHistogram *hist
) {
    Point *points = tree->data->points;
    for (int64_t i = 0; i < tree->data->size; ++i) {
        double dist_sq = EUCLIDEAN_DIST_SQ(point, points + i);
        // increment corresponding bin by weight
        hist_insert_dist_sq(hist, dist_sq, point->weight * points[i].weight);
    }
}

double balltree_count_radius(
    const BallTree *tree,
    const Point *point,
    double radius
) {
    return bnode_count_radius(tree->root, point, radius);
}

void balltree_count_range(
    const BallTree *tree,
    const Point *point,
    DistHistogram *hist
) {
    bnode_count_range(tree->root, point, hist);
}

double balltree_dualcount_radius(
    const BallTree *tree1,
    const BallTree *tree2,
    double radius
) {
    return bnode_dualcount_radius(tree1->root, tree2->root, radius);
}

void balltree_dualcount_range(
    const BallTree *tree1,
    const BallTree *tree2,
    DistHistogram *hist
) {
    bnode_dualcount_range(tree1->root, tree2->root, hist);
}

StatsVector *balltree_collect_stats(const BallTree *tree) {
    int64_t num_nodes = balltree_count_nodes(tree);
    StatsVector *vec = statvec_new_reserve(num_nodes);
    if (vec == NULL) {
        return NULL;
    }

    if (bnode_collect_stats(tree->root, vec, 0) == BTR_FAILED) {
        statvec_free(vec);
        return NULL;
    }
    return vec;
}
