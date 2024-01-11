#include <stdlib.h>

#include "ballnode.h"
#include "balltree_macros.h"

StatsVector *statvec_new() {
    return statvec_new_reserve(32);
}

StatsVector *statvec_new_reserve(int reserve_size) {
    StatsVector *vec = (StatsVector *)malloc(sizeof(StatsVector));
    if (vec == NULL) {
        EMIT_ERR_MSG(MemoryError, "StatsVector allocation failed");
        return NULL;
    }
    vec->size = reserve_size;
    vec->end = 0;

    vec->stats = (NodeStats *)malloc(vec->size * sizeof(NodeStats));
    if (vec->stats == NULL) {
        EMIT_ERR_MSG(MemoryError, "StatsVector buffer allocation failed");
        free(vec);
        return NULL;
    }
    return vec;
}

void statvec_free(StatsVector *vec) {
    if (vec->stats) {
        free(vec->stats);
    }
    free(vec);
}

int statvec_resize(StatsVector *vec, int size) {
    if (size < 1) {
        EMIT_ERR_MSG(ValueError, "StatsVector size must be positive");
        return BTR_FAILED;
    }

    size_t n_bytes = size * sizeof(NodeStats);
    NodeStats *stats = (NodeStats *)realloc(vec->stats, n_bytes);
    if (stats == NULL) {
        EMIT_ERR_MSG(MemoryError, "StatsVector resizing failed");
        return BTR_FAILED;
    }

    vec->stats = stats;
    vec->size = size;
    vec->end = (vec->end > size) ? size : vec->end;
    return BTR_SUCCESS;
}

int statvec_append(StatsVector *vec, const NodeStats *stats) {
    if (vec->end >= vec->size) {
        // double the vector size if necessary
        if (statvec_resize(vec, vec->size * 2) == BTR_FAILED) {
            return BTR_FAILED;
        }
    }
    vec->stats[vec->end] = *stats;
    ++(vec->end);
    return BTR_SUCCESS;
}

int bnode_collect_stats(const BallNode *node, StatsVector *vec, int depth) {
    NodeStats stats = {
        .depth = depth,
        .n_points = ptslc_get_size(&node->data),
        .sum_weight = node->sum_weight,
        .x = node->center.x,
        .y = node->center.y,
        .z = node->center.z,
        .radius = node->radius
    };
    if (statvec_append(vec, &stats) == BTR_FAILED) {
        return BTR_FAILED;
    }

    if (node->left != NULL) {
        if (bnode_collect_stats(node->left, vec, depth + 1) == BTR_FAILED) {
            return BTR_FAILED;
        }
    }
    if (node->right != NULL) {
        if (bnode_collect_stats(node->right, vec, depth + 1) == BTR_FAILED) {
            return BTR_FAILED;
        }
    }
    return BTR_SUCCESS;
}

int bnode_count_nodes(const BallNode *node) {
    int count = 1;
    if (node->left != NULL) {
        count += bnode_count_nodes(node->left);
    }
    if (node->right != NULL) {
        count += bnode_count_nodes(node->right);
    }
    return count;
}

