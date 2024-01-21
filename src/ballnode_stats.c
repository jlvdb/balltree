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
        .num_points = node->num_points,
        .sum_weight = node->sum_weight,
        .x = node->ball.x,
        .y = node->ball.y,
        .z = node->ball.z,
        .radius = node->ball.radius
    };
    if (statvec_append(vec, &stats) == BTR_FAILED) {
        return BTR_FAILED;
    }

    if (node->childs.left != NULL) {
        if (bnode_collect_stats(node->childs.left, vec, depth + 1) == BTR_FAILED) {
            return BTR_FAILED;
        }
    }
    if (node->childs.right != NULL) {
        if (bnode_collect_stats(node->childs.right, vec, depth + 1) == BTR_FAILED) {
            return BTR_FAILED;
        }
    }
    return BTR_SUCCESS;
}

int bnode_count_nodes(const BallNode *node) {
    int count = 1;
    if (!BALLNODE_IS_LEAF(node)) {
        count += bnode_count_nodes(node->childs.left);
        count += bnode_count_nodes(node->childs.right);
    }
    return count;
}

