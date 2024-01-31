#include <stdint.h>
#include <stdlib.h>

#include "ballnode.h"
#include "balltree_macros.h"

StatsVector *statvec_new() {
    return statvec_new_reserve(32L);
}

StatsVector *statvec_new_reserve(int64_t reserve_capacity) {
    if (reserve_capacity < 1) {
        EMIT_ERR_MSG(MemoryError, "StatsVector capacity must be positive");
        return NULL;
    }

    StatsVector *vec = malloc(sizeof(StatsVector));
    if (vec == NULL) {
        EMIT_ERR_MSG(MemoryError, "StatsVector allocation failed");
        return NULL;
    }

    vec->stats = malloc(reserve_capacity * sizeof(NodeStats));
    if (vec->stats == NULL) {
        EMIT_ERR_MSG(MemoryError, "StatsVector.stats allocation failed");
        statvec_free(vec);
        return NULL;
    }
    vec->capacity = reserve_capacity;
    vec->size = 0;
    return vec;
}

void statvec_free(StatsVector *vec) {
    if (vec->stats) {
        free(vec->stats);
    }
    free(vec);
}

int statvec_resize(StatsVector *vec, int64_t capacity) {
    if (capacity < 1) {
        EMIT_ERR_MSG(ValueError, "StatsVector capacity must be positive");
        return BTR_FAILED;
    }

    size_t n_bytes = capacity * sizeof(NodeStats);
    NodeStats *stats = realloc(vec->stats, n_bytes);
    if (stats == NULL) {
        EMIT_ERR_MSG(MemoryError, "StatsVector resizing failed");
        return BTR_FAILED;
    }

    vec->stats = stats;
    vec->capacity = capacity;
    vec->size = (vec->size > capacity) ? capacity : vec->size;
    return BTR_SUCCESS;
}

int statvec_append(StatsVector *vec, const NodeStats *stats) {
    if (vec->size >= vec->capacity) {
        // double the vector size if necessary
        if (statvec_resize(vec, vec->capacity * 2) == BTR_FAILED) {
            return BTR_FAILED;
        }
    }
    vec->stats[vec->size] = *stats;
    ++(vec->size);
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

int64_t bnode_count_nodes(const BallNode *node) {
    int64_t count = 1;
    if (!BALLNODE_IS_LEAF(node)) {
        count += bnode_count_nodes(node->childs.left);
        count += bnode_count_nodes(node->childs.right);
    }
    return count;
}
