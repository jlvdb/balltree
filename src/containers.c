#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "containers.h"
#include "balltree_macros.h"

DistHistogram *hist_new(int64_t size, double *dist_edges) {
    if (size < 1) {
        EMIT_ERR_MSG(ValueError, "DistHistogram requires at least 1 edges");
        return NULL;
    }
    size_t n_bytes_double = size * sizeof(double);

    DistHistogram *hist = malloc(sizeof(DistHistogram));
    if (hist == NULL) {
        EMIT_ERR_MSG(MemoryError, "DistHistogram allocation failed");
        return NULL;
    }
    hist->size = size;

    // allocate and initialise hist.sum_weight buffer
    double *sum_weight = calloc(size, sizeof(double));
    if (sum_weight == NULL) {
        EMIT_ERR_MSG(MemoryError, "DistHistogram.sum_weight allocation failed");
        hist_free(hist);
        return NULL;
    }
    hist->sum_weight = sum_weight;

    // allocate hist.dist buffer
    double *dist = malloc(n_bytes_double);
    if (dist == NULL) {
        EMIT_ERR_MSG(MemoryError, "DistHistogram.dist allocation failed");
        hist_free(hist);
        return NULL;
    }
    hist->dist = dist;

    // allocate hist.dist_sq buffer
    double *dist_sq = malloc(n_bytes_double);
    if (dist_sq == NULL) {
        EMIT_ERR_MSG(MemoryError, "DistHistogram.dist_sq allocation failed");
        hist_free(hist);
        return NULL;
    }
    hist->dist_sq = dist_sq;

    // initialise dist/_sq buffers
    double edge = 0.0;
    for (int64_t i = 0; i < size; ++i) {
        edge = dist_edges[i];
        dist[i] = edge;
        dist_sq[i] = edge * edge;
    }
    hist->dist_max = edge;
    hist->dist_sq_max = edge * edge;
    return hist;
}

void hist_free(DistHistogram *hist) {
    if (hist->sum_weight != NULL) {
        free(hist->sum_weight);
    }
    if (hist->dist != NULL) {
        free(hist->dist);
    }
    if (hist->dist_sq != NULL) {
        free(hist->dist_sq);
    }
    free(hist);
}

int64_t hist_insert_dist_sq(DistHistogram *hist, double dist_sq, double weight) {
    int64_t bin_idx = -1;
    if (dist_sq <= hist->dist_sq_max) {
        for (bin_idx = 0; bin_idx <= hist->size; ++bin_idx) {
            if (dist_sq <= hist->dist_sq[bin_idx]) {
                hist->sum_weight[bin_idx] += weight;
                break;
            }
        }
    }
    return bin_idx;
}

KnnQueue *knque_new(int64_t capacity) {
    if (capacity < 1) {
        EMIT_ERR_MSG(ValueError, "KnnQueue capacity must be positive");
        return NULL;
    }

    KnnQueue *queue = malloc(sizeof(KnnQueue));
    if (queue == NULL) {
        EMIT_ERR_MSG(MemoryError, "KnnQueue allocation failed");
        return NULL;
    }

    QueueItem *items = malloc(capacity * sizeof(QueueItem));
    if (items == NULL) {
        EMIT_ERR_MSG(MemoryError, "KnnQueue.items allocation failed");
        knque_free(queue);
        return NULL;
    }

    queue->items = items;
    queue->capacity = capacity;
    knque_clear(queue);
    return queue;
}

void knque_free(KnnQueue *queue) {
    if (queue->items != NULL) {
        free(queue->items);
    }
    free(queue);
}

void knque_clear(KnnQueue *queue) {
    queue->size = 0;
    queue->distance_max = INFINITY;
    for (int64_t i = 0; i < queue->capacity; ++i) {
        queue->items[i] = (QueueItem){
            .index = -1,
            .distance = INFINITY,
        };
    }
}

int knque_insert(KnnQueue *queue, int64_t item_index, double distance) {
    QueueItem *items = queue->items;
    if (distance >= knque_get_max_dist(queue)) {
        return 1;  // item not in queue
    }

    // find insertion index, note that (distance < distance_last_element)
    int64_t idx_insert = queue->size;
    while (idx_insert > 0 && distance < items[idx_insert - 1].distance) {
        --idx_insert;
    }
    // very first if statement guarantees (idx_insert < queue->capacity)

    // make room and insert item, drop last item if at capacity
    int queue_is_full = knque_is_full(queue);
    int64_t idx = queue_is_full ? (queue->capacity - 1) : queue->size;
    while (idx > idx_insert) {
        items[idx] = items[idx - 1];
        --idx;   
    }
    items[idx_insert].index = item_index;
    items[idx_insert].distance = distance;

    // update state of queue
    if (!queue_is_full) {
        ++(queue->size);
    }
    return 0;
}
