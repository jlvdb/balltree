#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "queue.h"
#include "balltree_macros.h"

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
    int queue_is_full = queue->size == queue->capacity;
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
