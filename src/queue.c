#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "queue.h"
#include "balltree_macros.h"

static inline int knque_is_full(const KnnQueue *queue);


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

    QueueItem *items = calloc(capacity, sizeof(QueueItem));
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
    queue->dist_sq_max = INFINITY;
    for (long i = 0; i < queue->capacity; ++i) {
        queue->items[i] = (QueueItem){
            .value = QUEUEITEM_DEFAULT,
            .dist_sq = INFINITY,
        };
    }
}

static inline int knque_is_full(const KnnQueue *queue) {
    return queue->capacity == queue->size;
}

int knque_insert(KnnQueue *queue, const QUEUEITEM_T *value, double dist_sq) {
    if (dist_sq >= queue->dist_sq_max) {
        return 1;
    }
    // insertion index must now be at least at end of queue

    // the actual location is reached when the element at the preceeding index
    // is smaller than the item to be inserted
    QueueItem *items = queue->items;
    long idx_last = (queue->size > 0) ? (queue->size - 1) : 0;
    long idx_insert = idx_last;
    for (; idx_insert >= 0; --idx_insert) {
        double dist_sq_next = items[idx_insert - 1].dist_sq;
        if (dist_sq_next <= dist_sq) {
            break;
        }
    }

    // shift up all items by one to make room for new item, drop last item if
    // the queue is full
    int is_full = knque_is_full(queue);
    if (!is_full) {
        // increase size and shift up last element
        items[queue->size] = items[idx_last];
        ++(queue->size);
    }
    for (long idx_dest = idx_last; idx_dest > idx_insert; --idx_dest) {
        items[idx_dest] = items[idx_dest - 1];
    }
    items[idx_insert].value = *value;
    items[idx_insert].dist_sq = dist_sq;

    if (is_full) {  // last element in queue has changed
        queue->dist_sq_max = items[queue->size - 1].dist_sq;
    }

    return 0;
}
