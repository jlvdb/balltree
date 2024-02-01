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
    queue->distance_max = INFINITY;
    for (int64_t i = 0; i < queue->capacity; ++i) {
        queue->items[i] = (QueueItem){
            .index = -1,
            .distance = INFINITY,
        };
    }
}

int knque_insert(KnnQueue *queue, int64_t item_index, double distance) {
    if (distance >= queue->distance_max) {
        return 1;
    }
    // insertion index must now be at least at end of queue

    // the actual location is reached when the element at the preceeding index
    // is smaller than the item to be inserted
    QueueItem *items = queue->items;
    int64_t idx_last = (queue->size > 0) ? (queue->size - 1) : 0;
    int64_t idx_insert = idx_last;
    for (; idx_insert >= 0; --idx_insert) {
        if (items[idx_insert - 1].distance <= distance) {
            break;
        }
    }

    // shift up all items by one to make room for new item, drop last item if
    // the queue is full
    int is_full = queue->capacity == queue->size;
    if (!is_full) {
        // increase size and shift up last element
        items[queue->size] = items[idx_last];
        ++(queue->size);
    }
    for (int64_t idx_dest = idx_last; idx_dest > idx_insert; --idx_dest) {
        items[idx_dest] = items[idx_dest - 1];
    }
    items[idx_insert].index = item_index;
    items[idx_insert].distance = distance;

    if (is_full) {  // last element in queue has changed
        queue->distance_max = items[queue->size - 1].distance;
    }

    return 0;
}
