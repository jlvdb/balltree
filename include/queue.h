#ifndef QUEUE_H
#define QUEUE_H

#include <stdint.h>

#include "point.h"

typedef struct {
    int64_t index;
    double distance;
} QueueItem;

typedef struct {
    int64_t capacity;
    int64_t size;
    QueueItem *items;
    double distance_max;
} KnnQueue;

KnnQueue *knque_new(int64_t);
void knque_free(KnnQueue *);
void knque_clear(KnnQueue *queue);
int knque_insert(KnnQueue *, int64_t, double);

inline int knque_is_full(KnnQueue *queue) {
    return queue->capacity == queue->size;
}

inline double knque_get_max_dist(KnnQueue *queue) {
    double limit = queue->distance_max;
    double dynamic = queue->items[queue->capacity - 1].distance;
    return (limit < dynamic) ? limit : dynamic;
}

#endif /* QUEUE_H */
