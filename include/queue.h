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

#endif /* QUEUE_H */
