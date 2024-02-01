#ifndef QUEUE_H
#define QUEUE_H

#include <stdint.h>

#include "point.h"

# define QUEUEITEM_T Point
# define QUEUEITEM_DEFAULT (Point){0.0, 0.0, 0.0, 0.0}

typedef struct {
    QUEUEITEM_T value;
    double dist_sq;
} QueueItem;

typedef struct {
    int64_t capacity;
    int64_t size;
    QueueItem *items;
    double dist_sq_max;
} KnnQueue;

KnnQueue *knque_new(int64_t);
void knque_free(KnnQueue *);
void knque_clear(KnnQueue *queue);
int knque_insert(KnnQueue *, const QUEUEITEM_T *, double);

#endif /* QUEUE_H */
