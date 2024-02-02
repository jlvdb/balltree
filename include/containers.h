#ifndef CONTAINERS_H
#define CONTAINERS_H

#include <stdint.h>

#include "point.h"

#define MIN(a, b) ((a) < (b)) ? (a) : (b)

typedef struct {
    int64_t size;
    double *sum_weight;
    double *dist;
    double dist_max;
    double *dist_sq;
    double dist_sq_max;
} DistHistogram;

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

DistHistogram *hist_new(int64_t, double *);
void hist_free(DistHistogram *);
int64_t hist_insert_dist_sq(DistHistogram *, double, double);

KnnQueue *knque_new(int64_t);
void knque_free(KnnQueue *);
void knque_clear(KnnQueue *queue);
int knque_insert(KnnQueue *, int64_t, double);

inline double knque_get_max_dist(KnnQueue *queue) {
    return MIN(queue->distance_max, queue->items[queue->capacity - 1].distance);
}

inline int knque_is_full(KnnQueue *queue) {
    return queue->capacity == queue->size;
}

#endif /* CONTAINERS_H */
