#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "histogram.h"
#include "balltree_macros.h"

DistHistogram *hist_new(long size, double *dist_edges) {
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
    for (long i = 0; i < size; ++i) {
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

long hist_insert(DistHistogram *hist, double distance, double weight) {
    double dist_sq = distance * distance;
    return HISTOGRAM_INSERT_DIST_SQ(hist, dist_sq, weight);
}