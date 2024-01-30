#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "histogram.h"
#include "balltree_macros.h"

Histogram *hist_new(long size, double *bin_edges) {
    if (size < 1) {
        EMIT_ERR_MSG(ValueError, "Histogram requires at least 1 edges");
        return NULL;
    }

    Histogram *hist = malloc(sizeof(Histogram));
    if (hist == NULL) {
        EMIT_ERR_MSG(MemoryError, "Histogram allocation failed");
        return NULL;
    }

    size_t n_bytes = size * sizeof(double);
    double *edges = malloc(n_bytes);
    if (edges == NULL) {
        EMIT_ERR_MSG(MemoryError, "Histogram edges allocation failed");
        hist_free(hist);
        return NULL;
    }
    for (long i = 0; i < size; ++i) {
        double edge = bin_edges[i];
        edges[i] = edge * edge;  // TODO: danger!
    }

    double *sum_weight = calloc(size, sizeof(double));
    if (sum_weight == NULL) {
        EMIT_ERR_MSG(MemoryError, "Histogram sum_weight allocation failed");
        hist_free(hist);
        return NULL;
    }

    hist->size = size;
    hist->edges = edges;
    hist->sum_weight = sum_weight;
    hist->max = edges[size - 1];
    return hist;
}

void hist_free(Histogram *hist) {
    if (hist->edges != NULL) {
        free(hist->edges);
    }
    if (hist->sum_weight != NULL) {
        free(hist->sum_weight);
    }
    free(hist);
}

long hist_insert(Histogram *hist, double value, double weight) {
    return HISTOGRAM_INSERT(hist, value, weight);
}
