#include <stdio.h>
#include <stdlib.h>

#include "histogram.h"
#include "balltree_macros.h"

Histogram *hist_new(long num_bins) {
    if (num_bins < 1) {
        EMIT_ERR_MSG(ValueError, "Histogram num_bins must be positive");
        return NULL;
    }

    Histogram *hist = malloc(sizeof(Histogram));
    if (hist == NULL) {
        EMIT_ERR_MSG(MemoryError, "Histogram allocation failed");
        return NULL;
    }

    double *edges = calloc(num_bins + 1, sizeof(double));
    if (edges == NULL) {
        EMIT_ERR_MSG(MemoryError, "Histogram edges allocation failed");
        hist_free(hist);
        return NULL;
    }
    double *sum_weight = calloc(num_bins, sizeof(double));
    if (sum_weight == NULL) {
        EMIT_ERR_MSG(MemoryError, "Histogram sum_weight allocation failed");
        hist_free(hist);
        return NULL;
    }

    hist->num_bins = num_bins;
    hist->edges = edges;
    hist->sum_weight = sum_weight;
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
    long num_bins = hist->num_bins;
    if (value <= hist->edges[0]) {
        return -1;
    } else if (value > hist->edges[num_bins + 1]) {
        return num_bins;
    }
    for (long edge_idx = 1; edge_idx <= num_bins; ++edge_idx) {
        if (value <= hist->edges[edge_idx]) {
            long bin_idx = edge_idx - 1;
            hist->sum_weight[bin_idx] += weight;
            return bin_idx;
        }
    }
    return -2;  // this should not be reachable
}
