#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "histogram.h"
#include "balltree_macros.h"

Histogram *hist_new(long num_edges, double *bin_edges) {
    if (num_edges < 2) {
        EMIT_ERR_MSG(ValueError, "Histogram requires at least 2 edges");
        return NULL;
    }

    Histogram *hist = malloc(sizeof(Histogram));
    if (hist == NULL) {
        EMIT_ERR_MSG(MemoryError, "Histogram allocation failed");
        return NULL;
    }

    size_t n_bytes = num_edges * sizeof(double);
    double *edges = malloc(n_bytes);
    if (edges == NULL) {
        EMIT_ERR_MSG(MemoryError, "Histogram edges allocation failed");
        hist_free(hist);
        return NULL;
    }
    for (long i = 0; i < num_edges; ++i) {
        double edge = bin_edges[i];
        edges[i] = edge * edge;
    }

    double *sum_weight = calloc(num_edges - 1, sizeof(double));
    if (sum_weight == NULL) {
        EMIT_ERR_MSG(MemoryError, "Histogram sum_weight allocation failed");
        hist_free(hist);
        return NULL;
    }

    hist->num_bins = num_edges - 1;
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
