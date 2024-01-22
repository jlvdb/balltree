#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "binning.h"

#define EQUAL_WIDTH_REL_TOL 1e-12

struct Binning binning_from_edges(const double *edges, int n_bins, int closed_left)
{
    int equal_width = true;
    double last_width = edges[1] - edges[0];
    for (int i = 1; i < n_bins; ++i) {
        double left = edges[i];
        double right = edges[i+1];
        double width = right - left;
        if (fabs((width - last_width) / last_width) > EQUAL_WIDTH_REL_TOL) {
            equal_width = false;
            break;
        }
    }
    return (struct Binning){
        .n_bins = n_bins,
        .closed_left = closed_left,
        .equal_width = equal_width,
        .edges = edges,
    };
}

const double* binning_get_edges(const struct Binning *bins)
{
    return bins->edges;
}

double* binning_get_centers(const struct Binning *bins)
{
    double *centers = malloc(bins->n_bins * sizeof(double));
    if (!centers) {
        return NULL;
    }
    for (int i = 0; i < bins->n_bins; ++i) {
        centers[i] = (bins->edges[i] + bins->edges[i+1]) * 0.5;
    }
    return centers;
}

double* binning_get_widths(const struct Binning *bins)
{
    double *widths = malloc(bins->n_bins * sizeof(double));
    if (!widths) {
        return NULL;
    }
    for (int i = 0; i < bins->n_bins; ++i) {
        widths[i] = bins->edges[i+1] - bins->edges[i];
    }
    return widths;
}

int* binning_count(const struct Binning *bins, double *restrict data, int size)
{
    int *count = calloc(bins->n_bins, sizeof(int));
    if (!count) {
        return NULL;
    }
    const double *restrict edges = bins->edges;
    if (bins->equal_width) {
        double min_value = bins->edges[0];
        double bin_width = bins->edges[1] - min_value;
        for (int i = 0; i < size; ++i) {
            int bin_index = (data[i] - min_value) / bin_width;
            ++count[bin_index];
        }
    } else if (bins->closed_left) {
        for (int i = 0; i < size; ++i) {
            double value = data[i];
            for (int bin_index = 0; bin_index < bins->n_bins; ++bin_index) {
                if ((edges[bin_index] <= value) && (value < edges[bin_index+1])) {
                    ++count[bin_index];
                    break;
                }
            }
        }
    } else {
        for (int i = 0; i < size; ++i) {
            double value = data[i];
            for (int bin_index = 0; bin_index < bins->n_bins; ++bin_index) {
                if ((edges[bin_index] < value) && (value <= edges[bin_index+1])) {
                    ++count[bin_index];
                    break;
                }
            }
        }
    }
    return count;
}

int* binning_apply(const struct Binning *bins, double *restrict data, int size)
{
    int *count = malloc(sizeof(int) * size);
    if (!count) {
        return NULL;
    }
    const double *restrict edges = bins->edges;
    int n_bins = bins->n_bins;
    if (bins->equal_width) {
        double min_value = bins->edges[0];
        double bin_width = bins->edges[1] - min_value;
        for (int i = 0; i < size; ++i) {
            int bin_index = (data[i] - min_value) / bin_width;
            if (bin_index < 0 || bin_index >= n_bins) {
                count[i] = -1;
            } else {
                count[i] = bin_index;
            }
        }
    } else if (bins->closed_left) {
        for (int i = 0; i < size; ++i) {
            double value = data[i];
            count[i] = -1;
            for (int bin_index = 0; bin_index < n_bins; ++bin_index) {
                if ((edges[bin_index] <= value) && (value < edges[bin_index+1])) {
                    count[i] = bin_index;
                    break;
                }
            }
        }
    } else {
        for (int i = 0; i < size; ++i) {
            double value = data[i];
            count[i] = -1;
            for (int bin_index = 0; bin_index < n_bins; ++bin_index) {
                if ((edges[bin_index] < value) && (value <= edges[bin_index+1])) {
                    count[i] = bin_index;
                    break;
                }
            }
        }
    }
    return count;
}

/*
int main()
{
    const int nbins = 8;
    double edges[nbins+1] = {0.000, 1.000, 2.000, 3.000, 4.000, 5.000, 6.000, 7.000, 8.0001};
    struct Binning bins = binning_from_edges(edges, nbins, 1);

    double data[10] = {1,1,2,2,3,3,3,3,3,4};
    int *counts = binning_count(&bins, data, 10);
    for (int i = 0; i < bins.n_bins; ++i) {
        printf("%02d -> %2d\n", i, counts[i]);
    }
}
*/
