#ifndef HISTOGRAM_H
#define HISTOGRAM_H

typedef struct {
    long num_bins;
    double *edges;
    double *sum_weight;
} Histogram;

Histogram *hist_new(long);
void hist_free(Histogram *);
long hist_insert(Histogram *, double, double);

#endif /* HISTOGRAM_H */
