#ifndef HISTOGRAM_H
#define HISTOGRAM_H

typedef struct {
    long size;
    double *edges;
    double *sum_weight;
    double max;
} Histogram;

#define HISTOGRAM_INSERT(hist, value, weight) \
    ({ \
        long __bin_idx = -1; \
        if ((value) <= (hist)->max) { \
            for (__bin_idx = 0; __bin_idx <= hist->size; ++__bin_idx) { \
                if ((value) <= (hist)->edges[__bin_idx]) { \
                    (hist)->sum_weight[__bin_idx] += (weight); \
                } \
            } \
        } \
        __bin_idx; \
    })

Histogram *hist_new(long, double *);
void hist_free(Histogram *);
long hist_insert(Histogram *, double, double);

#endif /* HISTOGRAM_H */
