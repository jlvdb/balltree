#ifndef HISTOGRAM_H
#define HISTOGRAM_H

typedef struct {
    long size;
    double *sum_weight;
    double *dist;
    double dist_max;
    double *dist_sq;
    double dist_sq_max;
} DistHistogram;

#define HISTOGRAM_INSERT_DIST_SQ(hist, dist_sq, weight) \
    ({ \
        long __bin_idx = -1; \
        if ((dist_sq) <= (hist)->dist_sq_max) { \
            for (__bin_idx = 0; __bin_idx <= hist->size; ++__bin_idx) { \
                if ((dist_sq) <= (hist)->dist_sq[__bin_idx]) { \
                    (hist)->sum_weight[__bin_idx] += (weight); \
                    break; \
                } \
            } \
        } \
        __bin_idx; \
    })

DistHistogram *hist_new(long, double *);
void hist_free(DistHistogram *);
long hist_insert(DistHistogram *, double, double);

#endif /* HISTOGRAM_H */
