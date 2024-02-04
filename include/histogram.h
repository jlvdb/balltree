#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <stdint.h>

typedef struct {
    int64_t size;
    double *sum_weight;
    double *dist;
    double dist_max;
    double *dist_sq;
    double dist_sq_max;
} DistHistogram;

DistHistogram *hist_new(int64_t, double *);
void hist_free(DistHistogram *);
int64_t hist_insert_dist_sq(DistHistogram *, double, double);

#endif /* HISTOGRAM_H */
