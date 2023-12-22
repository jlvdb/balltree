#ifndef BINNING_H
#define BINNING_H

struct Binning {
    int n_bins;
    short closed_left;
    short equal_width;
    const double *edges;
};

struct Binning binning_from_edges(const double*, int, int);
const double* binning_get_edges(const struct Binning*);
double* binning_get_centers(const struct Binning*);
double* binning_get_widths(const struct Binning*);
int* binning_count(const struct Binning*, double *restrict, int);
int* binning_apply(const struct Binning*, double *restrict, int);

#endif /* BINNING_H */
