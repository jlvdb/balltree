#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "point.h"
#include "balltree_macros.h"

static inline double rand_uniform(double low, double high);

PointBuffer *ptbuf_new(int64_t size) {
    if (size < 1) {
        EMIT_ERR_MSG(ValueError, "PointBuffer size must be positive");
        return NULL;
    }

    PointBuffer *buffer = malloc(sizeof(PointBuffer));
    if (buffer == NULL) {
        EMIT_ERR_MSG(MemoryError, "PointBuffer allocation failed");
        return NULL;
    }

    size_t n_bytes = size * sizeof(Point);
    Point *points = malloc(n_bytes);
    if (points == NULL) {
        EMIT_ERR_MSG(MemoryError, "PointBuffer.points allocation failed");
        ptbuf_free(buffer);
        return NULL;
    }
    // initialise range index
    for (int64_t i = 0; i < size; ++i) {
        points[i].index = i;
    }

    buffer->size = size;
    buffer->points = points;
    return buffer;
}

void ptbuf_free(PointBuffer *buffer) {
    if (buffer->points != NULL) {
        free(buffer->points);
    }
    free(buffer);
}

int ptbuf_resize(PointBuffer *buffer, int64_t size) {
    if (size < 1) {
        EMIT_ERR_MSG(ValueError, "PointBuffer size must be positive");
        return BTR_FAILED;
    }

    size_t n_bytes = size * sizeof(Point);
    Point *points = realloc(buffer->points, n_bytes);
    if (points == NULL) {
        EMIT_ERR_MSG(MemoryError, "PointBuffer resizing failed");
        return BTR_FAILED;
    }
    // update range index
    if (size > buffer->size) {
        for (int64_t i = buffer->size; i < size; ++i) {
            points[i].index = i;
        }
    }

    buffer->size = size;
    buffer->points = points;
    return BTR_SUCCESS;
}

PointBuffer *ptbuf_copy(const PointBuffer *buffer) {
    PointBuffer *copy = ptbuf_new(buffer->size);
    if (copy == NULL) {
        return NULL;
    }

    size_t n_bytes = buffer->size * sizeof(Point);
    memcpy(copy->points, buffer->points, n_bytes);
    return copy;
}

static inline double rand_uniform(double low, double high) {
    double rand_uniform_normalised = (double)rand() / RAND_MAX;
    return rand_uniform_normalised * (high - low) + low;
}

PointBuffer *ptbuf_gen_random(double low, double high, int64_t num_points) {
    PointBuffer *buffer = ptbuf_new(num_points);
    if (buffer == NULL) {
        return NULL;
    }

    for (int64_t i = 0; i < num_points; ++i) {
        Point *point = buffer->points + i;
        point->x = rand_uniform(low, high);
        point->y = rand_uniform(low, high);
        point->z = rand_uniform(low, high);
        point->weight = 1.0;
    }
    return buffer;
}

PointSlice *ptslc_from_buffer(const PointBuffer *buffer) {
    PointSlice *slice = malloc(sizeof(PointSlice));
    if (slice == NULL) {
        EMIT_ERR_MSG(MemoryError, "PointSlice allocation failed");
        return NULL;
    }

    slice->start = buffer->points;
    slice->end = buffer->points + buffer->size;
    return slice;
}

int64_t ptslc_get_size(const PointSlice *slice) {
    return slice->end - slice->start;
}

double ptslc_sumw_in_radius_sq(
    const PointSlice *slice,
    const Point *ref_point,
    double rad_sq
) {
    double sumw = 0.0;
    for (const Point *point = slice->start; point < slice->end; ++point) {
        double dist_sq = EUCLIDEAN_DIST_SQ(ref_point, point);
        // add point weight if condition is met otherwise zero
        int dist_mask = dist_sq <= rad_sq;
        sumw += point->weight * (double)dist_mask;
    }
    return sumw;
}

double ptslc_sumw_in_range_sq(
    const PointSlice *slice,
    const Point *ref_point,
    double rmin_sq,
    double rmax_sq
) {
    double sumw = 0.0;
    for (const Point *point = slice->start; point < slice->end; ++point) {
        double dist_sq = EUCLIDEAN_DIST_SQ(ref_point, point);
        // add point weight if condition is met otherwise zero
        int dist_mask = rmin_sq < dist_sq && dist_sq <= rmax_sq;
        sumw += point->weight * (double)dist_mask;
    }
    return sumw;
}

double ptslc_dualsumw_in_radius_sq(
    const PointSlice *slice1,
    const PointSlice *slice2,
    double rad_sq
) {
    double sumw = 0.0;
    for (const Point *point = slice1->start; point < slice1->end; ++point) {
        sumw += point->weight * ptslc_sumw_in_radius_sq(slice2, point, rad_sq);
    }
    return sumw;
}

double ptslc_dualsumw_in_range_sq(
    const PointSlice *slice1,
    const PointSlice *slice2,
    double rmin_sq,
    double rmax_sq
) {
    double sumw = 0.0;
    for (const Point *point = slice1->start; point < slice1->end; ++point) {
        sumw += point->weight * ptslc_sumw_in_range_sq(slice2, point, rmin_sq, rmax_sq);
    }
    return sumw;
}
