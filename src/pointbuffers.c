#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "point.h"
#include "balltree_macros.h"

static inline double rand_uniform(double low, double high);

PointBuffer *ptbuf_new(long size) {
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
        EMIT_ERR_MSG(MemoryError, "PointBuffer buffer allocation failed");
        free(buffer);
        return NULL;
    }

    buffer->size = size;
    buffer->points = points;
    return buffer;
}

PointBuffer *ptbuf_from_buffers(
    long size,
    double *x_vals,
    double *y_vals,
    double *z_vals
) {
    PointBuffer *buffer = ptbuf_new(size);
    if (buffer == NULL) {
        return NULL;
    }
    Point *points = buffer->points;
    for (long i = 0; i < size; ++i) {
        points[i] = point_create(x_vals[i], y_vals[i], z_vals[i]);
    }
    return buffer;
}

PointBuffer *ptbuf_from_buffers_weighted(
    long size,
    double *x_vals,
    double *y_vals,
    double *z_vals,
    double *weights
) {
    PointBuffer *buffer = ptbuf_from_buffers(size, x_vals, y_vals, z_vals);
    if (buffer == NULL) {
        return NULL;
    }
    Point *points = buffer->points;
    for (long i = 0; i < size; ++i) {
        points[i].weight = weights[i];
    }
    return buffer;
}

void ptbuf_free(PointBuffer *buffer) {
    if (buffer->points != NULL) {
        free(buffer->points);
    }
    free(buffer);
}

int ptbuf_resize(PointBuffer *buffer, long size) {
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

PointBuffer *ptbuf_gen_random(double low, double high, long num_points) {
    PointBuffer *buffer = ptbuf_new(num_points);
    if (buffer == NULL) {
        return NULL;
    }

    for (long i = 0; i < num_points; ++i) {
        double x = rand_uniform(low, high);
        double y = rand_uniform(low, high);
        double z = rand_uniform(low, high);
        buffer->points[i] = point_create(x, y, z);
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

long ptslc_get_size(const PointSlice *slice) {
    return slice->end - slice->start;
}
