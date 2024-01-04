#include <stdio.h>

#define SUCCESS 0
#define FAILED  1

#include "point.h"

PointBuffer *pointbuffer_new(int size) {
    if (size < 1) {
        fprintf(stderr, "ERROR: PointBuffer size must be positive\n");
        return NULL;
    }

    PointBuffer *buffer = (PointBuffer *)malloc(sizeof(PointBuffer));
    if (buffer == NULL) {
        fprintf(stderr, "ERROR: PointBuffer allocation failed\n");
        return NULL;
    }

    size_t n_bytes = size * sizeof(Point);
    Point *points = (Point *)malloc(n_bytes);
    if (points == NULL) {
        fprintf(stderr, "ERROR: PointBuffer buffer allocation failed\n");
        free(buffer);
        return NULL;
    }

    buffer->size = size;
    buffer->points = points;
    return buffer;
}

void pointbuffer_free(PointBuffer *buffer) {
    if (buffer->points != NULL) {
        free(buffer->points);
    }
    free(buffer);
}

int pointbuffer_resize(PointBuffer *buffer, int size) {
    if (size < 1) {
        fprintf(stderr, "ERROR: PointBuffer size must be positive\n");
        return FAILED;
    }

    size_t n_bytes = size * sizeof(Point);
    Point *points = (Point*)realloc(buffer->points, n_bytes);
    if (points == NULL) {
        fprintf(stderr, "ERROR: PointBuffer resizing failed\n");
        return FAILED;
    }

    buffer->size = size;
    buffer->points = points;
    return SUCCESS;
}

PointSlice *pointslice_from_buffer(const PointBuffer *buffer) {
    PointSlice *slice = (PointSlice*)malloc(sizeof(PointSlice));
    if (slice == NULL) {
        fprintf(stderr, "ERROR: PointSlice allocation failed\n");
        return NULL;
    }

    slice->start = 0;
    slice->end = buffer->size;
    slice->points = buffer->points;
    return slice;
}

int pointslice_get_size(const PointSlice *slice) {
    return slice->end - slice->start;
}
