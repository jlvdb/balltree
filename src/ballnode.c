#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "point.h"
#include "ballnode.h"
#include "balltree_macros.h"

typedef struct {
    double min;
    double max;
} Limits;

static Limits limits_new();
static void limits_update(Limits *limits, double value);
static double limits_get_range(Limits *limits);

static inline void point_swap(Point *p1, Point *p2);
static inline double point_get_coord(const Point *point, enum Axis axis);
static double ptslc_sum_weights(const PointSlice *);
static void ball_update_radius(Ball *ball, const PointSlice *slice);
static Ball ball_from_ptslc(const PointSlice *slice);
static enum Axis ptslc_get_maxspread_axis(const PointSlice *);
static Point *ptslc_partition(PointSlice *slice, Point *pivot, enum Axis axis);
static Point *ptslc_quickselect(PointSlice *slice, Point *partition, enum Axis axis);
static Point *ptslc_partition_maxspread_axis(PointSlice *slice);


static Limits limits_new() {
    return (Limits){INFINITY, -INFINITY};
}

static void limits_update(Limits *limits, double value) {
    if (value < limits->min) {
        limits->min = value;
    } else if (value > limits->max) {
        limits->max = value;
    }
}

static double limits_get_range(Limits *limits) {
    return limits->max - limits->min;
}

static inline void point_swap(Point *p1, Point *p2) {
    Point temp = *p1;
    *p1 = *p2;
    *p2 = temp;
}

static inline double point_get_coord(const Point *point, enum Axis axis) {
    return *((double*)point + axis);
}

static double ptslc_sum_weights(const PointSlice *slice) {
    double sumw = 0.0;
    for (const Point *point = slice->start; point < slice->end; ++point) {
        sumw += point->weight;
    }
    return sumw;
}

static void ball_update_radius(Ball *ball, const PointSlice *slice) {
    double dist_squared_max = 0.0;
    for (const Point *point = slice->start; point < slice->end; ++point) {
        double dist_squared = EUCLIDEAN_DIST_SQ(point, ball);
        if (dist_squared > dist_squared_max) {
            dist_squared_max = dist_squared;
        }
    }
    ball->radius = sqrt(dist_squared_max);
}

static Ball ball_from_ptslc(const PointSlice *slice) {
    double center_x = 0.0;
    double center_y = 0.0;
    double center_z = 0.0;
    size_t total = 0;

    for (const Point *point = slice->start; point < slice->end; ++point) {
        ++total;
        double scale = (double)total;
        center_x += (point->x - center_x) / scale;
        center_y += (point->y - center_y) / scale;
        center_z += (point->z - center_z) / scale;
    }

    Ball ball = {
        .x = center_x,
        .y = center_y,
        .z = center_z,
        .radius = -1.0,
    };
    ball_update_radius(&ball, slice);
    return ball;
}

static enum Axis ptslc_get_maxspread_axis(const PointSlice *slice) {
    Limits x_lim = limits_new();
    Limits y_lim = limits_new();
    Limits z_lim = limits_new();

    for (const Point *point = slice->start; point < slice->end; ++point) {
        limits_update(&x_lim, point->x);
        limits_update(&y_lim, point->y);
        limits_update(&z_lim, point->z);
    }

    double x_range = limits_get_range(&x_lim);
    double y_range = limits_get_range(&y_lim);
    double z_range = limits_get_range(&z_lim);
    if (x_range > y_range && x_range > z_range) {
        return (enum Axis)X;
    } else if (y_range > z_range) {
        return (enum Axis)Y;
    } else {
        return (enum Axis)Z;
    }
}

static Point *ptslc_partition(PointSlice *slice, Point *pivot, enum Axis axis) {
    Point *last = slice->end - 1;

    double pivot_value = point_get_coord(pivot, axis);
    point_swap(pivot, last);

    Point *partition = slice->start;
    for (Point *point = partition; point < last; ++point) {
        if (point_get_coord(point, axis) < pivot_value) {
            if (partition != point) {
                point_swap(point, partition);
            }
            ++partition;
        }
    }

    point_swap(last, partition);
    return partition;
}

static Point *ptslc_quickselect(PointSlice *slice, Point *partition, enum Axis axis) {
    if (slice->start < slice->end) {
        int64_t pivot_offset = (slice->end - slice->start) / 2;
        Point *pivot = slice->start + pivot_offset;
        pivot = ptslc_partition(slice, pivot, axis);

        // case: the paritioning element falls into the lower value range
        if (pivot < partition) {
            PointSlice subslice = {
                .start = pivot + 1,
                .end = slice->end,
            };
            pivot = ptslc_quickselect(&subslice, partition, axis);
        }
        
        // case: the paritioning element falls into the higher value range
        else if (pivot > partition) {
            PointSlice subslice = {
                .start = slice->start,
                .end = pivot,
            };
            pivot = ptslc_quickselect(&subslice, partition, axis);
        }

        return pivot;
    }
    return NULL;
}

static Point *ptslc_partition_maxspread_axis(PointSlice *slice) {
    enum Axis split_axis = ptslc_get_maxspread_axis(slice);
    int64_t median_offset = (slice->end - slice->start) / 2;
    Point *median = slice->start + median_offset;
    median = ptslc_quickselect(slice, median, split_axis);
    if (median == NULL) {
        EMIT_ERR_MSG(ValueError, "could not determine median element for partitioning");
    }
    return median;
}

BallNode *bnode_build(PointSlice *slice, int leafsize) {
    int64_t num_points = ptslc_get_size(slice);

    BallNode *node = calloc(1, sizeof(BallNode));
    if (node == NULL) {
        EMIT_ERR_MSG(MemoryError, "BallTree node allocation failed");
        return NULL;
    }
    node->ball = ball_from_ptslc(slice);

    // case: leaf node
    if (num_points <= leafsize) {
        node->data = *slice;
        node->is_leaf = 1;
        node->num_points = num_points;
        node->sum_weight = ptslc_sum_weights(slice);
    }

    // case: regular node with childs    
    else {
        // partition points at median of axis with max. value range (split-axis)
        Point *split = ptslc_partition_maxspread_axis(slice);
        if (split == NULL) {
            goto error;
        }
        PointSlice child_slice;

        // create left child from set points of with lower split-axis values
        child_slice.start = slice->start;
        child_slice.end = split;
        node->childs.left = bnode_build(&child_slice, leafsize);
        if (node->childs.left == NULL) {
            goto error;
        }

        // create right child from set of points with larger split-axis values
        child_slice.start = split;
        child_slice.end = slice->end;
        node->childs.right = bnode_build(&child_slice, leafsize);
        if (node->childs.right == NULL) {
            goto error;
        }

        // use number of points and weights computed further down in the leaf nodes
        node->is_leaf = 0;
        node->num_points = node->childs.left->num_points +
                           node->childs.right->num_points;
        node->sum_weight = node->childs.left->sum_weight +
                           node->childs.right->sum_weight;
    }
    return node;

error:
    bnode_free(node);
    return NULL;
}

void bnode_free(BallNode *node) {
    if (!BALLNODE_IS_LEAF(node)) {
        bnode_free(node->childs.left);
        bnode_free(node->childs.right);
    }
    free(node);
}

int bnode_is_leaf(const BallNode *node) {
    return BALLNODE_IS_LEAF(node);
}

PointSlice bnode_get_ptslc(const BallNode *node) {
    if (BALLNODE_IS_LEAF(node)) {
        return node->data;
    } else {
        PointSlice left = bnode_get_ptslc(node->childs.left);
        PointSlice right = bnode_get_ptslc(node->childs.right);
        return (PointSlice){
            .start = left.start,
            .end = right.end,
        };
    }
}
