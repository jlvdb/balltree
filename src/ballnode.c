#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "point.h"
#include "ballnode.h"

#define TRUE  1
#define FALSE 0

static inline void point_swap(Point *p1, Point *p2);
static inline double point_get_coord(const Point *point, enum Axis axis);
static double ptslc_sum_weights(const PointSlice *);
static Point ptslc_get_mean(const PointSlice *);
static double ptslc_get_maxdist(const PointSlice *, Point *);
static enum Axis ptslc_get_maxspread_axis(const PointSlice *);
static int ptslc_partition(PointSlice *slice, int pivot_idx, enum Axis axis);
static int ptslc_quickselect(PointSlice *slice, int partition_idx, enum Axis axis);
static int ptslc_partition_maxspread_axis(PointSlice *slice);
static BallNode *bnode_init(PointBuffer *buffer, int start, int end);


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
    Point *points = slice->points;
    for (size_t i = slice->start; i < slice->end; ++i) {
        Point *point_i = points + i;
        sumw += point_i->weight;
    }
    return sumw;
}

static Point ptslc_get_mean(const PointSlice *slice) {
    double center_x = 0.0;
    double center_y = 0.0;
    double center_z = 0.0;
    int total = 0;

    Point *points = slice->points;
    for (int i = slice->start; i < slice->end; ++i) {
        ++total;
        double scale = (double)total;
        Point point = points[i];
        center_x += (point.x - center_x) / scale;
        center_y += (point.y - center_y) / scale;
        center_z += (point.z - center_z) / scale;
    }
    return point_create(center_x, center_y, center_z);
}

static double ptslc_get_maxdist(const PointSlice *slice, Point *center) {
    Point *points = slice->points;
    double dist_squared_max = 0.0;
    for (int i = slice->start; i < slice->end; ++i) {
        double dist_squared = point_dist_sq(points + i, center);
        if (dist_squared > dist_squared_max) {
            dist_squared_max = dist_squared;
        }
    }
    return sqrt(dist_squared_max);
}

static enum Axis ptslc_get_maxspread_axis(const PointSlice *slice) {
    double x_min = INFINITY;
    double y_min = x_min;
    double z_min = x_min;

    double x_max = -INFINITY;
    double y_max = x_max;
    double z_max = x_max;

    Point *points = slice->points;
    double xi, yi, zi;
    for (int i = slice->start; i < slice->end; ++i) {
        Point point = points[i];

        xi = point.x;
        if (xi < x_min) {
            x_min = xi;
        } else if (xi > x_max) {
            x_max = xi;
        }
        yi = point.y;
        if (yi < y_min) {
            y_min = yi;
        } else if (yi > y_max) {
            y_max = yi;
        }
        zi = point.z;
        if (zi < z_min) {
            z_min = zi;
        } else if (zi > z_max) {
            z_max = zi;
        }
    }

    double x_spread = x_max - x_min;
    double y_spread = y_max - y_min;
    double z_spread = z_max - z_min;
    if (x_spread > y_spread && x_spread > z_spread) {
        return (enum Axis)X;
    } else if (y_spread > z_spread) {
        return (enum Axis)Y;
    } else {
        return (enum Axis)Z;
    }
}

static int ptslc_partition(PointSlice *slice, int pivot_idx, enum Axis axis) {
    Point *points = slice->points;
    int last_idx = slice->end - 1;

    double pivot = point_get_coord(points + pivot_idx, axis);
    point_swap(points + pivot_idx, points + last_idx);

    int partition_idx = slice->start;
    for (int i = slice->start; i < last_idx; ++i) {
        if (point_get_coord(points + i, axis) < pivot) {
            if (partition_idx != i) {
                point_swap(points + i, points + partition_idx);
            }
            ++partition_idx;
        }
    }

    point_swap(points + last_idx, points + partition_idx);
    return partition_idx;
}

static int ptslc_quickselect(PointSlice *slice, int partition_idx, enum Axis axis) {
    if (slice->start < slice->end) {
        int pivot_idx = (slice->start + slice->end) / 2;
        pivot_idx = ptslc_partition(slice, pivot_idx, axis);

        if (pivot_idx < partition_idx) {
            PointSlice subslice = {
                .start = pivot_idx + 1,
                .end = slice->end,
                .points = slice->points,
            };
            pivot_idx = ptslc_quickselect(&subslice, partition_idx, axis);
        }
        
        else if (pivot_idx > partition_idx) {
            PointSlice subslice = {
                .start = slice->start,
                .end = pivot_idx,
                .points = slice->points,
            };
            pivot_idx = ptslc_quickselect(&subslice, partition_idx, axis);
        }

        return pivot_idx;
    }
    return -1;
}

static int ptslc_partition_maxspread_axis(PointSlice *slice) {
    enum Axis split_axis = ptslc_get_maxspread_axis(slice);
    int median_idx = (slice->end + slice->start) / 2;
    return ptslc_quickselect(slice, median_idx, split_axis);
}

static BallNode *bnode_init(PointBuffer *buffer, int start, int end) {
    BallNode *node = (BallNode *)calloc(1, sizeof(BallNode));
    if (node == NULL) {
        fprintf(stderr, "ERROR: BallTree node allocation failed\n");
        return NULL;
    }

    node->data = (PointSlice){
        .start = start,
        .end = end,
        .points = buffer->points,
    };
    node->center = ptslc_get_mean(&node->data);
    node->radius = ptslc_get_maxdist(&node->data, &node->center);
    return node;
}

BallNode *bnode_build(PointBuffer *buffer, int start, int end, int leafsize) {
    BallNode *node = bnode_init(buffer, start, end);
    if (node == NULL) {
        return NULL;
    }
    PointSlice *slice = &node->data;

    if (end - start <= leafsize) {
        // finalize as leaf node
        node->sum_weight = ptslc_sum_weights(slice);
    }
    
    else {
        // partition points and create children
        int split_idx = ptslc_partition_maxspread_axis(slice);
        if (split_idx == -1) {
            fprintf(stderr, "ERROR: could not determine median element for partitioning\n");
            bnode_free(node);
            return NULL;
        }

        // create left child from the set points of with lower split-axis values
        node->left = bnode_build(buffer, start, split_idx, leafsize);
        if (node->left == NULL) {
            bnode_free(node);
            return NULL;
        }

        // create right child from the set of points with larger split-axis values
        node->right = bnode_build(buffer, split_idx, end, leafsize);
        if (node->right == NULL) {
            bnode_free(node);
            return NULL;
        }

        // use weights computed further down in the leaf nodes
        node->sum_weight = node->left->sum_weight + node->right->sum_weight;
    }
    return node;
}

void bnode_free(BallNode *node) {
    if (node->left != NULL) {
        bnode_free(node->left);
    }
    if (node->right != NULL) {
        bnode_free(node->right);
    }
    free(node);
}

int bnode_is_leaf(const BallNode *node) {
    return (node->left == NULL && node->right == NULL) ? TRUE : FALSE;
}
