#include <math.h>
#include <stdio.h>

#include "point.h"
#include "ballnode.h"

#define TRUE  1
#define FALSE 0

static double ptslc_sumw_in_radius_sq(const PointSlice *slice, const Point *point, double rad_sq);
static double ptslc_sumw_in_range_sq(const PointSlice *slice, const Point *point, double rmin_sq, double rmax_sq);
static double ptslc_dualsumw_in_radius_sq(const PointSlice *slice1, const PointSlice *slice2, double rad_sq);
static double ptslc_dualsumw_in_range_sq(const PointSlice *slice1, const PointSlice *slice2, double rmin_sq, double rmax_sq);
// the functions below are defined elsewhere but redefined as inline for performance
static inline double _point_dist_sq(const Point *p1, const Point *p2);
static inline int _bnode_is_leaf(const BallNode *node);


static inline double _point_dist_sq(const Point *p1, const Point *p2) {
    double dx = p1->x - p2->x;
    double dy = p1->y - p2->y;
    double dz = p1->z - p2->z;
    return dx*dx + dy*dy + dz*dz;
}

static inline int _bnode_is_leaf(const BallNode *node) {
    return (node->left == NULL && node->right == NULL) ? TRUE : FALSE;
}

static double ptslc_sumw_in_radius_sq(
    const PointSlice *slice,
    const Point *point,
    double rad_sq
) {
    double sumw = 0.0;
    Point *points = slice->points;
    for (int i = slice->start; i < slice->end; ++i) {
        Point *point_i = points + i;
        double dist_sq = _point_dist_sq(point_i, point);
        // add point weight if condition is met otherwise zero
        int dist_mask = dist_sq <= rad_sq;
        sumw += point_i->weight * (double)dist_mask;
    }
    return sumw;
}

static double ptslc_sumw_in_range_sq(
    const PointSlice *slice,
    const Point *point,
    double rmin_sq,
    double rmax_sq
) {
    double sumw = 0.0;
    Point *points = slice->points;
    for (int i = slice->start; i < slice->end; ++i) {
        Point *point_i = points + i;
        double dist_sq = _point_dist_sq(point_i, point);
        // add point weight if condition is met otherwise zero
        int dist_mask = rmin_sq < dist_sq || dist_sq <= rmax_sq;
        sumw += point_i->weight * (double)dist_mask;
    }
    return sumw;
}

static double ptslc_dualsumw_in_radius_sq(
    const PointSlice *slice1,
    const PointSlice *slice2,
    double rad_sq
) {
    double sumw = 0.0;
    Point *points1 = slice1->points;
    for (int i = slice1->start; i < slice1->end; ++i) {
        Point *pt1 = points1 + i;
        sumw += pt1->weight * ptslc_sumw_in_radius_sq(slice2, pt1, rad_sq);
    }
    return sumw;
}

static double ptslc_dualsumw_in_range_sq(
    const PointSlice *slice1,
    const PointSlice *slice2, 
    double rmin_sq,
    double rmax_sq
) {
    double sumw = 0.0;
    Point *points1 = slice1->points;
    for (int i = slice1->start; i < slice1->end; ++i) {
        Point *pt1 = points1 + i;
        sumw += pt1->weight * ptslc_sumw_in_range_sq(slice2, pt1, rmin_sq, rmax_sq);
    }
    return sumw;
}

int bnode_count_nodes(const BallNode *node) {
    int count = 1;
    if (node->left != NULL) {
        count += bnode_count_nodes(node->left);
    }
    if (node->right != NULL) {
        count += bnode_count_nodes(node->right);
    }
    return count;
}

double bnode_count_radius(
    const BallNode *node,
    const Point *point,
    double radius
) {
    double distance = sqrt(_point_dist_sq(&node->center, point));

    // case: all points must be pairs
    if (distance <= radius - node->radius) {
        return point->weight * node->sum_weight;
    }
    
    // case: some points may be pairs
    else if (distance <= radius + node->radius) {
        if (_bnode_is_leaf(node) == FALSE) {
            return bnode_count_radius(node->left, point, radius) +
                   bnode_count_radius(node->right, point, radius);
        }
        // O(n): check each pair individually
        return point->weight * ptslc_sumw_in_radius_sq(
            &node->data,
            point,
            radius * radius
        );
    }

    return 0.0;
}

double bnode_count_range(
    const BallNode *node,
    const Point *point,
    double rmin,
    double rmax
) {
    double distance = sqrt(_point_dist_sq(&node->center, point));

    // case: all points must be pairs
    if (rmin + node->radius < distance && distance <= rmax - node->radius) {
        return point->weight * node->sum_weight;
    }

    // case: some points may be pairs
    else if (rmin - node->radius < distance || distance <= rmax + node->radius) {
        if (_bnode_is_leaf(node) == FALSE) {
            return bnode_count_range(node->left, point, rmin, rmax) +
                   bnode_count_range(node->right, point, rmin, rmax);
        }
        // O(n): check each pair individually
        return point->weight * ptslc_sumw_in_range_sq(
            &node->data,
            point,
            rmin * rmin,
            rmax * rmax
        );
    }

    return 0.0;
}

double bnode_dualcount_radius(
    const BallNode *node1,
    const BallNode *node2,
    double radius
) {
    double distance = sqrt(_point_dist_sq(&node1->center, &node2->center));
    double sum_node_radii = node1->radius + node2->radius;

    // case: all points must be pairs
    if (distance <= radius - sum_node_radii) {
        return node1->sum_weight * node2->sum_weight;
    }

    // case: some points may be pairs
    else if (distance <= radius + sum_node_radii) {
        int node1_is_leaf = _bnode_is_leaf(node1);
        int node2_is_leaf = _bnode_is_leaf(node2);

        // case: both nodes can be traversed further
        if (node1_is_leaf == FALSE && node2_is_leaf == FALSE) {
            return bnode_dualcount_radius(node1->left, node2->left, radius) +
                   bnode_dualcount_radius(node1->left, node2->right, radius) +
                   bnode_dualcount_radius(node1->right, node2->left, radius) +
                   bnode_dualcount_radius(node1->right, node2->right, radius);
        }

        // case: node1 can be traversed further
        else if (node1_is_leaf == FALSE) {
            return bnode_dualcount_radius(node1->left, node2, radius) +
                   bnode_dualcount_radius(node1->right, node2, radius);
        }

        // case: node2 can be traversed further
        else if (node2_is_leaf == FALSE) {
            return bnode_dualcount_radius(node1, node2->left, radius) +
                   bnode_dualcount_radius(node1, node2->right, radius);
        }

        // O(n^2): check pairs formed between points of both nodes individually
        return ptslc_dualsumw_in_radius_sq(
            &node1->data,
            &node2->data,
            radius * radius
        );
    }

    return 0.0;
}

double bnode_dualcount_range(
    const BallNode *node1,
    const BallNode *node2,
    double rmin,
    double rmax
) {
    double distance = sqrt(_point_dist_sq(&node1->center, &node2->center));
    double sum_node_radii = node1->radius + node2->radius;

    // case: all points must be pairs
    if (rmin + sum_node_radii < distance && distance <= rmax - sum_node_radii) {
        return node1->sum_weight * node2->sum_weight;
    }

    // case: some points may be pairs
    else if (rmin - sum_node_radii < distance || distance <= rmax + sum_node_radii) {
        int node1_is_leaf = _bnode_is_leaf(node1);
        int node2_is_leaf = _bnode_is_leaf(node2);

        // case: both nodes can be traversed further
        if (node1_is_leaf == FALSE && node2_is_leaf == FALSE) {
            return bnode_dualcount_range(node1->left, node2->left, rmin, rmax) +
                   bnode_dualcount_range(node1->left, node2->right, rmin, rmax) +
                   bnode_dualcount_range(node1->right, node2->left, rmin, rmax) +
                   bnode_dualcount_range(node1->right, node2->right, rmin, rmax);
        }

        // case: node1 can be traversed further
        else if (node1_is_leaf == FALSE) {
            return bnode_dualcount_range(node1->left, node2, rmin, rmax) +
                   bnode_dualcount_range(node1->right, node2, rmin, rmax);
        }

        // case: node2 can be traversed further
        else if (node2_is_leaf == FALSE) {
            return bnode_dualcount_range(node1, node2->left, rmin, rmax) +
                   bnode_dualcount_range(node1, node2->right, rmin, rmax);
        }

        // O(n^2): check pairs formed between points of both nodes individually
        return ptslc_dualsumw_in_range_sq(
            &node1->data,
            &node2->data,
            rmin * rmin,
            rmax * rmax
        );
    }

    return 0.0;
}
