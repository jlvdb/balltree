#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

#include "point.h"
#include "histogram.h"
#include "ballnode.h"
#include "balltree_macros.h"

static inline void ptslc_sumw_in_hist_sq(const PointSlice *slice, const Point *ref_point, DistHistogram *hist);
static inline void ptslc_dualsumw_in_hist_sq(const PointSlice *slice1, const PointSlice *slice2, DistHistogram *hist);


static inline void ptslc_sumw_in_hist_sq(
    const PointSlice *slice,
    const Point *ref_point,
    DistHistogram *hist
) {
    for (const Point *point = slice->start; point < slice->end; ++point) {
        double dist_sq = EUCLIDEAN_DIST_SQ(ref_point, point);
        // increment corresponding bin by weight
        // TODO: improvements in ptslc_dualsumw_in_hist_sq allow to drop the
        // multiplication with the ref_point weight
        hist_insert_dist_sq(hist, dist_sq, point->weight * ref_point->weight);
    }
}

static inline void ptslc_dualsumw_in_hist_sq(
    const PointSlice *slice1,
    const PointSlice *slice2,
    DistHistogram *hist
) {
    // TODO: maybe a temporary histogram is better, scale with point weight and
    // clear sum_weight after each iteration and add counts to the total
    for (const Point *point = slice1->start; point < slice1->end; ++point) {
        ptslc_sumw_in_hist_sq(slice2, point, hist);
    }
}

double bnode_count_radius(
    const BallNode *node,
    const Point *point,
    double radius
) {
    double distance = sqrt(EUCLIDEAN_DIST_SQ(&node->ball, point));
    double node_radius = node->ball.radius;

    // case: node does not overlap with any bin
    if (distance > radius + node_radius) {
        return 0.0;
    }

    // case: node entirely overlaps with radius
    if (distance <= radius - node_radius) {
        return point->weight * node->sum_weight;
    }

    // case: node partialy overlaps with radius
    if (BALLNODE_IS_LEAF(node) == false) {
        return bnode_count_radius(node->childs.left, point, radius) +
                bnode_count_radius(node->childs.right, point, radius);
    }
    // O(n): check each pair individually
    return point->weight * ptslc_sumw_in_radius_sq(
        &node->data,
        point,
        radius * radius
    );
}

void bnode_count_range(
    const BallNode *node,
    const Point *point,
    DistHistogram *hist
) {
    double distance = sqrt(EUCLIDEAN_DIST_SQ(&node->ball, point));
    double node_radius = node->ball.radius;

    // case: node does not overlap with any bin
    if (distance > hist->dist_max + node_radius) {
        return;
    }

    // case: node may entirely fall into one bin
    double rmin = -INFINITY;  // ensure 0.0 is included at first iteration
    for (int64_t i = 0; i < hist-> size; ++i) {
        double rmax = hist->dist[i];
        if (rmin + node_radius < distance && distance <= rmax - node_radius) {
            hist->sum_weight[i] += point->weight * node->sum_weight;
            return;
        }
        rmin = rmax;
    }

    // case: node overlaps with multiple bins
    if (BALLNODE_IS_LEAF(node) == false) {
        bnode_count_range(node->childs.left, point, hist);
        bnode_count_range(node->childs.right, point, hist);
        return;
    }

    // O(n): check each pair individually
    ptslc_sumw_in_hist_sq(&node->data, point, hist);
}

double bnode_dualcount_radius(
    const BallNode *node1,
    const BallNode *node2,
    double radius
) {
    double distance = sqrt(EUCLIDEAN_DIST_SQ(&node1->ball, &node2->ball));
    double sum_node_radii = node1->ball.radius + node2->ball.radius;

    // case: nodes do not overlap within radius
    if (distance > radius + sum_node_radii) {
        return 0.0;
    }

    // case: nodes entirely overlap within radius
    if (distance <= radius - sum_node_radii) {
        return node1->sum_weight * node2->sum_weight;
    }

    // case: nodes partialy overlap within radius
    int node1_is_leaf = BALLNODE_IS_LEAF(node1);
    int node2_is_leaf = BALLNODE_IS_LEAF(node2);

    // case: both nodes can be traversed further
    if (node1_is_leaf == false && node2_is_leaf == false) {
        return bnode_dualcount_radius(node1->childs.left, node2->childs.left, radius) +
                bnode_dualcount_radius(node1->childs.left, node2->childs.right, radius) +
                bnode_dualcount_radius(node1->childs.right, node2->childs.left, radius) +
                bnode_dualcount_radius(node1->childs.right, node2->childs.right, radius);
    }

    // case: node1 can be traversed further
    else if (node1_is_leaf == false) {
        return bnode_dualcount_radius(node1->childs.left, node2, radius) +
                bnode_dualcount_radius(node1->childs.right, node2, radius);
    }

    // case: node2 can be traversed further
    else if (node2_is_leaf == false) {
        return bnode_dualcount_radius(node1, node2->childs.left, radius) +
                bnode_dualcount_radius(node1, node2->childs.right, radius);
    }

    // O(n^2): check pairs formed between points of both nodes individually
    return ptslc_dualsumw_in_radius_sq(
        &node1->data,
        &node2->data,
        radius * radius
    );
}

void bnode_dualcount_range(
    const BallNode *node1,
    const BallNode *node2,
    DistHistogram *hist
) {
    double distance = sqrt(EUCLIDEAN_DIST_SQ(&node1->ball, &node2->ball));
    double sum_node_radii = node1->ball.radius + node2->ball.radius;

    // case: nodes do not overlap within range of bins
    if (distance > hist->dist_max + sum_node_radii) {
        return;
    }

    // case: nodes may entirely fall into one bin
    double rmin = -INFINITY;  // ensure 0.0 is included at first iteration
    for (int64_t i = 0; i < hist-> size; ++i) {
        double rmax = hist->dist[i];
        if (rmin + sum_node_radii < distance && distance <= rmax - sum_node_radii) {
            hist->sum_weight[i] += node1->sum_weight * node2->sum_weight;
            return;
        }
        rmin = rmax;
    }

    // case: nodes overlaps with multiple bins
    int node1_is_leaf = BALLNODE_IS_LEAF(node1);
    int node2_is_leaf = BALLNODE_IS_LEAF(node2);

    // case: both nodes can be traversed further
    if (node1_is_leaf == false && node2_is_leaf == false) {
        bnode_dualcount_range(node1->childs.left, node2->childs.left, hist);
        bnode_dualcount_range(node1->childs.left, node2->childs.right, hist);
        bnode_dualcount_range(node1->childs.right, node2->childs.left, hist);
        bnode_dualcount_range(node1->childs.right, node2->childs.right, hist);
        return;
    }

    // case: node1 can be traversed further
    else if (node1_is_leaf == false) {
        bnode_dualcount_range(node1->childs.left, node2, hist);
        bnode_dualcount_range(node1->childs.right, node2, hist);
        return;
    }

    // case: node2 can be traversed further
    else if (node2_is_leaf == false) {
        bnode_dualcount_range(node1, node2->childs.left, hist);
        bnode_dualcount_range(node1, node2->childs.right, hist);
        return;
    }

    // O(n^2): check pairs formed between points of both nodes individually
    return ptslc_dualsumw_in_hist_sq(&node1->data, &node2->data, hist);
}
