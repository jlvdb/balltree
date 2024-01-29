#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#include "point.h"
#include "ballnode.h"
#include "balltree_macros.h"

double bnode_count_radius(
    const BallNode *node,
    const Point *point,
    double radius
) {
    double distance = sqrt(EUCLIDEAN_DIST_SQ(&node->ball, point));
    double node_radius = node->ball.radius;

    // case: all points must be pairs
    if (distance <= radius - node_radius) {
        return point->weight * node->sum_weight;
    }

    // case: some points may be pairs
    else if (distance <= radius + node_radius) {
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

    return 0.0;
}

double bnode_count_range(
    const BallNode *node,
    const Point *point,
    double rmin,
    double rmax
) {
    double distance = sqrt(EUCLIDEAN_DIST_SQ(&node->ball, point));
    double node_radius = node->ball.radius;

    // case: all points must be pairs
    if (rmin + node_radius < distance && distance <= rmax - node_radius) {
        return point->weight * node->sum_weight;
    }

    // case: some points may be pairs
    else if (rmin - node_radius < distance || distance <= rmax + node_radius) {
        if (BALLNODE_IS_LEAF(node) == false) {
            return bnode_count_range(node->childs.left, point, rmin, rmax) +
                   bnode_count_range(node->childs.right, point, rmin, rmax);
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
    double distance = sqrt(EUCLIDEAN_DIST_SQ(&node1->ball, &node2->ball));
    double sum_node_radii = node1->ball.radius + node2->ball.radius;

    // case: all points must be pairs
    if (distance <= radius - sum_node_radii) {
        return node1->sum_weight * node2->sum_weight;
    }

    // case: some points may be pairs
    else if (distance <= radius + sum_node_radii) {
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

    return 0.0;
}

double bnode_dualcount_range(
    const BallNode *node1,
    const BallNode *node2,
    double rmin,
    double rmax
) {
    double distance = sqrt(EUCLIDEAN_DIST_SQ(&node1->ball, &node2->ball));
    double sum_node_radii = node1->ball.radius + node2->ball.radius;

    // case: all points must be pairs
    if (rmin + sum_node_radii < distance && distance <= rmax - sum_node_radii) {
        return node1->sum_weight * node2->sum_weight;
    }

    // case: some points may be pairs
    else if (rmin - sum_node_radii < distance || distance <= rmax + sum_node_radii) {
        int node1_is_leaf = BALLNODE_IS_LEAF(node1);
        int node2_is_leaf = BALLNODE_IS_LEAF(node2);

        // case: both nodes can be traversed further
        if (node1_is_leaf == false && node2_is_leaf == false) {
            return bnode_dualcount_range(node1->childs.left, node2->childs.left, rmin, rmax) +
                   bnode_dualcount_range(node1->childs.left, node2->childs.right, rmin, rmax) +
                   bnode_dualcount_range(node1->childs.right, node2->childs.left, rmin, rmax) +
                   bnode_dualcount_range(node1->childs.right, node2->childs.right, rmin, rmax);
        }

        // case: node1 can be traversed further
        else if (node1_is_leaf == false) {
            return bnode_dualcount_range(node1->childs.left, node2, rmin, rmax) +
                   bnode_dualcount_range(node1->childs.right, node2, rmin, rmax);
        }

        // case: node2 can be traversed further
        else if (node2_is_leaf == false) {
            return bnode_dualcount_range(node1, node2->childs.left, rmin, rmax) +
                   bnode_dualcount_range(node1, node2->childs.right, rmin, rmax);
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
