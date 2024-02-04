Usage Examples
==============


Below are two examples that illustrate how to use the ball tree from `C` and
`Python` to count all pairs with a maximum separation of 0.2 that can be found
in a box of side length 2 with one million randomly generated data points.


Using the C library
-------------------

.. code-block:: c

    #include <stdio.h>
    #include <stdlib.h>

    #include "point.h"     // point_create, ptbuf_gen_random
    #include "balltree.h"  // balltree_build, balltree_count_radius

    int main(int argc, char** argv) {
        // uniform random points with (x, y, z) coordinates in range [-1, 1)
        int n_data = 1000000;
        srand(12345);  // seed random generator
        PointBuffer *buffer = ptbuf_gen_random(-1.0, 1.0, n_data);
        if (buffer == NULL) return 1;

        // build tree from points with default leaf size
        BallTree *tree = balltree_build(buffer);
        if (tree == NULL) return 1;

        // count neighbours of all points
        double query_radius = 0.2;
        double count;
        for (long i = 0; i < buffer->size; ++i) {
            Point *query_point = buffer->points + i;
            count += balltree_count_radius(tree, query_point, query_radius);
        }
        printf("pairs in r <= %.1f: %.0f\n", query_radius, count);

        // count neighbours of all points using the efficient dual-tree algorithm
        count = balltree_dualcount_radius(tree, tree, query_radius);
        printf("pairs in r <= %.1f: %.0f\n", query_radius, count);

        return 0;
    }


Using the Python wrapper
------------------------

.. code-block:: python

    import numpy as np
    from balltree import BallTree


    if __name__ == "__main__":
        # uniform random points with (x, y, z) coordinates in range [-1, 1)
        n_data = 1_000_000
        rng = np.random.default_rng(12345)
        points = rng.uniform(-1.0, 1.0, size=(n_data, 3))

        # build tree from points with default leaf size
        tree = BallTree(points)

        # count neighbours of all points
        query_radius = 0.2
        count = tree.count_radius(points, query_radius)
        print(f"pairs in r <= {query_radius:.1f}: {count:.0f}")

        # count neighbours of all points using the efficient dual-tree algorithm
        count = tree.dualcount_radius(tree, query_radius)
        print(f"pairs in r <= {query_radius:.1f}: {count:.0f}")
