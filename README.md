# balltree

Fast balltree implementation for 3-dim data with an Euclidean distance norm.
The base implementation is in C and there is a wrapper for python in development.

## `C` usage

```c
#include <stdio.h>
#include <stdlib.h>

#include "point.h"     // point_create, ptbuf_gen_random
#include "balltree.h"  // balltree_build, balltree_count_radius

int main(int argc, char** argv) {
    // uniform random points in range [-1, 1)
    int n_data = 1000000;
    srand(12345);  // seed random generator
    PointBuffer *points = ptbuf_gen_random(-1.0, 1.0, n_data);
    if (points == NULL) return 1;

    // build tree from points with default leaf size
    BallTree *tree = balltree_build(points);
    if (tree == NULL) return 1;

    // query a single point
    Point query_point = point_create(0.0, 0.0, 0.0);
    double query_radius = 0.2;
    double count = balltree_count_radius(tree, &query_point, query_radius);
    printf("pairs in r <= %.1f: %.0f\n", query_radius, count);
    return 0;
}
```

## `Python` usage

```python
import numpy as np
from balltree import BallTree

# uniform random points in range [-1, 1)
n_data = 1_000_000
rng = np.random.default_rng(12345)
x = rng.uniform(-1.0, 1.0, size=n_data)
y = rng.uniform(-1.0, 1.0, size=n_data)
z = rng.uniform(-1.0, 1.0, size=n_data)

# build tree from points with default leaf size
tree = BallTree(x, y, z)

# query a single point
query_point = [0.0, 0.0, 0.0]
query_radius = 0.2
count = tree.count_radius(query_point, query_radius)
print(f"pairs in r <= {query_radius:.1f}: {count:.0f}")
```
