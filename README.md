# balltree

Fast balltree implementation for 3-dim data with an Euclidean distance norm.
The base implementation is in C and there is a wrapper for python in development.

## `C` usage

```c
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "balltree.h"

PointBuffer *create_random_points(int num_points) {
    srand(time(NULL));
    PointBuffer *points = ptbuf_new(num_points);
    if (points == NULL) {
        return NULL;
    }
    for (int i = 0; i < num_points; ++i) {
        double x = (double)rand() / RAND_MAX;
        double y = (double)rand() / RAND_MAX;
        double z = (double)rand() / RAND_MAX;
        points[i] = point_create(x, y, z);
    }
    return points;
}

int main(int argc, char** argv) {
    // build tree from random points in range [1, 0)
    PointBuffer *buffer = create_random_points(1'000'000);
    if (buffer == NULL) {
        return 1;
    }
    BallTree *tree = balltree_build(buffer);
    if (tree == NULL) {
        return 1;
    }

    // query a single point
    Point query_point = point_create(0.5, 0.5, 0.5);
    double query_radius = 0.1;
    double count = balltree_count_radius(tree, &query_point, query_radius);
    printf("pairs in r <= %.1f: %.0f\n", query_radius, count);
    return 0;
}
```

## `Python` usage

```python
import numpy as np
from balltree import BallTree

# build tree from random points in range [1, 0)
n_data = 1_000_000
x = np.random.uniform(size=n_data)
y = np.random.uniform(size=n_data)
z = np.random.uniform(size=n_data)
tree = BallTree.from_data(x, y, z)

# query a single point
query_point = [0.5, 0.5, 0.5]
query_radius = 0.1
count = tree.count_radius(query_point, query_radius)
print(f"pairs in r <= {query_radius:.1f}: {count:.0f}")
```
