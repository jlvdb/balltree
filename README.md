# balltree

Fast balltree implementation for 3-dim (weighted) data with an Euclidean
distance norm. The base implementation is in `C` and there is a wrapper for
`Python`.

The tree is optimised towards spatial correlation function calculations since
the query routines are geared towards range queries, i.e. counting pairs with a
given (range of) separations. Fixed number nearest neighbour search is currently
not implemented.

Range queries are typically 25-30x faster than the corresponding implementation
in `scipy.spatial.KDTree` (see below).

#### Contents
[Installation](#installation)  
[Usage](#usage)  
[Comparison to scipy.spatial.KDTree](#comparison-to-scipyspatialkdtree)  
[Maintainers](#aintainers)  


## Installation

A `C` library can be built with the provided make file, the python wrapper is
automatically compiled and installed with `pip install .`.

The installation does not require any external `C` libraries, the python wrapper
requires the `Python.h` header (which should be included in a default python
installation) and `numpy` (including `numpy/arrayobject.h`).


## Usage

Below are two examples that illustrate how to use the ball tree from `C` and
`Python`.

### Using the `C` library

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

### Using the `Python` wrapper

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


## Comparison scipy.spatial.KDTree

The python package `scipy` implements a popular KDTree in
`scipy.spatial.KDTree`. The majority of this code is written in `Cython/C++`.

### Setup

- Dataset: `953,255` galaxies from the Baryon Oscillation Spectroscopic Survey,
  converted from sky coordinates *(right ascension, declination)* to points on the
  3D unit sphere *(x, y, z)*.
- Counting pairs formed between all objects within a fixed radius of `r <= 0.2`:
    - `balltree.count_radius(...)` (with unit weights)
    - `scipy.spatial.KDTree.query_ball_point(..., return_length=True)` (no weights)
- Counting the same pairs using the optimised dualtree algorithm.
    - `balltree.dualcount_radius(...)` (with unit weights)
    - `scipy.spatial.KDTree.count_neighbors(...)` (with unit weights)

### Results (single thread, AMD Epyc)

- Single point query using all points:
```
    balltree.count_radius:     found 24688969825 pairs in  26.737 sec
    KDTree.query_ball_point:   found 24688969825 pairs in 630.395 sec
```
- Using the dualtree algorithm:
```
    balltree.dualcount_radius: found 24688969825 pairs in  11.591 sec
    KDTree.count_neighbors:    found 24688969825 pairs in 321.993 sec
```

This corresponds to a **speed of of 25-30x** given test hardware and dataset.


## Maintainers

- Jan Luca van den Busch (author, Ruhr-Universität Bochum, Astronomisches Institut)
