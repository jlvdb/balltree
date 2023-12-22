import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free

from numpy.typing import NDArray


cdef extern from "../include/point.h":
    cdef struct Point:
        double x
        double y
        double z
        double weight

    cdef struct PointBuffer:
        int size
        Point *points

    cdef struct PointSlice:
        int start
        int end
        Point *points

    Point point_create(double, double, double)
    Point point_create_weighted(double, double, double, double);
    PointBuffer* pointbuffer_create(int)
    void pointbuffer_free(PointBuffer*)
    int pointbuffer_resize(PointBuffer*, int)


cdef extern from "../include/ballnode.h":
    cdef struct BallNode:
        Point center
        double radius
        BallNode *left
        BallNode *right
        PointSlice data


cdef extern from "../include/balltree.h":
    cdef struct BallTree:
        BallNode *root
        PointBuffer data
        int leafsize

    BallTree* balltree_build(PointBuffer*);
    BallTree* balltree_build_leafsize(PointBuffer*, int);
    void balltree_free(BallTree*);
    int balltree_count_nodes(BallTree*);
    int balltree_to_file(BallTree*, char*);
    BallTree* balltree_from_file(char*);

    double balltree_count_radius(BallTree*, Point*, double);
    double balltree_count_range(BallTree*, Point*, double, double);
    double balltree_dualcount_radius(BallTree*, BallTree*, double);


def BallTree_from_data(
    double[:] x,
    double[:] y,
    double[:] z,
    double[:] weight = None,
    int leafsize = -1,
) -> BallTreeWrapped:
    size = x.shape[0]
    if weight is None:
        weight = np.ones_like(x)

    if size != y.shape[0] or size != z.shape[0] or size != weight.shape[0]:
        raise ValueError("input arrays must have the same length")

    cdef PointBuffer *point_buffer = pointbuffer_create(size)
    if point_buffer is NULL:
        raise MemoryError
    for i in range(size):
        point_buffer.points[i].x = x[i]
        point_buffer.points[i].y = y[i]
        point_buffer.points[i].z = z[i]
        point_buffer.points[i].weight = weight[i]

    if leafsize == -1:
        tree = balltree_build(point_buffer)
    else:
        tree = balltree_build_leafsize(point_buffer, leafsize)
    pointbuffer_free(point_buffer)  # data is copied into the tree itself
    if tree is NULL:
        raise RuntimeError("failed to allocate or build tree structure")

    cdef BallTreeWrapped wrapped = BallTreeWrapped.__new__(BallTreeWrapped)
    wrapped._tree = tree
    return wrapped


def BallTree_from_file(str path) -> BallTreeWrapped:
    cdef BallTreeWrapped wrapped = BallTreeWrapped.__new__(BallTreeWrapped)
    wrapped._tree = balltree_from_file(path.encode("utf-8"))
    return wrapped


cdef class BallTreeWrapped:
    cdef BallTree *_tree

    def __dealloc__(self):
        if self._tree is not NULL:
            balltree_free(self._tree)
            self._tree = NULL

    def __init__(self):
        # Prevent accidental instantiation from normal Python code
        # since we cannot pass a struct pointer into a Python constructor.
        raise TypeError("class cannot be instantiated directly, use .from_*() methods")

    from_data = staticmethod(BallTree_from_data)

    from_file = staticmethod(BallTree_from_file)

    def to_file(self, path: str):
        balltree_to_file(self._tree, path.encode("utf-8"))

    def count_radius(self, point: NDArray[np.double], radius: float, weight: float | None = None) -> float:
        x, y, z = point
        if weight is None:
            weight = 1.0
        cdef Point qpoint = Point(x, y, z, weight)
        return balltree_count_radius(self._tree, &qpoint, radius)

    def dualcount_radius(self, other: BallTreeWrapped, radius: float) -> float:
        return balltree_dualcount_radius(self._tree, other._tree, radius)
