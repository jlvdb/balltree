import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free


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

    PointBuffer* pointbuffer_create(int)
    void pointbuffer_free(PointBuffer*)
    int pointbuffer_resize(PointBuffer*, int)


cdef extern from "../include/balltree.h":
    cdef struct BallNode:
        Point center
        double radius
        BallNode *left
        BallNode *right
        PointSlice data

    cdef struct BallTree:
        BallNode *root
        PointBuffer data
        int leafsize

    void balltree_free(BallTree*)
    BallTree* balltree_build(PointBuffer*, int)
    double balltree_count_radius(BallTree*, Point*, double)
    double balltree_count_range(BallTree*, Point*, double, double)
    double balltree_dualcount_radius(BallTree*, BallTree*, double)


cdef class BallTreeWrapped:
    cdef BallTree *_tree

    def __cinit__(
        self,
        np.ndarray[np.double_t, ndim=1] x,
        np.ndarray[np.double_t, ndim=1] y,
        np.ndarray[np.double_t, ndim=1] z,
        np.ndarray[np.double_t, ndim=1] weight=None,
        int leafsize=20,
    ):
        size = x.shape[0]
        if weight is None:
            # Create a weight array with double ones if not provided
            weight = np.ones_like(x)

        # Assuming x, y, and z have the same length and are 1-dimensional arrays
        if size != y.shape[0] or size != z.shape[0] or size != weight.shape[0]:
            raise ValueError("Input arrays must have the same length")

        # Create a PointBuffer and copy data
        cdef PointBuffer *point_buffer = pointbuffer_create(size)
        if point_buffer is NULL:
            raise MemoryError("Failed to allocate memory for PointBuffer")

        for i in range(point_buffer.size):
            point_buffer.points[i].x = x[i]
            point_buffer.points[i].y = y[i]
            point_buffer.points[i].z = z[i]
            point_buffer.points[i].weight = weight[i]

        # Build the BallTree
        self._tree = balltree_build(point_buffer, leafsize)
        pointbuffer_free(point_buffer)

    def __dealloc__(self):
        balltree_free(self._tree)

    def count_radius(self, x: float, y: float, z: float, weight: float, radius: float) -> float:
        cdef Point point = Point(x, y, z, weight)
        return balltree_count_radius(self._tree, &point, radius)

    def dualcount_radius(self, other: BallTreeWrapped, radius: float) -> float:
        return balltree_dualcount_radius(self._tree, other._tree, radius)
