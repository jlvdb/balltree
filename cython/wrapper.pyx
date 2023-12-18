import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free

from cython cimport cast
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
    int balltree_to_file(BallTree*, char*);
    BallTree* balltree_from_file(char*);


cdef class BallTreeWrapped:
    cdef BallTree *_tree

    def __cinit__(
        self,
        np.ndarray[np.double_t, ndim=1] x=None,
        np.ndarray[np.double_t, ndim=1] y=None,
        np.ndarray[np.double_t, ndim=1] z=None,
        np.ndarray[np.double_t, ndim=1] weight=None,
        int leafsize=20,
    ):
        size = x.shape[0]
        if weight is None:
            weight = np.ones_like(x)

        if size != y.shape[0] or size != z.shape[0] or size != weight.shape[0]:
            raise ValueError("input arrays must have the same length")

        cdef PointBuffer *point_buffer = pointbuffer_create(size)
        if point_buffer is NULL:
            raise MemoryError("Failed to allocate memory for PointBuffer")
        for i in range(size):
            point_buffer.points[i].x = x[i]
            point_buffer.points[i].y = y[i]
            point_buffer.points[i].z = z[i]
            point_buffer.points[i].weight = weight[i]

        self._tree = balltree_build(point_buffer, leafsize)
        pointbuffer_free(point_buffer)  # data is copied into the tree itself
        if self._tree is NULL:
            raise RuntimeError("failed to allocate and build tree structure")

    def __dealloc__(self):
        if self._tree is not NULL:
            balltree_free(self._tree)

    def count_radius(self, point: NDArray[np.double], radius: float, weight: float | None = None) -> float:
        x, y, z = point
        if weight is None:
            weight = 1.0
        cdef Point qpoint = Point(x, y, z, weight)
        return balltree_count_radius(self._tree, &qpoint, radius)

    def dualcount_radius(self, other: BallTreeWrapped, radius: float) -> float:
        return balltree_dualcount_radius(self._tree, other._tree, radius)

    def serialize(self, path: str):
        balltree_to_file(self._tree, path.encode("utf-8"))


cdef class BallTreeWrapped2:
    cdef BallTree *_tree

    def __cinit__(self, path: str):
        self._tree = balltree_from_file(path.encode("utf-8"))
        if self._tree is NULL:
            raise RuntimeError("failed to allocate and build tree structure")

    def __dealloc__(self):
        if self._tree is not NULL:
            balltree_free(self._tree)

    def count_radius(self, point: NDArray[np.double], radius: float, weight: float | None = None) -> float:
        x, y, z = point
        if weight is None:
            weight = 1.0
        cdef Point qpoint = Point(x, y, z, weight)
        return balltree_count_radius(self._tree, &qpoint, radius)

    def dualcount_radius(self, other: BallTreeWrapped, radius: float) -> float:
        return balltree_dualcount_radius(self._tree, other._tree, radius)

    def serialize(self, path: str):
        balltree_to_file(self._tree, path.encode("utf-8"))
