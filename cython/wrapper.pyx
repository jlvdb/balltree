import numpy as np
cimport numpy as np


cdef extern from "../include/point.h":
    cdef struct Point:
        double x
        double y
        double z
        double weight

    cdef enum Axis:
        X
        Y
        Z

    cdef struct PointBuffer:
        int size
        Point *points

    cdef struct PointSlice:
        int start
        int end
        Point *points

    Point create_point_weighted(double, double, double, double)
    Point create_point_unweighted(double, double, double)
    void print_point(const Point*)
    double points_distance(const Point*, const Point*)
    double points_distance2(const Point*, const Point*)

    PointBuffer* pointbuffer_create(int)
    void pointbuffer_free(PointBuffer*)
    int pointbuffer_resize(PointBuffer*, int)
    void print_pointbuffer(const PointBuffer*)
    double count_within_radius(const PointBuffer*, const Point*, double)
    double count_within_range(const PointBuffer*, const Point*, double, double)

    PointSlice* pointslice_from_buffer(const PointBuffer*)
    void pointslice_free(PointSlice*)
    void print_pointslice(const PointSlice*)
    int get_pointslice_size(const PointSlice*)


cdef extern from "../include/balltree.h":
    cdef struct BallTree:
        Point center
        double radius
        BallTree *left
        BallTree *right
        PointBuffer data

    int balltree_is_leaf(const BallTree*)
    void balltree_free(BallTree*)
    void balltree_print(const BallTree*)
    BallTree* balltree_build_recursive(PointSlice*, int)
    BallTree* balltree_build(PointBuffer*, int)
    double balltree_count_radius(BallTree*, Point*, double)
    double balltree_count_range(BallTree*, Point*, double, double)


cdef class TestClass:
    cdef Point point

    def __cinit__(self, x: double, y: double, z: double):
        self.point = Point(x, y, z, 1.0)
        print_point(&self.point)
