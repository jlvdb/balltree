#include <Python.h>
#include <numpy/arrayobject.h>

#include <stdio.h>

#include "point.h"
#include "balltree.h"
#include "balltree_macros.h"

// helper functions
static const char *PyString_to_char(PyObject* py_string);
static Point *PyIter_to_point(PyObject *point_iter, double weight);
static NpyIter *numpy_iter_point_arrays(PyArrayObject *x_obj, PyArrayObject *y_obj, PyArrayObject *z_obj, PyArrayObject *weight_obj);
static npy_intp check_arrays_dtype_shape(PyObject *x_obj, PyObject *y_obj, PyObject *z_obj, PyObject *weight_obj);
static PyObject *numpy_ones_1dim(npy_intp size);
static PointBuffer *ptbuf_from_numpy_array(PyArrayObject *x_obj, PyArrayObject *y_obj, PyArrayObject *z_obj, PyArrayObject *weight_obj);
static PyObject *ptbuf_get_numpy_view(PointBuffer *buffer);
static PyObject *statvec_get_numpy_array(StatsVector *vec);

typedef struct {
    PyObject_HEAD
    BallTree *balltree;
} PyBallTree;
static PyTypeObject PyBallTreeType;

// PyBallTree: constructors & deallocators
static PyObject *PyBallTree_from_data(PyTypeObject *type, PyObject *args, PyObject *kwds);
static PyObject *PyBallTree_from_random(PyTypeObject *type, PyObject *args, PyObject *kwds);
static PyObject *PyBallTree_from_file(PyTypeObject *type, PyObject *args);
static PyObject *PyBallTree_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
static void PyBallTree_dealloc(PyObject *self);
// PyBallTree: property implementation
static PyObject *PyBallTree_get_data(PyObject *self, void *closure);
static PyObject *PyBallTree_get_num_data(PyObject *self, void *closure);
static PyObject *PyBallTree_get_leafsize(PyObject *self, void *closure);
static PyObject *PyBallTree_get_sum_weight(PyObject *self, void *closure);
static PyObject *PyBallTree_get_center(PyObject *self, void *closure);
static PyObject *PyBallTree_get_radius(PyObject *self, void *closure);
// PyBallTree: method implementations
static PyObject *PyBallTree_str(PyObject *self);
static PyObject *PyBallTree_to_file(PyObject *self, PyObject *args);
static PyObject *PyBallTree_count_nodes(PyObject *self);
static PyObject *PyBallTree_get_node_data(PyObject *self);
static PyObject *PyBallTree_count_radius(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *PyBallTree_count_range(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *PyBallTree_dualcount_radius(PyObject *self, PyObject *args);
static PyObject *PyBallTree_dualcount_range(PyObject *self, PyObject *args);

// helper functions ////////////////////////////////////////////////////////////

static const char *PyString_to_char(PyObject* py_string) {
    if (!PyUnicode_Check(py_string)) {
        PyErr_SetString(PyExc_TypeError, "input must be a string");
        return NULL;
    }

    // Convert Python string to UTF-8 C string
    const char *string = PyUnicode_AsUTF8(py_string);
    if (string == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "failed to convert string to UTF-8");
        return NULL;
    }
    return string;
}

static Point *PyIter_to_point(PyObject *coord_iter, double weight) {
    // check the input to be an iterator with size 3
    const char err_msg[] = "expected a sequence with length 3";
    PyObject *coord_seq = PySequence_Fast(coord_iter, err_msg);
    if (coord_seq == NULL) {
        goto fail;
    }
    Py_ssize_t seq_length = PySequence_Fast_GET_SIZE(coord_seq);
    if (seq_length != 3) {
        PyErr_SetString(PyExc_ValueError, err_msg);
        goto fail;
    }

    // create the point instance from the python sequence
    Point *point = (Point *)malloc(sizeof(Point));
    if (point == NULL) {
        PyErr_SetString(PyExc_MemoryError, "failed to allocate Point");
        goto fail;
    }
    point->x = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(coord_seq, 0));
    point->y = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(coord_seq, 1));
    point->z = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(coord_seq, 2));
    point->weight = weight;
    Py_DECREF(coord_seq);
    if (PyErr_Occurred()) {
        PyErr_SetString(PyExc_ValueError, "cannot convert coordinates to float");
        return NULL;
    }

    return point;

fail:
    Py_DECREF(coord_seq);
    return NULL;
}

static NpyIter *numpy_iter_point_arrays(
    PyArrayObject *x_obj,
    PyArrayObject *y_obj,
    PyArrayObject *z_obj,
    PyArrayObject *weight_obj
) {
    const int n_arrays = 4;
    PyArrayObject *arrays[4] = {x_obj, y_obj, z_obj, weight_obj};
    npy_uint32 array_flags[4] = {
        NPY_ITER_READONLY,
        NPY_ITER_READONLY,
        NPY_ITER_READONLY,
        NPY_ITER_READONLY,
    };
    NpyIter *iter = NpyIter_MultiNew(
        n_arrays,
        arrays,
        NPY_ITER_EXTERNAL_LOOP,
        NPY_KEEPORDER,
        NPY_SAFE_CASTING,
        array_flags,
        NULL
    );
    if (iter == NULL) {
        PyErr_SetString(PyExc_IndexError, "failed to create interator for input arrays");
        return NULL;
    }
    return iter;
}

static npy_intp check_arrays_dtype_shape(
    PyObject *x_obj,
    PyObject *y_obj,
    PyObject *z_obj,
    PyObject *weight_obj
) {
    // check if all inputs are arrays
    int weights_provided = weight_obj != Py_None;
    if (PyArray_Check(x_obj) == 0) {
        PyErr_SetString(PyExc_TypeError, "'x' must be a numpy array");
        return -1;
    }
    if (PyArray_Check(y_obj) == 0) {
        PyErr_SetString(PyExc_TypeError, "'y' must be a numpy array");
        return -1;
    }
    if (PyArray_Check(z_obj) == 0) {
        PyErr_SetString(PyExc_TypeError, "'z' must be a numpy array");
        return -1;
    }
    if (weights_provided && PyArray_Check(weight_obj) == 0) {
        PyErr_SetString(PyExc_TypeError, "'weight' must be a numpy array or None");
        return -1;
    }

    // check input array dimensions
    if (PyArray_NDIM(x_obj) != 1 ||
        PyArray_NDIM(y_obj) != 1 ||
        PyArray_NDIM(z_obj) != 1 ||
        (weights_provided && PyArray_NDIM(weight_obj) != 1)
    ) {
        PyErr_SetString(PyExc_ValueError, "all arrays must be one-dimensional");
        return -1;
    }
    npy_intp size = PyArray_DIM(x_obj, 0);
    if (PyArray_DIM(y_obj, 0) != size ||
        PyArray_DIM(z_obj, 0) != size ||
        (weights_provided && PyArray_DIM(weight_obj, 0) != size)
    ) {
        PyErr_SetString(PyExc_ValueError, "all arrays must have the same shape");
        return -1;
    }

    return size;
}

static PyObject *numpy_ones_1dim(npy_intp size) {
    const npy_intp ndim = 1;
    npy_intp shape[1] = {size};
    PyObject *array = PyArray_EMPTY(ndim, shape, NPY_DOUBLE, 0);
    if (array == NULL) {
        PyErr_SetString(PyExc_MemoryError, "failed to allocate array");
        return NULL;
    }

    npy_double *buffer = PyArray_DATA(array);
    for (npy_intp i = 0; i < size; ++i) {
        buffer[i] = 1.0;
    }
    return array;
}

static PointBuffer *ptbuf_from_numpy_array(
    PyArrayObject *x_obj,
    PyArrayObject *y_obj,
    PyArrayObject *z_obj,
    PyArrayObject *weight_obj
) {
    NpyIter *multi_iter = numpy_iter_point_arrays(x_obj, y_obj, z_obj, weight_obj);
    if (multi_iter == NULL) {
        return NULL;
    }

    // set up the output buffer
    int size = NpyIter_GetIterSize(multi_iter);
    PointBuffer *buffer = ptbuf_new(size);
    if (buffer == NULL) {
        NpyIter_Deallocate(multi_iter);
        return NULL;
    }
    Point *points = buffer->points;

    // set up pointers for iteration, no require any defereferencing
    NpyIter_IterNextFunc *iternext = NpyIter_GetIterNext(multi_iter, NULL);
    if (iternext == NULL) {
        PyErr_SetString(PyExc_ValueError, "failed to create iterator for input arrays");
        NpyIter_Deallocate(multi_iter);
        ptbuf_free(buffer);
        return NULL;
    }
    char **dataptr = NpyIter_GetDataPtrArray(multi_iter);
    npy_intp *strideptr = NpyIter_GetInnerStrideArray(multi_iter);
    npy_intp *innersizeptr = NpyIter_GetInnerLoopSizePtr(multi_iter);

    // iterate and interpret values as doubles and copy values into PointBuffer
    int idx_ptbuf = 0;
    do {  // the outer loop should be just a single iteration
        char *ptr_x = dataptr[0];
        char *ptr_y = dataptr[1];
        char *ptr_z = dataptr[2];
        char *ptr_weight = dataptr[3];
        npy_intp stride_x = strideptr[0];
        npy_intp stride_y = strideptr[1];
        npy_intp stride_z = strideptr[2];
        npy_intp stride_weight = strideptr[3];

        // run inner loop of all elements
        npy_intp count = *innersizeptr;
        while (count--) {
            points[idx_ptbuf] = point_create_weighted(
                *(double *)ptr_x,
                *(double *)ptr_y,
                *(double *)ptr_z,
                *(double *)ptr_weight
            );
            ++idx_ptbuf;
            // move the iterator corresponding to the data size
            ptr_x += stride_x;
            ptr_y += stride_y;
            ptr_z += stride_z;
            ptr_weight += stride_weight;

        }
    } while (iternext(multi_iter));

    NpyIter_Deallocate(multi_iter);
    return buffer;
}

static PyObject *ptbuf_get_numpy_view(PointBuffer *buffer) {
    const npy_intp ndim = 1;
    npy_intp shape[1] = {buffer->size};

    // construct an appropriate dtype for Point
    PyObject *arr_dtype = Py_BuildValue(
        "[(ss)(ss)(ss)(ss)]",
        "x", "f8",
        "y", "f8",
        "z", "f8",
        "weight", "f8"
    );
    if (arr_dtype == NULL) {
        return NULL;
    
}
    // get the numpy API array descriptor
    PyArray_Descr *arr_descr;
    int result = PyArray_DescrConverter(arr_dtype, &arr_descr);  // PyArray_Descr **
    Py_DECREF(arr_dtype);
    if (result != NPY_SUCCEED) {
        return NULL;
    }

    // get an array view
    return PyArray_NewFromDescr(
        &PyArray_Type,
        arr_descr,  // reference stolen, no Py_DECREF needed
        ndim,
        shape,
        NULL,
        buffer->points,
        NPY_ARRAY_CARRAY,
        NULL
    );
}

static PyObject *statvec_get_numpy_array(StatsVector *vec) {
    const npy_intp ndim = 1;
    npy_intp shape[1] = {vec->end};

    // construct an appropriate dtype for StatsVector
    PyObject *arr_dtype = Py_BuildValue(
        "[(ss)(ss)(ss)(ss)(ss)(ss)(ss)]",
        "depth", "i4",
        "n_points", "i4",
        "sum_weight", "f8",
        "x", "f8",
        "y", "f8",
        "z", "f8",
        "radius", "f8"
    );
    if (arr_dtype == NULL) {
        return NULL;
    }

    // get the numpy API array descriptor
    PyArray_Descr *arr_descr;
    int result = PyArray_DescrConverter(arr_dtype, &arr_descr);  // PyArray_Descr **
    Py_DECREF(arr_dtype);
    if (result != NPY_SUCCEED) {
        return NULL;
    }

    // create an uninitialised array and copy the data into it
    PyObject *array = PyArray_Empty(ndim, shape, arr_descr, 0);
    Py_DECREF(arr_descr);
    if (array == NULL) {
        return NULL;
    }
    void *ptr = PyArray_DATA(array);
    memcpy(ptr, vec->stats, sizeof(NodeStats) * vec->end);
    return array;
}

// PyBallTree: constructors & deallocators /////////////////////////////////////

static PyObject *PyBallTree_from_data(
    PyTypeObject *type,
    PyObject *args,
    PyObject *kwargs
) {
    static char *kwlist[] = {"x", "y", "z", "weight", "leafsize", NULL};
    PyObject *x_obj, *y_obj, *z_obj, *weight_obj = Py_None;
    int leafsize = DEFAULT_LEAFSIZE;
    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, "OOO|Oi", kwlist,
            &x_obj, &y_obj, &z_obj, &weight_obj, &leafsize)
    ) {
        return NULL;
    }
    int weights_provided = weight_obj != Py_None;

    // check the input data and create a numpy array with ones if weight is None
    npy_intp size = check_arrays_dtype_shape(x_obj, y_obj, z_obj, weight_obj);
    if (size == -1) {
        return NULL;
    }
    if (size == 0) {
        PyErr_SetString(PyExc_ValueError, "need at least one data point");
        return NULL;
    }
    if (!weights_provided) {
        weight_obj = numpy_ones_1dim(size);
    }

    PointBuffer *buffer = ptbuf_from_numpy_array(
        (PyArrayObject *)x_obj,
        (PyArrayObject *)y_obj,
        (PyArrayObject *)z_obj,
        (PyArrayObject *)weight_obj
    );
    if (!weights_provided) {
        Py_DECREF(&weight_obj);
    }
    if (buffer == NULL) {
        return NULL;
    }

    // build the balltree
    BallTree *tree = balltree_build_nocopy(buffer, leafsize);  // buffer owner
    if (tree == NULL) {
        return NULL;
    }
    PyBallTree *self = (PyBallTree *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->balltree = tree;
    }
    return (PyObject *)self;
}

PyDoc_STRVAR(
    from_random__doc__,
    "from_random(low, high, size, leafsize=None)\n"
    "--\n\n"  // Required + Python convention
    "Build a ball tree from randomly generated points.\n"
    "..."
);

static PyObject *PyBallTree_from_random(
    PyTypeObject *type,
    PyObject *args,
    PyObject *kwargs
) {
    static char *kwlist[] = {"low", "high", "size", "leafsize", NULL};
    double low, high;
    int size;
    int leafsize = DEFAULT_LEAFSIZE;
    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, "ddi|i", kwlist,
            &low, &high, &size, &leafsize)
    ) {
        return NULL;
    }

    // generate random data points
    PointBuffer *buffer = ptbuf_gen_random(low, high, size);
    if (buffer == NULL) {
        return NULL;
    }

    // build the balltree
    BallTree *tree = balltree_build_leafsize(buffer, leafsize);
    if (tree == NULL) {
        return NULL;
    }
    ptbuf_free(buffer);  // buffer is copied into tree
    PyBallTree *self = (PyBallTree *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->balltree = tree;
    }
    return (PyObject *)self;
}

static PyObject *PyBallTree_from_file(PyTypeObject *type, PyObject *args) {
    PyObject *py_string;
    if (!PyArg_ParseTuple(args, "O!", &PyUnicode_Type, &py_string)) {
        return NULL;
    }
    const char *path = PyString_to_char(py_string);  // buffer is managed by py_string
    if (path == NULL) {
        return NULL;
    }
    // build the balltree
    BallTree *tree = balltree_from_file(path);
    if (tree == NULL) {
        return NULL;
    }
    PyBallTree *self = (PyBallTree *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->balltree = tree;
    }
    return (PyObject *)self;
}

static PyObject *PyBallTree_new(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
    return PyBallTree_from_data(type, args, kwargs);
}

static void PyBallTree_dealloc(PyObject *self) {
    PyBallTree *pytree = (PyBallTree *)self;
    if (pytree->balltree != NULL) {
        balltree_free(pytree->balltree);
    }
    Py_TYPE(self)->tp_free(self);
}

// PyBallTree: property implementations ////////////////////////////////////////

static PyObject *PyBallTree_get_data(PyObject *self, void *closure) {
    PyBallTree *pytree = (PyBallTree *)self;
    return ptbuf_get_numpy_view(&pytree->balltree->data);
}

static PyObject *PyBallTree_get_num_data(PyObject *self, void *closure) {
    PyBallTree *pytree = (PyBallTree *)self;
    return PyLong_FromLong(pytree->balltree->data.size);
}

static PyObject *PyBallTree_get_leafsize(PyObject *self, void *closure) {
    PyBallTree *pytree = (PyBallTree *)self;
    return PyLong_FromLong(pytree->balltree->leafsize);
}

static PyObject *PyBallTree_get_sum_weight(PyObject *self, void *closure) {
    PyBallTree *pytree = (PyBallTree *)self;
    return PyFloat_FromDouble(pytree->balltree->root->sum_weight);
}

static PyObject *PyBallTree_get_center(PyObject *self, void *closure) {
    PyBallTree *pytree = (PyBallTree *)self;
    return PyTuple_Pack(
        3,
        PyFloat_FromDouble(pytree->balltree->root->center.x),
        PyFloat_FromDouble(pytree->balltree->root->center.y),
        PyFloat_FromDouble(pytree->balltree->root->center.z)
    );
}

static PyObject *PyBallTree_get_radius(PyObject *self, void *closure) {
    PyBallTree *pytree = (PyBallTree *)self;
    return PyFloat_FromDouble(pytree->balltree->root->radius);
}

// PyBallTree: method implementations //////////////////////////////////////////

static PyObject *PyBallTree_str(PyObject *self) {
    PyBallTree *pytree = (PyBallTree *)self;
    BallTree *tree = pytree->balltree;
    BallNode *node = tree->root;

    // constuct the string representation
    char buffer[256];
    int n_bytes = snprintf(
        buffer,
        sizeof(buffer),
        "PyBallTree(num_data=%d, radius=%.3f, center={%+.3f, %+.3f, %+.3f})",
        tree->data.size,
        node->radius,
        node->center.x,
        node->center.y,
        node->center.z
    );
    if (n_bytes < 0 || n_bytes >= (int)sizeof(buffer)) {
        PyErr_SetString(PyExc_RuntimeError, "failed to format the string");
        return NULL;
    }

    // convert to python str
    PyObject *py_string = PyUnicode_DecodeUTF8(buffer, n_bytes, "strict");
    if (py_string == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "failed to convert to python unicode");
        return NULL;
    }
    return py_string;
}

PyDoc_STRVAR(
    to_file__doc__,
    "to_file(path)\n"
    "--\n\n"  // Required + Python convention
    "Serialise the tree to a binary file.\n"
    "..."
);

static PyObject *PyBallTree_to_file(PyObject *self, PyObject *args) {
    PyObject *py_string;
    if (!PyArg_ParseTuple(args, "O!", &PyUnicode_Type, &py_string)) {
        return NULL;
    }
    const char *path = PyString_to_char(py_string);  // buffer is managed by py_string
    if (path == NULL) {
        return NULL;
    }

    PyBallTree *pytree = (PyBallTree *)self;
    if (balltree_to_file(pytree->balltree, path) == BTR_FAILED) {
        return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject *PyBallTree_count_nodes(PyObject *self) {
    PyBallTree *pytree = (PyBallTree *)self;
    int count = balltree_count_nodes(pytree->balltree);
    return PyLong_FromLong(count);
}

static PyObject *PyBallTree_get_node_data(PyObject *self) {
    PyBallTree *pytree = (PyBallTree *)self;
    StatsVector *vec = balltree_collect_stats(pytree->balltree);
    if (vec == NULL) {
        return NULL;
    }

    PyObject *array = statvec_get_numpy_array(vec);
    statvec_free(vec);
    return array;
}

static PyObject *PyBallTree_count_radius(
    PyObject *self,
    PyObject *args,
    PyObject *kwargs
) {
    static char *kwlist[] = {"point", "radius", "weight", NULL};
    PyObject *coord_iter;  // TODO: refcount?
    double radius;
    double weight = 1.0;  // weight is optional
    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, "Od|d", kwlist,
            &coord_iter, &radius, &weight)
    ) {
        return NULL;
    }
    Point *point = PyIter_to_point(coord_iter, weight);
    if (point == NULL) {
        return NULL;
    }

    PyBallTree *pytree = (PyBallTree *)self;
    double count = balltree_count_radius(pytree->balltree, point, radius);
    free(point);
    return PyFloat_FromDouble(count);
}

static PyObject *PyBallTree_count_range(
    PyObject *self,
    PyObject *args,
    PyObject *kwargs
) {
    static char *kwlist[] = {"point", "rmin", "rmax", "weight", NULL};
    PyObject *coord_iter;  // TODO: refcount?
    double rmin, rmax;
    double weight = 1.0;  // weight is optional
    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, "Odd|d", kwlist,
            &coord_iter, &rmin, &rmax, &weight)
    ) {
        return NULL;
    }
    Point *point = PyIter_to_point(coord_iter, weight);
    if (point == NULL) {
        return NULL;
    }

    PyBallTree *pytree = (PyBallTree *)self;
    double count = balltree_count_range(pytree->balltree, point, rmin, rmax);
    free(point);
    return PyFloat_FromDouble(count);
}

static PyObject *PyBallTree_dualcount_radius(PyObject *self, PyObject *args) {
    PyBallTree *other_tree;
    double radius;
    if (!PyArg_ParseTuple(args, "O!d", &PyBallTreeType, &other_tree, &radius)) {
        return NULL;
    }

    PyBallTree *pytree = (PyBallTree *)self;
    double count = balltree_dualcount_radius(
        pytree->balltree,
        other_tree->balltree,
        radius
    );
    return PyFloat_FromDouble(count);
}

static PyObject *PyBallTree_dualcount_range(PyObject *self, PyObject *args) {
    PyBallTree *other_tree;
    double rmin, rmax;
    if (!PyArg_ParseTuple(args, "O!dd", &PyBallTreeType, &other_tree, &rmin, &rmax)) {
        return NULL;
    }

    PyBallTree *pytree = (PyBallTree *)self;
    double count = balltree_dualcount_range(
        pytree->balltree,
        other_tree->balltree,
        rmin,
        rmax
    );
    return PyFloat_FromDouble(count);
}

// PyBallTree definition ///////////////////////////////////////////////////////

static PyGetSetDef PyBallTree_getset[] = {
    {
        "data",
        PyBallTree_get_data,
        NULL,
        "get a view of the underlying data",
        NULL
    },
    {
        "num_data",
        PyBallTree_get_num_data,
        NULL,
        "get the number of data points stored in the tree",
        NULL
    },
    {
        "leafsize",
        PyBallTree_get_leafsize,
        NULL,
        "get the leafsize of the tree",
        NULL
    },
    {
        "sum_weight",
        PyBallTree_get_sum_weight,
        NULL,
        "get the sum of points weights of the tree",
        NULL
    },
    {
        "center",
        PyBallTree_get_center,
        NULL,
        "get the center point associated with the root node",
        NULL
    },
    {
        "radius",
        PyBallTree_get_radius,
        NULL,
        "get radius associated with the root node",
        NULL
    },
    {NULL, NULL, NULL, NULL, NULL}
};

static PyMethodDef PyBallTree_methods[] = {
    // constructors
    {
        "from_random",
        (PyCFunction)PyBallTree_from_random,
        METH_CLASS | METH_VARARGS | METH_KEYWORDS,
        from_random__doc__
    },
    {
        "from_file",
        (PyCFunction)PyBallTree_from_file,
        METH_CLASS | METH_VARARGS,
        "Deserialize a PyBallTree instance from a file"
    },
    // regular methods
    {
        "to_file",
        (PyCFunction)PyBallTree_to_file,
        METH_VARARGS,
        to_file__doc__
    },
    {
        "count_nodes",
        (PyCFunction)PyBallTree_count_nodes,
        METH_NOARGS,
        "Count nodes contained in tree"
    },
    {
        "get_node_data",
        (PyCFunction)PyBallTree_get_node_data,
        METH_NOARGS,
        "Get the primary node information"
    },
    {
        "count_radius",
        (PyCFunction)PyBallTree_count_radius,
        METH_VARARGS | METH_KEYWORDS,
        "Count pairs within a radius"
    },
    {
        "count_range",
        (PyCFunction)PyBallTree_count_range,
        METH_VARARGS | METH_KEYWORDS,
        "Count pairs within a range"
    },
    {
        "dualcount_radius",
        (PyCFunction)PyBallTree_dualcount_radius,
        METH_VARARGS,
        "Count pairs within a radius with the dualtree algorithm"
    },
    {
        "dualcount_range",
        (PyCFunction)PyBallTree_dualcount_range,
        METH_VARARGS,
        "Count pairs within a range with the dualtree algorithm"
    },
    {NULL, NULL, 0, NULL}
};

static PyTypeObject PyBallTreeType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "balltree.BallTree",
    .tp_doc = "Python wrapper for C BallTree",
    .tp_basicsize = sizeof(PyBallTree),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = PyBallTree_new,
    .tp_dealloc = PyBallTree_dealloc,
    .tp_methods = PyBallTree_methods,
    .tp_getset = PyBallTree_getset,
    .tp_str = PyBallTree_str,
};

// module configuration ////////////////////////////////////////////////////////

static struct PyModuleDef pyballtree = {
    PyModuleDef_HEAD_INIT,
    .m_name = "balltree",
    .m_doc = "Fast balltree implementation",
    .m_size = -1,
};

PyMODINIT_FUNC PyInit_balltree(void) {
    // Import NumPy API and check for errors
    import_array();
    if (PyErr_Occurred()) {
        PyErr_SetString(PyExc_ImportError, "failed to import NumPy array module");
        return NULL;
    }

    PyObject *module = PyModule_Create(&pyballtree);
    if (module == NULL) {
        return NULL;
    }

    if (PyType_Ready(&PyBallTreeType) == -1) {
        return NULL;
    }
    Py_INCREF(&PyBallTreeType);
    if (PyModule_AddObject(module, "BallTree", (PyObject *)&PyBallTreeType) == -1) {
        Py_DECREF(&PyBallTreeType);
        return NULL;
    }

    return module;
}
