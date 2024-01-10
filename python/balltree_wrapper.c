#include <Python.h>
#include <numpy/arrayobject.h>

#include <stdio.h>

#include "point.h"
#include "balltree.h"
#define SET_PYERR_STRING
#include "balltree_macros.h"

struct PyBallTree;

static const char *PyString_to_char(PyObject* py_string);
static Point *PyIter_to_point(PyObject *point_iter, double weight);
static PyObject *ptbuf_get_numpy_view(PointBuffer *buffer);
static PointBuffer *ptbuf_from_numpy_array(PyArrayObject *x_obj, PyArrayObject *y_obj, PyArrayObject *z_obj, PyArrayObject *weight_obj);
static PyObject *statvec_get_numpy_array(StatsVector *vec);

// PyBallTree: constructors & deallocators
static PyObject *PyBallTree_from_data(PyTypeObject *type, PyObject *args, PyObject *kwds);
static PyObject *PyBallTree_from_random(PyTypeObject *type, PyObject *args, PyObject *kwds);
static PyObject *PyBallTree_from_file(PyTypeObject *type, PyObject *args);
static PyObject *pyballtree_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
static void pyballtree_dealloc(PyObject *self);
// PyBallTree: property implementation
static PyObject *pyballtree_get_data(PyObject *self, void *closure);
static PyObject *pyballtree_get_num_data(PyObject *self, void *closure);
static PyObject *pyballtree_get_leafsize(PyObject *self, void *closure);
static PyObject *pyballtree_get_sum_weight(PyObject *self, void *closure);
static PyObject *pyballtree_get_center(PyObject *self, void *closure);
static PyObject *pyballtree_get_radius(PyObject *self, void *closure);
// PyBallTree: method implementations
static PyObject* pyballtree_str(PyObject *self);
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

static Point *PyIter_to_point(PyObject *point_iter, double weight) {
    // Convert the iterable to a fast sequence object
    PyObject *point_sequence = PySequence_Fast(
        point_iter,
        "expected a sequence with length 3"
    );
    if (point_sequence == NULL) {
        return NULL;
    }

    Py_ssize_t seq_length = PySequence_Fast_GET_SIZE(point_sequence);
    if (seq_length != 3) {
        PyErr_SetString(PyExc_ValueError, "iterable must have length 3");
        Py_DECREF(point_sequence);
        return NULL;
    }

    // create the point instance
    Point *point = (Point *)malloc(sizeof(Point));
    if (point == NULL) {
        return NULL;
    }
    point->x = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(point_sequence, 0));
    point->y = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(point_sequence, 1));
    point->z = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(point_sequence, 2));
    if (PyErr_Occurred()) {
        PyErr_SetString(PyExc_ValueError, "invalid float values in the iterable");
        Py_DECREF(point_sequence);
        return NULL;
    }
    point->weight = weight;

    Py_DECREF(point_sequence);
    return point;
}

static PointBuffer *ptbuf_from_numpy_array(
    PyArrayObject *x_obj,
    PyArrayObject *y_obj,
    PyArrayObject *z_obj,
    PyArrayObject *weight_obj
) {
    // set up iterator to iterate input arrays simultaneously
    const int nop = 4;
    PyArrayObject *op[nop] = {x_obj, y_obj, z_obj, weight_obj};
    npy_uint32 op_flags[nop] = {
        NPY_ITER_READONLY,
        NPY_ITER_READONLY,
        NPY_ITER_READONLY,
        NPY_ITER_READONLY,
    };
    NpyIter *iter = NpyIter_MultiNew(
        nop,
        op,
        NPY_ITER_EXTERNAL_LOOP,
        NPY_KEEPORDER,
        NPY_SAFE_CASTING,
        op_flags,
        NULL
    );
    if (iter == NULL) {
        return NULL;
    }

    // set up output buffer
    int size = NpyIter_GetIterSize(iter);
    PointBuffer *buffer = ptbuf_new(size);
    if (buffer == NULL) {
        return NULL;
    }
    Point *points = buffer->points;

    // set up pointers for iteration
    NpyIter_IterNextFunc *iternext = NpyIter_GetIterNext(iter, NULL);
    if (iternext == NULL) {
        NpyIter_Deallocate(iter);
        ptbuf_free(buffer);
        return NULL;
    }
    char **dataptr = NpyIter_GetDataPtrArray(iter);
    npy_intp *strideptr = NpyIter_GetInnerStrideArray(iter);
    npy_intp *innersizeptr = NpyIter_GetInnerLoopSizePtr(iter);

    int i = 0;
    do {
        char *ptr_x = dataptr[0];
        char *ptr_y = dataptr[1];
        char *ptr_z = dataptr[2];
        char *ptr_weight = dataptr[3];

        npy_intp stride_x = strideptr[0];
        npy_intp stride_y = strideptr[1];
        npy_intp stride_z = strideptr[2];
        npy_intp stride_weight = strideptr[3];
        npy_intp count = *innersizeptr;

        while (count--) {
            points[i] = point_create_weighted(
                *(double *)ptr_x,
                *(double *)ptr_y,
                *(double *)ptr_z,
                *(double *)ptr_weight
            );

            ptr_x += stride_x;
            ptr_y += stride_y;
            ptr_z += stride_z;
            ptr_weight += stride_weight;

            ++i;
        }
    } while (iternext(iter));

    NpyIter_Deallocate(iter);
    return buffer;
}

static PyObject *ptbuf_get_numpy_view(PointBuffer *buffer) {
    npy_intp dims[1] = {buffer->size};

    // construct an appropriate dtype for Point
    PyObject *arr_dtype = Py_BuildValue(
        "[(ss)(ss)(ss)(ss)]",
        "x", "f8",
        "y", "f8",
        "z", "f8",
        "weight", "f8"
    );
    if (arr_dtype == NULL)
        return NULL;

    // get an array view
    PyArray_Descr *arr_descr;
    int result = PyArray_DescrConverter(arr_dtype, &arr_descr);
    Py_DECREF(arr_dtype);
    if (result != NPY_SUCCEED)
        return NULL;
    return PyArray_NewFromDescr(
        &PyArray_Type,
        arr_descr,
        1,     // no. of dimensions
        dims,  // shape
        NULL,  // strides
        buffer->points,
        NPY_ARRAY_CARRAY,
        NULL
    );
}

static PyObject *statvec_get_numpy_array(StatsVector *vec) {
    // construct an appropriate dtype for Point
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
    if (arr_dtype == NULL)
        return NULL;

    // create a NumPy array containg a copy of the data
    npy_intp dims[1] = {vec->end};
    PyArray_Descr *arr_descr;
    int result = PyArray_DescrConverter(arr_dtype, &arr_descr);
    Py_DECREF(arr_dtype);
    if (result != NPY_SUCCEED)
        return NULL;
    
    PyObject *array = PyArray_Empty(1, dims, arr_descr, 0);
    if (array == NULL) {
        return NULL;
    }

    void *ptr = PyArray_DATA(array);
    memcpy(ptr, vec->stats, sizeof(NodeStats) * vec->end);
    return array;
}

// PyBallTree definition ///////////////////////////////////////////////////////

typedef struct {
    PyObject_HEAD
    BallTree *balltree;
} PyBallTree;

static PyMethodDef pyballtree_methods[] = {
    // constructors
    {
        "from_random",
        (PyCFunction)PyBallTree_from_random,
        METH_CLASS | METH_VARARGS | METH_KEYWORDS,
        "Build a PyBallTree instance"
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
        "Serialize a PyBallTree instance to a file"
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

static PyGetSetDef pyballtree_getset[] = {
    {
        "data",
        pyballtree_get_data,
        NULL,
        "get a view of the underlying data",
        NULL
    },
    {
        "num_data",
        pyballtree_get_num_data,
        NULL,
        "get the number of data points stored in the tree",
        NULL
    },
    {
        "leafsize",
        pyballtree_get_leafsize,
        NULL,
        "get the leafsize of the tree",
        NULL
    },
    {
        "sum_weight",
        pyballtree_get_sum_weight,
        NULL,
        "get the sum of points weights of the tree",
        NULL
    },
    {
        "center",
        pyballtree_get_center,
        NULL,
        "get the center point associated with the root node",
        NULL
    },
    {
        "radius",
        pyballtree_get_radius,
        NULL,
        "get radius associated with the root node",
        NULL
    },
    {NULL, NULL, NULL, NULL, NULL}
};

static PyTypeObject PyBallTreeType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "balltree.BallTree",
    .tp_doc = "Python wrapper for C BallTree",
    .tp_basicsize = sizeof(PyBallTree),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = pyballtree_new,
    .tp_dealloc = pyballtree_dealloc,
    .tp_methods = pyballtree_methods,
    .tp_getset = pyballtree_getset,
    .tp_str = pyballtree_str,
};

// PyBallTree: constructors & deallocators /////////////////////////////////////

static PyObject *PyBallTree_from_data(
    PyTypeObject *type,
    PyObject *args,
    PyObject *kwds
) {
    // todo: leafsize
    PyObject *x_obj, *y_obj, *z_obj, *weight_obj = Py_None;
    int leafsize = DEFAULT_LEAFSIZE;
    static char *kwlist[] = {"x", "y", "z", "weight", "leafsize", NULL};
    if (!PyArg_ParseTupleAndKeywords(
            args, kwds, "OOO|Oi", kwlist,
            &x_obj, &y_obj, &z_obj, &weight_obj, &leafsize)
    ) {
        return NULL;
    }
    int weights_provided = weight_obj != Py_None;

    // check if all inputs are arrays
    if (PyArray_Check(x_obj) == 0) {
        PyErr_SetString(PyExc_TypeError, "'x' must be a numpy array");
        return NULL;
    }
    if (PyArray_Check(y_obj) == 0) {
        PyErr_SetString(PyExc_TypeError, "'y' must be a numpy array");
        return NULL;
    }
    if (PyArray_Check(z_obj) == 0) {
        PyErr_SetString(PyExc_TypeError, "'z' must be a numpy array");
        return NULL;
    }
    if (weight_obj != Py_None && PyArray_Check(weight_obj) == 0) {
        PyErr_SetString(PyExc_TypeError, "'weight' must be a numpy array or None");
        return NULL;
    }

    // check input array dimensions
    if (
        PyArray_NDIM(x_obj) != 1 ||
        PyArray_NDIM(y_obj) != 1 ||
        PyArray_NDIM(z_obj) != 1 ||
        (weights_provided && PyArray_NDIM(weight_obj) != 1)
    ) {
        PyErr_SetString(PyExc_ValueError, "all arrays must be one-dimensional");
        return NULL;
    }

    // create a numpy array with ones if weight is None
    npy_intp shape[1] = {PyArray_DIM(x_obj, 0)};
    if (!weights_provided) {
        weight_obj = PyArray_EMPTY(1, shape, NPY_DOUBLE, 0);
        if (!weight_obj) {
            return NULL;
        }
        // Fill the array with ones
        npy_double *weight_data = PyArray_DATA(weight_obj);
        for (npy_intp i = 0; i < shape[0]; ++i) {
            weight_data[i] = 1.0;
        }
    }

    PointBuffer *buffer = ptbuf_from_numpy_array(
        (PyArrayObject *)x_obj,
        (PyArrayObject *)y_obj,
        (PyArrayObject *)z_obj,
        (PyArrayObject *)weight_obj
    );
    if (buffer == NULL) {
        return NULL;
    }

    // create the object
    PyBallTree *self;
    self = (PyBallTree *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->balltree = balltree_build_nocopy(buffer, leafsize);
        if (self->balltree == NULL) {
            Py_DECREF(self);
            return NULL;
        }
    }
    return (PyObject *)self;
}

static PyObject *PyBallTree_from_random(
    PyTypeObject *type,
    PyObject *args,
    PyObject *kwds
) {
    // parse the python arguments
    double low, high;
    int size;
    int leafsize = DEFAULT_LEAFSIZE;

    static char *kwlist[] = {"low", "high", "size", "leafsize", NULL};
    if (!PyArg_ParseTupleAndKeywords(
            args, kwds, "ddi|i", kwlist,
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
    BallTree *tree;
    tree = balltree_build_leafsize(buffer, leafsize);
    if (tree == NULL) {
        return NULL;
    }
    ptbuf_free(buffer);  // buffer is copied into tree

    // create the object
    PyBallTree *self;
    self = (PyBallTree *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->balltree = tree;
        if (self->balltree == NULL) {
            Py_DECREF(self);
            return NULL;
        }
    }
    return (PyObject *)self;
}

static PyObject *PyBallTree_from_file(PyTypeObject *type, PyObject *args) {
    // Parse arguments
    PyObject *py_string;
    if (!PyArg_ParseTuple(args, "O", &py_string)) {
        return NULL;
    }
    const char *path = PyString_to_char(py_string);
    if (path == NULL) {
        return NULL;
    }

    // create the object
    PyBallTree *self;
    self = (PyBallTree *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->balltree = balltree_from_file(path);
        if (self->balltree == NULL) {
            Py_DECREF(self);
            return NULL;
        }
    }
    return (PyObject *)self;
}

static PyObject *pyballtree_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    return PyBallTree_from_data(type, args, kwds);
}

static void pyballtree_dealloc(PyObject *self) {
    PyBallTree *pytree = (PyBallTree *)self;
    if (pytree->balltree != NULL) {
        balltree_free(pytree->balltree);
    }
    Py_TYPE(self)->tp_free(self);
}

// PyBallTree: property implementations ////////////////////////////////////////

static PyObject *pyballtree_get_data(PyObject *self, void *closure) {
    PyBallTree *pytree = (PyBallTree *)self;
    return ptbuf_get_numpy_view(&pytree->balltree->data);
}

static PyObject *pyballtree_get_num_data(PyObject *self, void *closure) {
    PyBallTree *pytree = (PyBallTree *)self;
    return PyLong_FromLong(pytree->balltree->data.size);
}

static PyObject *pyballtree_get_leafsize(PyObject *self, void *closure) {
    PyBallTree *pytree = (PyBallTree *)self;
    return PyLong_FromLong(pytree->balltree->leafsize);
}

static PyObject *pyballtree_get_sum_weight(PyObject *self, void *closure) {
    PyBallTree *pytree = (PyBallTree *)self;
    return PyFloat_FromDouble(pytree->balltree->root->sum_weight);
}

static PyObject *pyballtree_get_center(PyObject *self, void *closure) {
    PyBallTree *pytree = (PyBallTree *)self;
    return PyTuple_Pack(
        3,
        PyFloat_FromDouble(pytree->balltree->root->center.x),
        PyFloat_FromDouble(pytree->balltree->root->center.y),
        PyFloat_FromDouble(pytree->balltree->root->center.z)
    );
}

static PyObject *pyballtree_get_radius(PyObject *self, void *closure) {
    PyBallTree *pytree = (PyBallTree *)self;
    return PyFloat_FromDouble(pytree->balltree->root->radius);
}

// PyBallTree: method implementations //////////////////////////////////////////

static PyObject* pyballtree_str(PyObject *self) {
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
        PyErr_SetString(PyExc_RuntimeError, "failed to format the string.");
        return NULL;
    }

    // convert to python str
    PyObject *py_string = PyUnicode_DecodeUTF8(buffer, n_bytes, "strict");
    if (!py_string) {
        PyErr_SetString(PyExc_RuntimeError, "failed to convert the string to Python.");
        return NULL;
    }
    return py_string;
}

static PyObject *PyBallTree_to_file(PyObject *self, PyObject *args) {
    PyBallTree *pytree = (PyBallTree *)self;

    // Parse arguments
    PyObject *py_string;
    if (!PyArg_ParseTuple(args, "O", &py_string)) {
        return NULL;
    }
    const char *path = PyString_to_char(py_string);
    if (path == NULL) {
        return NULL;
    }

    balltree_to_file(pytree->balltree, path);
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
    PyBallTree *pytree = (PyBallTree *)self;

    // Parse the Python arguments
    PyObject *coord_iter;
    double radius;
    double weight = 1.0;  // weight is optional
    static char *kwlist[] = {"point", "radius", "weight", NULL};
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

    double count = balltree_count_radius(pytree->balltree, point, radius);
    free(point);
    return PyFloat_FromDouble(count);
}

static PyObject *PyBallTree_count_range(
    PyObject *self,
    PyObject *args,
    PyObject *kwargs
) {
    PyBallTree *pytree = (PyBallTree *)self;

    // Parse the Python arguments
    PyObject *coord_iter;
    double rmin, rmax;
    double weight = 1.0;  // weight is optional
    static char *kwlist[] = {"point", "rmin", "rmax", "weight", NULL};
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

    double count = balltree_count_range(pytree->balltree, point, rmin, rmax);
    free(point);
    return PyFloat_FromDouble(count);
}

static PyObject *PyBallTree_dualcount_radius(PyObject *self, PyObject *args) {
    PyBallTree *pytree = (PyBallTree *)self;

    // Parse the Python arguments
    PyBallTree *other_tree;
    double radius;
    if (!PyArg_ParseTuple(args, "O!d", &PyBallTreeType, &other_tree, &radius)) {
        return NULL;
    }

    double count = balltree_dualcount_radius(
        pytree->balltree,
        other_tree->balltree,
        radius
    );
    return PyFloat_FromDouble(count);
}

static PyObject *PyBallTree_dualcount_range(PyObject *self, PyObject *args) {
    PyBallTree *pytree = (PyBallTree *)self;

    // Parse the Python arguments
    PyBallTree *other_tree;
    double rmin, rmax;
    if (!PyArg_ParseTuple(args, "O!dd", &PyBallTreeType, &other_tree, &rmin, &rmax)) {
        return NULL;
    }

    double count = balltree_dualcount_range(
        pytree->balltree,
        other_tree->balltree,
        rmin,
        rmax
    );
    return PyFloat_FromDouble(count);
}

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

    if (PyType_Ready(&PyBallTreeType) < 0) {
        return NULL;
    }
    Py_INCREF(&PyBallTreeType);

    PyObject *module = PyModule_Create(&pyballtree);
    if (module == NULL) {
        return NULL;
    }
    PyModule_AddObject(module, "BallTree", (PyObject *)&PyBallTreeType);

    return module;
}
