#include <Python.h>
#include <numpy/arrayobject.h>

#include <stdio.h>

#include "point.h"
#include "balltree.h"

struct PyBallTree;

static Point *PyIter_to_point(PyObject *point_iter, double weight);
PyObject *ptbuf_get_numpy_view(PointBuffer *buffer);

// PyBallTree: constructors & deallocators
static PyObject *PyBallTree_from_random(PyTypeObject *type, PyObject *args, PyObject *kwds);
static PyObject *pyballtree_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
static void pyballtree_dealloc(PyObject *self);
// PyBallTree: property implementation
static PyObject *pyballtree_get_data(PyObject *self, void *closure);
static PyObject *pyballtree_get_leafsize(PyObject *self, void *closure);
static PyObject *pyballtree_get_sum_weight(PyObject *self, void *closure);
static PyObject *pyballtree_get_center(PyObject *self, void *closure);
static PyObject *pyballtree_get_radius(PyObject *self, void *closure);
// PyBallTree: method implementations
static PyObject *PyBallTree_count_nodes(PyObject *self, PyObject *args);
static PyObject *PyBallTree_count_radius(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *PyBallTree_count_range(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *PyBallTree_dualcount_radius(PyObject *self, PyObject *args);
static PyObject *PyBallTree_dualcount_range(PyObject *self, PyObject *args);

// helper functions ////////////////////////////////////////////////////////////

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

PyObject *ptbuf_get_numpy_view(PointBuffer *buffer) {
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

// PyBallTree definition ///////////////////////////////////////////////////////

typedef struct {
    PyObject_HEAD
    BallTree *balltree;
} PyBallTree;

static PyMethodDef pyballtree_methods[] = {
    {
        "from_random",
        (PyCFunction)PyBallTree_from_random,
        METH_CLASS | METH_VARARGS | METH_KEYWORDS,
        "Build a PyBallTree instance"
    },
    {
        "count_nodes",
        (PyCFunction)PyBallTree_count_nodes,
        METH_NOARGS,
        "Count nodes contained in tree"
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
        METH_VARARGS | METH_KEYWORDS,
        "Count pairs within a radius with the dualtree algorithm"
    },
    {
        "dualcount_range",
        (PyCFunction)PyBallTree_dualcount_range,
        METH_VARARGS | METH_KEYWORDS,
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
    .tp_name = "pyballtree.PyBallTree",
    .tp_doc = "Python wrapper for C BallTree",
    .tp_basicsize = sizeof(PyBallTree),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = pyballtree_new,
    .tp_dealloc = pyballtree_dealloc,
    .tp_methods = pyballtree_methods,
    .tp_getset = pyballtree_getset,
};

// PyBallTree: constructors & deallocators /////////////////////////////////////

static PyObject *PyBallTree_from_random(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    // parse the python arguments
    double low, high;
    int size;
    int leafsize = -1;

    static char *kwlist[] = {"low", "high", "size", "leafsize", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ddi|i", kwlist, &low, &high, &size, &leafsize)) {
        return NULL;
    }

    // generate random data points
    PointBuffer *buffer = ptbuf_gen_random(low, high, size);
    if (buffer == NULL) {
        return NULL;
    }

    // build the balltree
    BallTree *tree;
    if (leafsize == -1) {
        tree = balltree_build(buffer);
    } else {
        tree = balltree_build_leafsize(buffer, leafsize);
    }
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

static PyObject *pyballtree_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    PyBallTree *self;
    self = (PyBallTree *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->balltree = NULL;  // TODO
        if (self->balltree == NULL) {
            Py_DECREF(self);
            return NULL;
        }
    }
    return (PyObject *)self;
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

static PyObject *PyBallTree_count_nodes(PyObject *self, PyObject *args) {
    PyBallTree *pytree = (PyBallTree *)self;
    int count = balltree_count_nodes(pytree->balltree);
    return PyLong_FromLong(count);
}

static PyObject *PyBallTree_count_radius(PyObject *self, PyObject *args, PyObject *kwargs) {
    PyBallTree *pytree = (PyBallTree *)self;

    // Parse the Python arguments
    PyObject *coord_iter;
    double radius;
    double weight = 1.0;  // weight is optional
    static char *kwlist[] = {"point", "radius", "weight", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Od|d", kwlist, &coord_iter, &radius, &weight)) {
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

static PyObject *PyBallTree_count_range(PyObject *self, PyObject *args, PyObject *kwargs) {
    PyBallTree *pytree = (PyBallTree *)self;

    // Parse the Python arguments
    PyObject *coord_iter;
    double rmin, rmax;
    double weight = 1.0;  // weight is optional
    static char *kwlist[] = {"point", "rmin", "rmax", "weight", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Odd|d", kwlist, &coord_iter, &rmin, &rmax, &weight)) {
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

    double count = balltree_dualcount_radius(pytree->balltree, other_tree->balltree, radius);
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

    double count = balltree_dualcount_range(pytree->balltree, other_tree->balltree, rmin, rmax);
    return PyFloat_FromDouble(count);
}

// module configuration ////////////////////////////////////////////////////////

static struct PyModuleDef pyballtree = {
    PyModuleDef_HEAD_INIT,
    .m_name = "pyballtree",
    .m_doc = "Fast balltree implementation",
    .m_size = -1,
};

PyMODINIT_FUNC PyInit_pyballtree(void) {
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
    PyModule_AddObject(module, "PyBallTree", (PyObject *)&PyBallTreeType);

    return module;
}
