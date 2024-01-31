#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>

#include <stdint.h>
#include <stdio.h>

#include "point.h"
#include "balltree.h"
#include "balltree_macros.h"

PyDoc_STRVAR(
    // .. py:class::
    balltree_doc,
    "BallTree(xyz: ArrayLike, weight: ArrayLike | None = None, leafsize: int = default_leafsize)\n\n"
    "A wrapper for the `C` implementation of BallTree.\n\n"
    "The tree is implemented for 3-dim coordinates and an Euclidean distance\n"
    "metric. It provides different alorigthms to count neighbours within a\n"
    "single or multiple query points.\n\n"
    "The data point(s) ``xyz`` can be a numpy array of shape (3,) or (N, 3),\n"
    "or an equivalent python object. The optional ``weights`` can be a float\n"
    "or a 1-dim sequence of matching length, the optional ``leafsize``\n"
    "determines when the tree query algorithms switch from traversal to brute\n"
    "force."
);

typedef double (*count_radius_func)(const BallTree *, const Point *, double);
typedef void (*count_range_func)(const BallTree *, const Point *, DistHistogram *);

typedef struct {
    NpyIter *iter;
    NpyIter_IterNextFunc *next;
    char **dataptr;
    npy_intp *stride;
    npy_intp *size;
    npy_intp idx;
} NpyIterHelper;

typedef struct {
    int64_t size;
    PyArrayObject *xyz_arr;
    NpyIterHelper *xyz_iter;
    PyArrayObject *weight_arr;
    double *weight_buffer;
} InputIterData;

typedef struct {
    PyObject_HEAD
    BallTree *balltree;
} PyBallTree;
static PyTypeObject PyBallTreeType;

// helper functions
static NpyIterHelper *npyiterhelper_new(PyArrayObject *xyz_arr);
static void npyiterhelper_free(NpyIterHelper *iter);
static inline int iter_get_next_xyz(NpyIterHelper *iter, double *x, double *y, double *z);
static const char *PyString_to_char(PyObject* py_string);
static PyArrayObject *ensure_numpy_array_double(PyObject *obj);
PyArrayObject *numpy_array_add_dim(PyArrayObject* array);
static PyArrayObject *ensure_numpy_array_1dim_double(PyObject *weight_obj);
static PyArrayObject *xyz_ensure_2dim_double(PyObject *xyz_obj);
static PyArrayObject *weight_matched_1dim_double(PyObject *weight_obj, npy_intp length);
static InputIterData *inputiterdata_new(PyObject *xyz_obj, PyObject *weight_obj);
void inputiterdata_free(InputIterData *);
static PointBuffer *ptbuf_from_PyObjects(PyObject *xyz_obj, PyObject *weight_obj);
static PyObject *ptbuf_get_numpy_view(PointBuffer *buffer);
static PyObject *statvec_get_numpy_array(StatsVector *vec);
static DistHistogram *disthistogram_from_PyObject(PyObject *edges_obj);
static PyObject *PyObject_from_disthistogram(DistHistogram *hist);
static PyObject *PyBallTree_accumulate_radius(PyBallTree *self, count_radius_func accumulator, PyObject *xyz_obj, double radius, PyObject *weight_obj);
static PyObject *PyBallTree_accumulate_range(PyBallTree *self, count_range_func accumulator, PyObject *xyz_obj, PyObject *edges_obj, PyObject *weight_obj);

// PyBallTree: constructors & deallocators
static PyObject *PyBallTree_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
static void PyBallTree_dealloc(PyObject *self);
static int PyBallTree_init(PyBallTree *self, PyObject *args, PyObject *kwds);
static PyObject *PyBallTree_from_random(PyTypeObject *cls, PyObject *args, PyObject *kwds);
static PyObject *PyBallTree_from_file(PyTypeObject *cls, PyObject *args);
// PyBallTree: property implementation
static PyObject *PyBallTree_get_data(PyBallTree *self, void *closure);
static PyObject *PyBallTree_get_num_points(PyBallTree *self, void *closure);
static PyObject *PyBallTree_get_leafsize(PyBallTree *self, void *closure);
static PyObject *PyBallTree_get_sum_weight(PyBallTree *self, void *closure);
static PyObject *PyBallTree_get_center(PyBallTree *self, void *closure);
static PyObject *PyBallTree_get_radius(PyBallTree *self, void *closure);
// PyBallTree: method implementations
static PyObject *PyBallTree_str(PyObject *self);
static PyObject *PyBallTree_to_file(PyBallTree *self, PyObject *args);
static PyObject *PyBallTree_count_nodes(PyBallTree *self);
static PyObject *PyBallTree_get_node_data(PyBallTree *self);
static PyObject *PyBallTree_brute_radius(PyBallTree *self, PyObject *args, PyObject *kwargs);
static PyObject *PyBallTree_count_radius(PyBallTree *self, PyObject *args, PyObject *kwargs);
static PyObject *PyBallTree_dualcount_radius(PyBallTree *self, PyObject *args);
static PyObject *PyBallTree_brute_range(PyBallTree *self, PyObject *args, PyObject *kwargs);
static PyObject *PyBallTree_count_range(PyBallTree *self, PyObject *args, PyObject *kwargs);
static PyObject *PyBallTree_dualcount_range(PyBallTree *self, PyObject *args);


// helper functions ////////////////////////////////////////////////////////////

static NpyIterHelper *npyiterhelper_new(PyArrayObject *xyz_arr) {
    NpyIterHelper *iterhelper = malloc(sizeof(NpyIterHelper));
    if (iterhelper == NULL) {
        PyErr_SetString(PyExc_MemoryError, "failed to allocate NpyIterHelper");
        return NULL;
    }

    NpyIter *iter = NpyIter_New(
        xyz_arr,
        NPY_ITER_READONLY | NPY_ITER_EXTERNAL_LOOP,
        NPY_KEEPORDER,
        NPY_NO_CASTING,
        NULL
    );
    if (iter == NULL) {
        free(iterhelper);
        return NULL;
    }
    iterhelper->iter = iter;

    iterhelper->next = NpyIter_GetIterNext(iter, NULL);
    iterhelper->dataptr = NpyIter_GetDataPtrArray(iter);
    iterhelper->stride = NpyIter_GetInnerStrideArray(iter);
    iterhelper->size = NpyIter_GetInnerLoopSizePtr(iter);
    iterhelper->idx = 0;
    return iterhelper;
}

static void npyiterhelper_free(NpyIterHelper *iter) {
    if (iter->iter != NULL) {
        NpyIter_Deallocate(iter->iter);
    }
    free(iter);
}

static inline int iter_get_next_xyz(NpyIterHelper *iter, double *x, double *y, double *z) {
    // NOTE: If a slice with negative increment is used, this loop will still
    // iterate in the order the data is stored in memory. Nevertheless the
    // correct elements are selected, just in "reversed" order.
    if (iter->idx >= *iter->size) {  // loop while idx < size
        if (!iter->next(iter->iter)) {  // update inner loop pointers
            return 0;
        }
        iter->idx = 0;
    }

    double *xyz = (double *)*iter->dataptr;
    // inner loop count of (x, y, z) array is guaranteed to be multiple of 3
    *x = xyz[(iter->idx)++];
    *y = xyz[(iter->idx)++];
    *z = xyz[(iter->idx)++];
    return 1;
}

static const char *PyString_to_char(PyObject* py_string) {
    if (!PyUnicode_Check(py_string)) {
        PyErr_SetString(PyExc_TypeError, "input must be of type 'str'");
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

static PyArrayObject *ensure_numpy_array_double(PyObject *obj) {
    PyArrayObject *array, *array_typed;
    const npy_int type_expected = NPY_DOUBLE;
    const char type_err_msg[] = "cannot convert 'xyz' to array with type double";

    // ensure that array type, either by creating a new object or by incrementing
    // the ref. coun of on input object
    if (!PyArray_Check(obj)) {
        // convert input object to numpy array of type double
        array = (PyArrayObject *)PyArray_FROM_OTF(obj, type_expected, NPY_ARRAY_CARRAY_RO);
        if (array == NULL) {
            PyErr_SetString(PyExc_TypeError, type_err_msg);
            return NULL;
        }
    } else {
        array = (PyArrayObject *)obj;
        Py_INCREF(array);  // consistency with branch above
    }

    // ensure array is of type double
    if (PyArray_TYPE(array) != type_expected) {
        array_typed = (PyArrayObject *)PyArray_Cast(array, type_expected);
        if (array == NULL) {
            PyErr_SetString(PyExc_TypeError, type_err_msg);
            return NULL;
        }
    } else {
        array_typed = array;
        Py_INCREF(array_typed);
    }
    Py_DECREF(array);

    return array_typed;
}

PyArrayObject *numpy_array_add_dim(PyArrayObject* array) {
    if (!PyArray_Check(array)) {
        PyErr_SetString(PyExc_TypeError, "input is not a numpy array");
        return NULL;
    }
    npy_intp ndim = PyArray_NDIM(array);
    npy_intp *shape = PyArray_SHAPE(array);
    PyArrayObject *reshaped;

    // create new shape with additional dimension of length 1
    npy_intp ndim_new = ndim + 1;
    npy_intp *shape_new = malloc(ndim_new * sizeof(npy_intp));
    if (shape_new == NULL) {
        PyErr_SetString(PyExc_MemoryError, "failed to allocate new shape");
        return NULL;
    }
    shape_new[0] = 1;
    for (npy_intp i = 0; i < ndim; ++i) {
        shape_new[i + 1] = shape[i];
    }

    // reshape the array or create a new one
    PyArray_Dims dims_new = {shape_new, ndim_new};
    reshaped = (PyArrayObject *)PyArray_Newshape(array, &dims_new, NPY_CORDER);
    free(shape_new);
    if (reshaped == NULL) {
        PyErr_SetString(PyExc_MemoryError, "failed to reshape array");
    }
    return reshaped;
}

static PyArrayObject *ensure_numpy_array_1dim_double(PyObject *weight_obj) {
    PyObject *weight_seq;
    PyArrayObject *weight_arr;
    const npy_int ndim_expected = 1;

    // ensure that input is a sequence at minimum
    if (PyArray_IsAnyScalar(weight_obj)) {
        weight_seq = Py_BuildValue("(O)", weight_obj);
        if (weight_seq == NULL) {
            return NULL;
        }
    } else {
        weight_seq = weight_obj;
        Py_INCREF(weight_seq);  // consistency with branch above
    }

    weight_arr = ensure_numpy_array_double(weight_seq);
    Py_DECREF(weight_seq);
    if (weight_arr == NULL) {
        return NULL;
    }

    // check that array has correct number of dimensions
    npy_int ndim = PyArray_NDIM(weight_arr);
    if (ndim != ndim_expected) {
        PyErr_SetString(PyExc_ValueError, "'weight' must be scalar or of shape (N,)");
        Py_DECREF(weight_arr);
        return NULL;
    }
    return weight_arr;
}

static PyArrayObject *xyz_ensure_2dim_double(PyObject *xyz_obj) {
    PyArrayObject *xyz_arr, *xyz_arr_2dim;
    const npy_int dim1_expected = 3;
    const char shape_err_msg[] = "'xyz' must be of shape (3,) or (N, 3)";

    xyz_arr = ensure_numpy_array_double(xyz_obj);
    if (xyz_arr == NULL) {
        return NULL;
    }

    // ensure that array has correct number of dimensions
    npy_int ndim = PyArray_NDIM(xyz_arr);
    if (ndim == 1) {
        xyz_arr_2dim = numpy_array_add_dim(xyz_arr);
        Py_DECREF(xyz_arr);
        if (xyz_arr_2dim == NULL) {
            return NULL;
        }
    } else if (ndim == 2) {
        xyz_arr_2dim = xyz_arr;
        // implicitly: Py_INCREF(xyz_arr_2dim); Py_DECREF(xyz_arr);
    } else {
        PyErr_SetString(PyExc_ValueError, shape_err_msg);
        Py_DECREF(xyz_arr);
        return NULL;
    }

    // check that array has correct shape
    if (PyArray_DIM(xyz_arr_2dim, 1) != dim1_expected) {
        PyErr_SetString(PyExc_ValueError, shape_err_msg);
        Py_DECREF(xyz_arr_2dim);
        return NULL;
    }
    return xyz_arr_2dim;
}

static PyArrayObject *weight_matched_1dim_double(PyObject *weight_obj, npy_intp length) {
    PyArrayObject *weight_arr;

    // create array of correct length initialised to 1.0
    if (weight_obj == Py_None) {
        npy_intp shape[1] = {length};
        weight_arr = (PyArrayObject *)PyArray_EMPTY(1, shape, NPY_DOUBLE, 0);
        if (weight_arr == NULL) {
            PyErr_SetString(PyExc_MemoryError, "failed to allocate default weight array");
            return NULL;
        }

        double *weight_buffer = PyArray_DATA(weight_arr);
        for (npy_intp i = 0; i < length; ++i) {
            weight_buffer[i] = 1.0;
        }
    }

    // try converting input to numpy array with correct type and shape
    else {
        weight_arr = ensure_numpy_array_1dim_double(weight_obj);
        if (weight_arr == NULL) {
            return NULL;
        }
        if (PyArray_DIM(weight_arr, 0) != length) {
            PyErr_SetString(PyExc_ValueError, "'xyz' and 'weight' must have the same length");
            Py_DECREF(weight_arr);
            return NULL;
        }
    }

    return weight_arr;
}

static InputIterData *inputiterdata_new(PyObject *xyz_obj, PyObject *weight_obj) {
    InputIterData *data = calloc(1, sizeof(InputIterData));
    if (data == NULL) {
        return NULL;
    }

    // ensure that inputs are numpy arrays of correct type and shape
    data->xyz_arr = xyz_ensure_2dim_double(xyz_obj);
    if (data->xyz_arr == NULL) {
        inputiterdata_free(data);
        return NULL;
    }
    data->size = PyArray_DIM(data->xyz_arr, 0);
    if (data->size == 0) {
        PyErr_SetString(PyExc_ValueError, "'xyz' needs to contain at least one element");
        inputiterdata_free(data);
        return NULL;
    }
    data->weight_arr = weight_matched_1dim_double(weight_obj, data->size);
    if (data->weight_arr == NULL) {
        inputiterdata_free(data);
        return NULL;
    }

    // get an iterator for input data
    data->xyz_iter = npyiterhelper_new(data->xyz_arr);
    if (data->xyz_iter == NULL) {
        inputiterdata_free(data);
        return NULL;
    }
    data->weight_buffer = PyArray_DATA(data->weight_arr);

    return data;
}

void inputiterdata_free(InputIterData *data) {
    Py_XDECREF(data->xyz_arr);
    if (data->xyz_iter != NULL) {
        npyiterhelper_free(data->xyz_iter);
    }
    Py_XDECREF(data->weight_arr);
    free(data);
}

static PointBuffer *ptbuf_from_PyObjects(PyObject *xyz_obj, PyObject *weight_obj) {
    InputIterData *data = inputiterdata_new(xyz_obj, weight_obj);
    if (data == NULL) {
        return NULL;
    }
    // build and fill point buffer
    PointBuffer *buffer = ptbuf_new(data->size);
    if (buffer == NULL) {
        inputiterdata_free(data);
        return NULL;
    }
    int64_t idx = 0;
    double x, y, z;
    while (iter_get_next_xyz(data->xyz_iter, &x, &y, &z)) {
        buffer->points[idx] = (Point){x, y, z, data->weight_buffer[idx]};
        ++idx;
    }
    inputiterdata_free(data);
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
    int result = PyArray_DescrConverter(arr_dtype, &arr_descr);
    Py_DECREF(arr_dtype);
    if (result != NPY_SUCCEED) {
        return NULL;
    }

    // get an array view
    PyObject *view = PyArray_NewFromDescr(
        &PyArray_Type,
        arr_descr,  // reference stolen, no Py_DECREF needed
        ndim,
        shape,
        NULL,
        buffer->points,
        NPY_ARRAY_CARRAY_RO,
        NULL
    );
    if (view == NULL) {
        Py_DECREF(arr_descr);  // not sure if that is necessary
    }
    return view;
}

static DistHistogram *disthistogram_from_PyObject(PyObject *edges_obj) {
    PyArrayObject *edges_arr = ensure_numpy_array_1dim_double(edges_obj);
    if (edges_arr == NULL) {
        return NULL;
    }
    int64_t num_edges = (int64_t)PyArray_DIM(edges_arr, 0);
    double *edges_buffer = PyArray_DATA(edges_arr);
    DistHistogram *hist = hist_new(num_edges, edges_buffer);
    Py_DECREF(edges_arr);
    return hist;
}

static PyObject *PyObject_from_disthistogram(DistHistogram *hist) {
    npy_intp dims = {hist->size};
    PyObject *pycount = PyArray_SimpleNew(1, &dims, NPY_DOUBLE);
    if (pycount == NULL) {
        return NULL;
    }
    double *count_buffer = PyArray_DATA(pycount);
    for (int64_t i = 0; i < hist->size; ++i) {
        count_buffer[i] = hist->sum_weight[i];
    }
    return pycount;
}

static PyObject *statvec_get_numpy_array(StatsVector *vec) {
    const npy_intp ndim = 1;
    npy_intp shape[1] = {vec->size};

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
    memcpy(ptr, vec->stats, sizeof(NodeStats) * vec->size);
    return array;
}

static PyObject *PyBallTree_accumulate_radius(
    PyBallTree *self,
    count_radius_func accumulator,
    PyObject *xyz_obj,
    double radius,
    PyObject *weight_obj
) {
    InputIterData *data = inputiterdata_new(xyz_obj, weight_obj);
    if (data == NULL) {
        return NULL;
    }
    // count neighbours for all inputs
    double count = 0.0;
    int64_t idx = 0;
    Point point;
    while (iter_get_next_xyz(data->xyz_iter, &point.x, &point.y, &point.z)) {
        point.weight = data->weight_buffer[idx];
        count += accumulator(self->balltree, &point, radius);
        ++idx;
    }
    inputiterdata_free(data);
    return PyFloat_FromDouble(count);
}

static PyObject *PyBallTree_accumulate_range(
    PyBallTree *self,
    count_range_func accumulator,
    PyObject *xyz_obj,
    PyObject *edges_obj,
    PyObject *weight_obj
) {
    InputIterData *data = inputiterdata_new(xyz_obj, weight_obj);
    if (data == NULL) {
        return NULL;
    }
    DistHistogram *hist = disthistogram_from_PyObject(edges_obj);
    if (hist == NULL) {
        inputiterdata_free(data);
        return NULL;
    }
    // count neighbours for all inputs
    int64_t idx = 0;
    Point point;
    while (iter_get_next_xyz(data->xyz_iter, &point.x, &point.y, &point.z)) {
        point.weight = data->weight_buffer[idx];
        accumulator(self->balltree, &point, hist);
        ++idx;
    }
    PyObject *pycounts = PyObject_from_disthistogram(hist);
    inputiterdata_free(data);
    hist_free(hist);
    return pycounts;
}

// PyBallTree: constructors & deallocators /////////////////////////////////////

static PyObject *PyBallTree_new(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
    PyBallTree *self = (PyBallTree *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->balltree = NULL;
    }
    return (PyObject *)self;
}

static void PyBallTree_dealloc(PyObject *self) {
    PyBallTree *pytree = (PyBallTree *)self;
    if (pytree->balltree != NULL) {
        balltree_free(pytree->balltree);
    }
    Py_TYPE(self)->tp_free(self);
}

static int PyBallTree_init(
    PyBallTree *self,
    PyObject *args,
    PyObject *kwargs
) {
    static char *kwlist[] = {"xyz", "weight", "leafsize", NULL};
    PyObject *xyz_obj, *weight_obj = Py_None;
    int leafsize = DEFAULT_LEAFSIZE;
    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, "O|Oi", kwlist,
            &xyz_obj, &weight_obj, &leafsize)
    ) {
        return -1;
    }
    PointBuffer *buffer = ptbuf_from_PyObjects(xyz_obj, weight_obj);
    if (buffer == NULL) {
        return -1;
    }

    // build the balltree (takes ownership of buffer)
    BallTree *tree = balltree_build_nocopy(buffer, leafsize);
    if (tree == NULL) {
        ptbuf_free(buffer);
        return -1;
    }
    tree->data_owned = 1;  // ownership transfer, buffer deallocated with tree
    self->balltree = tree;
    return 0;
}

PyDoc_STRVAR(
    // .. py:method::
    //     :classmethod:
    from_random_doc,
    "from_random(cls, low: float, high: float, size: int, leafsize: int = default_leafsize) -> BallTree\n"
    "--\n\n"
    "Build a new BallTree instance from randomly generated points.\n\n"
    "The (x, y, z) coordinates are generated uniformly in the interval\n"
    "[`low`, `high`), `size` controlls the number of points generated. The\n"
    "optional `leafsize` determines when the tree query algorithms switch from\n"
    "traversal to brute force."
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

    // build the balltree (takes ownership of buffer)
    BallTree *tree = balltree_build_nocopy(buffer, leafsize);
    if (tree == NULL) {
        ptbuf_free(buffer);
        return NULL;
    }
    tree->data_owned = 1;  // ownership transfer, buffer deallocated with tree

    PyBallTree *self = (PyBallTree *)type->tp_alloc(type, 0);
    if (self == NULL) {
        balltree_free(tree);
    } else {
        self->balltree = tree;
    }
    return (PyObject *)self;
}

PyDoc_STRVAR(
    // .. py:method::
    //     :classmethod:
    from_file_doc,
    "from_file(cls, path: str) -> BallTree\n"
    "--\n\n"
    "Load back a class instance from a file created by the ``to_file()`` method."
);

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
    if (self == NULL) {
        balltree_free(tree);
    } else {
        self->balltree = tree;
    }
    return (PyObject *)self;
}

// PyBallTree: property implementations ////////////////////////////////////////

PyDoc_STRVAR(
    // .. py:property:: data
    //    :type: NDArray
    get_data_doc,
    "Get a view of the data points stored in the tree.\n\n"
    "Note that the points are reordered when the tree is constructed."
);

static PyObject *PyBallTree_get_data(PyBallTree *self, void *closure) {
    PyObject *view = ptbuf_get_numpy_view(self->balltree->data);
    if (view == NULL) {
        return NULL;
    }
    // make sure that BallTree instance with the underlying buffer stays alive
    // as long as view is in use
    Py_INCREF(self);  // next line steals one
    if (PyArray_SetBaseObject((PyArrayObject *)view, (PyObject *)self)) {
        Py_DECREF(self);
        Py_DECREF(view);
        return NULL;
    }
    return view;
}

PyDoc_STRVAR(
    // .. py:attribute:: num_points
    //    :type: int
    get_num_points_doc,
    "Get the number of data points stored in the tree."
);

static PyObject *PyBallTree_get_num_points(PyBallTree *self, void *closure) {
    return PyLong_FromLong(self->balltree->data->size);
}

PyDoc_STRVAR(
    // .. py:attribute:: leafsize
    //    :type: int
    get_leafsize_doc,
    "Get the leaf size of the tree."
);

static PyObject *PyBallTree_get_leafsize(PyBallTree *self, void *closure) {
    return PyLong_FromLong(self->balltree->leafsize);
}

PyDoc_STRVAR(
    // .. py:attribute:: sum_weight
    //    :type: float
    get_sum_weight_doc,
    "Get the sum over the data point weights stored in the tree."
);

static PyObject *PyBallTree_get_sum_weight(PyBallTree *self, void *closure) {
    return PyFloat_FromDouble(self->balltree->root->sum_weight);
}

PyDoc_STRVAR(
    // .. py:attribute:: center
    //    :type: tuple(float, float, float)
    get_center_doc,
    "Get the coordinates of the root node of the tree."
);

static PyObject *PyBallTree_get_center(PyBallTree *self, void *closure) {
    PyObject *tuple = NULL;
    PyObject *x = PyFloat_FromDouble(self->balltree->root->ball.x);
    PyObject *y = PyFloat_FromDouble(self->balltree->root->ball.y);
    PyObject *z = PyFloat_FromDouble(self->balltree->root->ball.z);

    if (x != NULL && y != NULL && z != NULL) {
        tuple = PyTuple_Pack(3, x, y, z);  // adds additional references to x, y, z
    }
    Py_XDECREF(x);
    Py_XDECREF(y);
    Py_XDECREF(z);
    return tuple;
}

PyDoc_STRVAR(
    // .. py:attribute:: radius
    //    :type: float
    get_radius_doc,
    "Get the radius of the root node of the tree."
);

static PyObject *PyBallTree_get_radius(PyBallTree *self, void *closure) {
    return PyFloat_FromDouble(self->balltree->root->ball.radius);
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
        "BallTree(num_points=%ld, radius=%lf, center=(%lf, %lf, %lf))",
        tree->data->size,
        node->ball.radius,
        node->ball.x,
        node->ball.y,
        node->ball.z
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
    // .. py:method::
    to_file_doc,
    "to_file(self, path: str) -> None\n"
    "--\n\n"
    "Store a representation of the tree instance in a binary file."
);

static PyObject *PyBallTree_to_file(PyBallTree *self, PyObject *args) {
    PyObject *py_string;
    if (!PyArg_ParseTuple(args, "O!", &PyUnicode_Type, &py_string)) {
        return NULL;
    }
    const char *path = PyString_to_char(py_string);  // buffer is managed by py_string
    if (path == NULL) {
        return NULL;
    }

    if (balltree_to_file(self->balltree, path) == BTR_FAILED) {
        return NULL;
    }
    Py_RETURN_NONE;
}

PyDoc_STRVAR(
    // .. py:method::
    count_nodes_doc,
    "count_nodes(self) -> int\n"
    "--\n\n"
    "Get a count of all nodes of the tree, including the root node."
);

static PyObject *PyBallTree_count_nodes(PyBallTree *self) {
    int count = balltree_count_nodes(self->balltree);
    return PyLong_FromLong(count);
}

PyDoc_STRVAR(
    // .. py:method::
    get_node_data_doc,
    "get_node_data(self) -> NDArray\n"
    "--\n\n"
    "Collect the meta data of all tree nodes in a numpy array.\n\n"
    "The array fields record `depth` (starting from the root node),\n"
    "`num_points`, `sum_weight`, `x`, `y`, `z` (node center) and node `radius`."
);

static PyObject *PyBallTree_get_node_data(PyBallTree *self) {
    StatsVector *vec = balltree_collect_stats(self->balltree);
    if (vec == NULL) {
        return NULL;
    }

    PyObject *array = statvec_get_numpy_array(vec);
    statvec_free(vec);
    return array;
}

PyDoc_STRVAR(
    // .. py:method::
    brute_radius_doc,
    "brute_radius(self, xyz: ArrayLike, radius: float, weight: ArrayLike | None = None) -> float\n"
    "--\n\n"
    "Count neighbours within a given radius using brute force.\n\n"
    "The query point(s) ``xyz`` can be a numpy array of shape (3,) or (N, 3),\n"
    "or an equivalent python object. The optional ``weights`` can be a float\n"
    "or a 1-dim sequence of matching length."
);

static PyObject *PyBallTree_brute_radius(
    PyBallTree *self,
    PyObject *args,
    PyObject *kwargs
) {
    static char *kwlist[] = {"xyz", "radius", "weight", NULL};
    PyObject *xyz_obj, *weight_obj = Py_None;
    double radius;
    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, "Od|O", kwlist,
            &xyz_obj, &radius, &weight_obj)
    ) {
        return NULL;
    }
    return PyBallTree_accumulate_radius(
        self,
        balltree_brute_radius,
        xyz_obj,
        radius,
        weight_obj
    );
}

PyDoc_STRVAR(
    // .. py:method::
    count_radius_doc,
    "count_radius(self, xyz: ArrayLike, radius: float, weight: ArrayLike | None = None) -> float\n"
    "--\n\n"
    "Count neighbours within a given radius using tree traversal.\n\n"
    "The query point(s) ``xyz`` can be a numpy array of shape (3,) or (N, 3),\n"
    "or an equivalent python object. The optional ``weights`` can be a float\n"
    "or a 1-dim sequence of matching length."
);

static PyObject *PyBallTree_count_radius(
    PyBallTree *self,
    PyObject *args,
    PyObject *kwargs
) {
    static char *kwlist[] = {"xyz", "radius", "weight", NULL};
    PyObject *xyz_obj, *weight_obj = Py_None;
    double radius;
    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, "Od|O", kwlist,
            &xyz_obj, &radius, &weight_obj)
    ) {
        return NULL;
    }
    return PyBallTree_accumulate_radius(
        self,
        balltree_count_radius,
        xyz_obj,
        radius,
        weight_obj
    );
}

PyDoc_STRVAR(
    // .. py:method::
    dualcount_radius_doc,
    "dualcount_radius(self, tree: BallTree, radius: float) -> NDArray\n"
    "--\n\n"
    "Count all pairs within a given radius from points in another tree.\n\n"
    "The pairs between the two trees are computed with the efficient dualtree\n"
    "algorithm."
);

static PyObject *PyBallTree_dualcount_radius(PyBallTree *self, PyObject *args) {
    PyBallTree *other_tree;
    double radius;
    if (!PyArg_ParseTuple(args, "O!d", &PyBallTreeType, &other_tree, &radius)) {
        return NULL;
    }
    double count = balltree_dualcount_radius(
        self->balltree,
        other_tree->balltree,
        radius
    );
    return PyFloat_FromDouble(count);
}

PyDoc_STRVAR(
    // .. py:method::
    brute_range_doc,
    "brute_range(self, xyz: ArrayLike, radii: ArrayLike, weight: ArrayLike | None = None) -> NDArray\n"
    "--\n\n"
    "Count neighbours within a sequence of radii using brute force.\n\n"
    "The query point(s) ``xyz`` can be a numpy array of shape (3,) or (N, 3),\n"
    "or an equivalent python object. The ``radii`` must either be a float or\n"
    "monotonic sequence. The optional ``weights`` can be a float or a 1-dim\n"
    "sequence of matching length.\n\n"
    "Returns an array of counts. The first element contains the count of all\n"
    "neighbours ``0 <= r <= r_1``, the following values contain the\n"
    "incremental counts ``r_i-1 <= r <= r_i``."
);

static PyObject *PyBallTree_brute_range(
    PyBallTree *self,
    PyObject *args,
    PyObject *kwargs
) {
    static char *kwlist[] = {"xyz", "radii", "weight", NULL};
    PyObject *xyz_obj, *edges_obj, *weight_obj = Py_None;
    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, "OO|O", kwlist,
            &xyz_obj, &edges_obj, &weight_obj)
    ) {
        return NULL;
    }
    return PyBallTree_accumulate_range(
        self,
        balltree_brute_range,
        xyz_obj,
        edges_obj,
        weight_obj
    );
}

PyDoc_STRVAR(
    // .. py:method::
    count_range_doc,
    "count_range(self, xyz: ArrayLike, radii: ArrayLike, weight: ArrayLike | None = None) -> NDArray\n"
    "--\n\n"
    "Count neighbours within a sequence of radii using tree traversal.\n\n"
    "The query point(s) ``xyz`` can be a numpy array of shape (3,) or (N, 3),\n"
    "or an equivalent python object. The ``radii`` must either be a float or\n"
    "monotonic sequence. The optional ``weights`` can be a float or a 1-dim\n"
    "sequence of matching length.\n\n"
    "Returns an array of counts. The first element contains the count of all\n"
    "neighbours ``0 <= r <= r_1``, the following values contain the\n"
    "incremental counts ``r_i-1 <= r <= r_i``."
);

static PyObject *PyBallTree_count_range(
    PyBallTree *self,
    PyObject *args,
    PyObject *kwargs
) {
    static char *kwlist[] = {"xyz", "radii", "weight", NULL};
    PyObject *xyz_obj, *edges_obj, *weight_obj = Py_None;
    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, "OO|O", kwlist,
            &xyz_obj, &edges_obj, &weight_obj)
    ) {
        return NULL;
    }
    return PyBallTree_accumulate_range(
        self,
        balltree_count_range,
        xyz_obj,
        edges_obj,
        weight_obj
    );
}

PyDoc_STRVAR(
    // .. py:method::
    dualcount_range_doc,
    "dualcount_range(self, other: BallTree, radii: ArrayLike) -> NDArray\n"
    "--\n\n"
    "Count all pairs within a sequence of radii from points in another tree.\n\n"
    "The pairs between the two trees are computed with the efficient dualtree\n"
    "algorithm. The ``radii`` must either be a float or monotonic sequence.\n"
    "The optional ``weights`` can be a float or a 1-dim.\n\n"
    "Returns an array of counts. The first element contains the count of all\n"
    "neighbours ``0 <= r <= r_1``, the following values contain the\n"
    "incremental counts ``r_i-1 <= r <= r_i``."
);

static PyObject *PyBallTree_dualcount_range(PyBallTree *self, PyObject *args) {
    PyBallTree *other_tree;
    PyObject *edges_obj;
    if (!PyArg_ParseTuple(args, "O!O", &PyBallTreeType, &other_tree, &edges_obj)) {
        return NULL;
    }

    DistHistogram *hist = disthistogram_from_PyObject(edges_obj);
    if (hist == NULL) {
        return NULL;
    }
    balltree_dualcount_range(self->balltree, other_tree->balltree, hist);
    PyObject *pycounts = PyObject_from_disthistogram(hist);
    hist_free(hist);
    return pycounts;
}

// PyBallTree definition ///////////////////////////////////////////////////////

static PyGetSetDef PyBallTree_getset[] = {
    {
        .name = "data",
        .get = (getter)PyBallTree_get_data,
        .doc = get_data_doc,
    },
    {
        .name = "num_points",
        .get = (getter)PyBallTree_get_num_points,
        .doc = get_num_points_doc,
    },
    {
        .name = "leafsize",
        .get = (getter)PyBallTree_get_leafsize,
        .doc = get_leafsize_doc,
    },
    {
        .name = "sum_weight",
        .get = (getter)PyBallTree_get_sum_weight,
        .doc = get_sum_weight_doc,
    },
    {
        .name = "center",
        .get = (getter)PyBallTree_get_center,
        .doc = get_center_doc,
    },
    {
        .name = "radius",
        .get = (getter)PyBallTree_get_radius,
        .doc = get_radius_doc,
    },
    {NULL, NULL, NULL, NULL, NULL}
};

static PyMethodDef PyBallTree_methods[] = {
    // constructors
    {
        .ml_name = "from_random",
        .ml_meth = (PyCFunctionWithKeywords)PyBallTree_from_random,
        .ml_flags = METH_CLASS | METH_VARARGS | METH_KEYWORDS,
        .ml_doc = from_random_doc
    },
    {
        .ml_name = "from_file",
        .ml_meth = (PyCFunction)PyBallTree_from_file,
        .ml_flags = METH_CLASS | METH_VARARGS,
        .ml_doc = from_file_doc
    },
    // regular methods
    {
        .ml_name = "to_file",
        .ml_meth = (PyCFunction)PyBallTree_to_file,
        .ml_flags = METH_VARARGS,
        .ml_doc = to_file_doc
    },
    {
        .ml_name = "count_nodes",
        .ml_meth = (PyCFunction)PyBallTree_count_nodes,
        .ml_flags = METH_NOARGS,
        .ml_doc = count_nodes_doc
    },
    {
        .ml_name = "get_node_data",
        .ml_meth = (PyCFunction)PyBallTree_get_node_data,
        .ml_flags = METH_NOARGS,
        .ml_doc = get_node_data_doc
    },
    {
        .ml_name = "brute_radius",
        .ml_meth = (PyCFunctionWithKeywords)PyBallTree_brute_radius,
        .ml_flags = METH_VARARGS | METH_KEYWORDS,
        .ml_doc = brute_radius_doc
    },
    {
        .ml_name = "count_radius",
        .ml_meth = (PyCFunctionWithKeywords)PyBallTree_count_radius,
        .ml_flags = METH_VARARGS | METH_KEYWORDS,
        .ml_doc = count_radius_doc
    },
    {
        .ml_name = "dualcount_radius",
        .ml_meth = (PyCFunction)PyBallTree_dualcount_radius,
        .ml_flags = METH_VARARGS,
        .ml_doc = dualcount_radius_doc
    },
    {
        .ml_name = "brute_range",
        .ml_meth = (PyCFunctionWithKeywords)PyBallTree_brute_range,
        .ml_flags = METH_VARARGS | METH_KEYWORDS,
        .ml_doc = brute_range_doc
    },
    {
        .ml_name = "count_range",
        .ml_meth = (PyCFunctionWithKeywords)PyBallTree_count_range,
        .ml_flags = METH_VARARGS | METH_KEYWORDS,
        .ml_doc = count_range_doc
    },
    {
        .ml_name = "dualcount_range",
        .ml_meth = (PyCFunction)PyBallTree_dualcount_range,
        .ml_flags = METH_VARARGS,
        .ml_doc = dualcount_range_doc
    },
    {NULL, NULL, 0, NULL}
};

static PyTypeObject PyBallTreeType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "balltree.BallTree",
    .tp_doc = balltree_doc,
    .tp_basicsize = sizeof(PyBallTree),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = PyBallTree_new,
    .tp_dealloc = PyBallTree_dealloc,
    .tp_init = (initproc)PyBallTree_init,
    .tp_methods = PyBallTree_methods,
    .tp_getset = PyBallTree_getset,
    .tp_str = PyBallTree_str,
};

// module configuration ////////////////////////////////////////////////////////

static struct PyModuleDef pyballtree = {
    PyModuleDef_HEAD_INIT,
    .m_name = "balltree.balltree",
    .m_doc = "Fast balltree implementation",
    .m_size = -1,
};

PyMODINIT_FUNC PyInit_balltree(void) {
    if (PyType_Ready(&PyBallTreeType) < 0) {
        return NULL;
    }

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

    Py_INCREF(&PyBallTreeType);
    if (PyModule_AddObject(module, "BallTree", (PyObject *)&PyBallTreeType) < 0) {
        Py_DECREF(&PyBallTreeType);
        Py_DECREF(module);
        return NULL;
    }
    if (PyModule_AddIntConstant(module, "default_leafsize", DEFAULT_LEAFSIZE) < -1) {
        Py_DECREF(module);
        return NULL;
    }
    return module;
}
