#ifndef BALLTREE_MACROS_H
#define BALLTREE_MACROS_H

#ifdef SET_PYERR_STRING
#include <Python.h>
#endif

#include <stdio.h>

// behaviour in python environment: format and set the python error string
#ifdef SET_PYERR_STRING
#define EMIT_ERR_MSG(ErrorName, msg_format, ...) \
    PyErr_Format(PyExc_##ErrorName, msg_format, ##__VA_ARGS__)
#else
// default behaviour: format and print an error message to stderr
#define EMIT_ERR_MSG(ErrorName, msg_format, ...) \
    fprintf(stderr, #ErrorName ": " msg_format "\n", ##__VA_ARGS__)
#endif

#define BTR_SUCCESS 0
#define BTR_FAILED  1

#endif  // BALLTREE_MACROS_H
