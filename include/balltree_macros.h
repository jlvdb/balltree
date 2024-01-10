#ifndef BALLTREE_MACROS_H
#define BALLTREE_MACROS_H

#ifdef SET_PYERR_STRING
#include <Python.h>
#endif

#include <stdio.h>

// behaviour in python environment: format and set the python error string
#ifdef SET_PYERR_STRING
#define EMIT_ERR_MSG(ErrorName, msg_format, ...) \
    do { \
        char *msg_buffer; \
        if (asprintf(&msg_buffer, msg_format, ##__VA_ARGS__) != -1) { \
            PyErr_SetString(PyExc_##ErrorName, msg_buffer); \
            free(msg_buffer); \
        } else { \
            PyErr_SetString(PyExc_MemoryError, "failed to allocate memory for error message"); \
        } \
    } while (0)
#else
// default behaviour: format and print an error message to stderr
#define EMIT_ERR_MSG(ErrorName, msg_format, ...) \
    do { \
        fprintf(stderr, #ErrorName ": "); \
        fprintf(stderr, msg_format, ##__VA_ARGS__); \
        fprintf(stderr, "\n"); \
    } while (0)
#endif

#define BTR_SUCCESS 0
#define BTR_FAILED  1

#endif  // BALLTREE_MACROS_H
