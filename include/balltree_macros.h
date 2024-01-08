#ifndef BALLTREE_MACROS_H
#define BALLTREE_MACROS_H

#include <stdio.h>

#define PRINT_ERR_MSG(format, ...) do { \
    fprintf(stderr, "ERROR: "); \
    fprintf(stderr, format, ##__VA_ARGS__); \
    fprintf(stderr, "\n"); \
} while (0)

#define BTR_SUCCESS 0
#define BTR_FAILED  1

#endif  // BALLTREE_MACROS_H
