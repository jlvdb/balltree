#ifndef ERROR_HANDLING_H
#define ERROR_HANDLING_H

#include <stdio.h>

#define PRINT_ERR_MSG(format, ...) do { \
    fprintf(stderr, "ERROR: "); \
    fprintf(stderr, format, ##__VA_ARGS__); \
    fprintf(stderr, "\n"); \
} while (0)

#endif  // ERROR_HANDLING_H
