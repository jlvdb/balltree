#include <stdio.h>
#include <time.h>
#include <locale.h>

#include "point.h"
#include "balltree.h"

struct PointBuffer* load_data_from_file()
{
    struct PointBuffer *buffer = pointbuffer_create(100);
    if (!buffer) {
        fprintf(stderr, "ERROR: memory allocation failed");
        pointbuffer_free(buffer);
        return NULL;
    }

    static const char filepath[] = "testing/points.txt";
    FILE *file = fopen(filepath, "r");
    if (file == NULL) {
        fprintf(stderr, "ERROR: failed to open file: %s", filepath);
        pointbuffer_free(buffer);
        return NULL;
    }

    int n_records = 0;
    struct Point point;
    while (fscanf(file, "%lf %lf %lf", &point.x, &point.y, &point.z) == 3) {
        if (n_records == buffer->size) {
            if (!pointbuffer_resize(buffer, buffer->size * 2)) {
                fprintf(stderr, "ERROR: failed to expand buffer");
                pointbuffer_free(buffer);
                fclose(file);
                return NULL;
            }
        }
        buffer->points[n_records] = point;
        ++n_records;
    }
    fclose(file);

    if (n_records == 0) {
        fprintf(stderr, "ERROR: could not read any records from file");
        pointbuffer_free(buffer);
        return NULL;
    }
    if (!pointbuffer_resize(buffer, n_records)) {
        fprintf(stderr, "ERROR: memory reallocation failed");
        pointbuffer_free(buffer);
        return NULL;
    }
    return buffer;
}

int main(int argc, char** argv)
{
    struct Point query_point = {0.0, 0.0, 0.0};
    double query_radius = 0.1;
    int leafsize = 4096;

    struct PointBuffer *buffer;
    struct BallTree *tree;
    double count;

    clock_t time;
    double elapsed;
    setlocale(LC_NUMERIC, "");

    // read the input file contents, show the elapsed time
    time = clock();
    buffer = load_data_from_file();
    if (!buffer) {
        return 1;
    }
    elapsed = (double)(clock() - time) / CLOCKS_PER_SEC;
    int n_records = buffer->size;
    printf("read %'d records in %.3lf sec\n", n_records, elapsed);

    // build the ball tree, show the elapsed time
    time = clock();
    tree = balltree_build(buffer, leafsize);
    if (!tree) {
        pointbuffer_free(buffer);
        return 1;
    }
    elapsed = (double)(clock() - time) / CLOCKS_PER_SEC;
    printf("built tree in %.3lf sec\n", elapsed);

    // query point at fixed radius, show the elapsed time
    time = clock();
    count = balltree_count_radius(tree, &query_point, query_radius);
    elapsed = (double)(clock() - time) / CLOCKS_PER_SEC * 1000.0;
    printf("found %.0lf pairs in %.3lf ms\n", count, elapsed);

    // bruteforce query point at fixed radius, show the elapsed time
    time = clock();
    count = count_within_radius(buffer, &query_point, query_radius);
    elapsed = (double)(clock() - time) / CLOCKS_PER_SEC * 1000.0;
    printf("brute %.0lf pairs in %.3lf ms\n", count, elapsed);

    balltree_free(tree);
    pointbuffer_free(buffer);
    return 0;
}
