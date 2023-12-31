#include <stdio.h>
#include <time.h>
#include <locale.h>

#include "point.h"
#include "balltree.h"

struct PointBuffer* load_data_from_file()
{
    struct PointBuffer *buffer = pointbuffer_create(100);
    if (!buffer) {
        fprintf(stderr, "ERROR: memory allocation failed\n");
        pointbuffer_free(buffer);
        return NULL;
    }

    static const char filepath[] = "testing/points.txt";
    FILE *file = fopen(filepath, "r");
    if (file == NULL) {
        fprintf(stderr, "ERROR: failed to open file: %s\n", filepath);
        pointbuffer_free(buffer);
        return NULL;
    }

    int n_records = 0;
    struct Point point = {0.0, 0.0, 0.0, 0.5};
    while (fscanf(file, "%lf %lf %lf", &point.x, &point.y, &point.z) == 3) {
        if (n_records == buffer->size) {
            if (!pointbuffer_resize(buffer, buffer->size * 2)) {
                fprintf(stderr, "ERROR: failed to expand buffer\n");
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
        fprintf(stderr, "ERROR: could not read any records from file\n");
        pointbuffer_free(buffer);
        return NULL;
    }
    if (!pointbuffer_resize(buffer, n_records)) {
        fprintf(stderr, "ERROR: memory reallocation failed\n");
        pointbuffer_free(buffer);
        return NULL;
    }
    return buffer;
}

int main(int argc, char** argv)
{
    struct Point query_point = {0.0, 0.0, 0.0, 0.5};
    double query_radius = 0.2;
    int leafsize = 20;

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
    tree = balltree_build_leafsize(buffer, leafsize);
    if (!tree) {
        return 1;
    }
    elapsed = (double)(clock() - time) / CLOCKS_PER_SEC;
    printf("built tree in %.3lf sec\n", elapsed);
    printf("radius=%.3lf\n", query_radius);

    // query point at fixed radius, show the elapsed time
    int imax = 1;
    while (imax <= 100) {
        time = clock();
        for (int i = 0; i < imax ; ++i)
            count = balltree_count_radius(tree, &query_point, query_radius);
        elapsed = (double)(clock() - time) / CLOCKS_PER_SEC * 1000.0;
        printf("%3dx found %9.0lf pairs in %7.3lf ms\n", imax, count, elapsed);
        imax *= 10;
    }

    // query all points, show the elapsed time
    time = clock();
    count = 0.0;
    for (int i = 0; i < buffer->size; ++i) {
        count += balltree_count_radius(tree, buffer->points + i, query_radius);
    }
    elapsed = (double)(clock() - time) / CLOCKS_PER_SEC;
    printf(" all found %9.0lf pairs in %7.3lf sec\n", count, elapsed);

    // query with itself, show the elapsed time
    time = clock();
    count = balltree_dualcount_radius(tree, tree, query_radius);
    elapsed = (double)(clock() - time) / CLOCKS_PER_SEC;
    printf("self found %9.0lf pairs in %7.3lf sec\n", count, elapsed);

    // dump and restore
    time = clock();
    if (!balltree_to_file(tree, "testing/tree.dump")) {
        pointbuffer_free(buffer);
        balltree_free(tree);
        return 1;
    }
    elapsed = (double)(clock() - time) / CLOCKS_PER_SEC * 1000.0;
    printf("dumped   in %7.3lf ms\n", elapsed);
    time = clock();
    struct BallTree *tree2 = balltree_from_file("testing/tree.dump");
    if (!tree2) {
        pointbuffer_free(buffer);
        balltree_free(tree); 
        return 1;  
    }
    elapsed = (double)(clock() - time) / CLOCKS_PER_SEC * 1000.0;
    printf("restored in %7.3lf ms\n", elapsed);

    pointbuffer_free(buffer);
    balltree_free(tree2);
    balltree_free(tree);
    return 0;
}
