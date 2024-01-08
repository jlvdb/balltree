#include <stdio.h>
#include <time.h>
#include <locale.h>

#include "point.h"
#include "balltree.h"
#include "balltree_macros.h"

PointBuffer *load_data_from_file();

int main(int argc, char** argv) {
    Point query_point = {0.0, 0.0, 0.0, 0.5};
    double query_radius = 0.2;
    int leafsize = 20;

    PointBuffer *buffer;
    BallTree *tree;
    double count;

    clock_t time;
    double elapsed;
    setlocale(LC_NUMERIC, "");

    // read the input file contents, show the elapsed time
    time = clock();
    buffer = load_data_from_file();
    if (buffer == NULL) {
        return 1;
    }
    elapsed = (double)(clock() - time) / CLOCKS_PER_SEC;
    int n_records = buffer->size;
    printf("read %'d records in %.3lf sec\n", n_records, elapsed);

    // build the ball tree, show the elapsed time
    time = clock();
    tree = balltree_build_leafsize(buffer, leafsize);
    if (tree == NULL) {
        ptbuf_free(buffer);
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
    ptbuf_free(buffer);  // no longer needed

    // query with itself, show the elapsed time
    time = clock();
    count = balltree_dualcount_radius(tree, tree, query_radius);
    elapsed = (double)(clock() - time) / CLOCKS_PER_SEC;
    printf("self found %9.0lf pairs in %7.3lf sec\n", count, elapsed);

    // dump and restore
    time = clock();
    if (balltree_to_file(tree, "testing/tree.dump") != 0) {
        balltree_free(tree);
        return 1;
    }
    elapsed = (double)(clock() - time) / CLOCKS_PER_SEC * 1000.0;
    printf("dumped   in %7.3lf ms\n", elapsed);
    time = clock();
    BallTree *tree2 = balltree_from_file("testing/tree.dump");
    if (tree2 == NULL) {
        balltree_free(tree); 
        return 1;  
    }
    elapsed = (double)(clock() - time) / CLOCKS_PER_SEC * 1000.0;
    printf("restored in %7.3lf ms\n", elapsed);

    balltree_free(tree);
    count = balltree_dualcount_radius(tree2, tree2, query_radius);
    printf("self found %9.0lf pairs\n", count);

    balltree_free(tree2);
    return 0;
}

PointBuffer *load_data_from_file() {
    PointBuffer *buffer = ptbuf_new(256);
    if (buffer == NULL) {
        PRINT_ERR_MSG("memory allocation failed\n");
        ptbuf_free(buffer);
        return NULL;
    }

    static const char filepath[] = "testing/points.txt";
    FILE *file = fopen(filepath, "r");
    if (file == NULL) {
        PRINT_ERR_MSG("failed to open file: %s\n", filepath);
        ptbuf_free(buffer);
        return NULL;
    }

    int n_records = 0;
    Point point = {0.0, 0.0, 0.0, 0.5};
    while (fscanf(file, "%lf %lf %lf", &point.x, &point.y, &point.z) == 3) {
        if (n_records == buffer->size) {
            if (ptbuf_resize(buffer, buffer->size * 2) != 0) {
                PRINT_ERR_MSG("failed to expand buffer\n");
                ptbuf_free(buffer);
                fclose(file);
                return NULL;
            }
        }
        buffer->points[n_records] = point;
        ++n_records;
    }
    fclose(file);

    if (n_records == 0) {
        PRINT_ERR_MSG("could not read any records from file\n");
        ptbuf_free(buffer);
        return NULL;
    }
    ptbuf_resize(buffer, n_records);
    return buffer;
}
