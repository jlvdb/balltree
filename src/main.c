#include <stdio.h>
#include <time.h>
#include <locale.h>

#include "point.h"
#include "balltree.h"


int main(int argc, char** argv)
{
    clock_t start_time, end_time;
    double difftime;
    setlocale(LC_NUMERIC, "");

    double xi, yi, zi;
    int n_records = 0;
    struct PointBuffer points = pointbuffer_create(100);
    if (points.size < 0) {
        perror("memory allocation failed");
        free(points.points);
        return 1;
    }

    start_time = clock();
    FILE *file = fopen("testing/points.txt", "r");
    if (file == NULL) {
        perror("failed to open file");
        free(points.points);
        return 1;
    }
    // load the data from file while dynamically growing the point buffer
    while (fscanf(file, "%lf %lf %lf", &xi, &yi, &zi) == 3) {
        if (n_records == points.size) {
            if (!pointbuffer_resize(&points, points.size * 2)) {
                perror("memory allocation failed");
                free(points.points);
                fclose(file);
                return 1;
            }
        }
        struct Point point = {xi, yi, zi};
        points.points[n_records] = point;
        ++n_records;
    }
    fclose(file);
    // truncate the buffer to the actual size
    if (n_records == 0) {
        perror("did not read any records from file");
        free(points.points);
        return 1;
    }
    if (!pointbuffer_resize(&points, n_records)) {
        perror("memory reallocation failed");
        free(points.points);
        return 1;
    }
    end_time = clock();
    difftime = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("read %'d records in %.3lf sec\n", n_records, difftime);


    struct Point qpoint = {0.0, 0.0, 0.0};
    double radius = 0.1;
    int leafsize = 4096;


    start_time = clock();
    struct PointSlice slice = pointslice_from_buffer(points);
    struct BallTree *tree = balltree_build(&slice, leafsize);
    if (!tree) {
        perror("tree building failed");
        free(points.points);
        return 1;
    }
    end_time = clock();
    difftime = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("built tree in %.3lf sec\n", difftime);

    start_time = clock();
    double count = balltree_count_radius(tree, &qpoint, radius);
    end_time = clock();
    difftime = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000.0;
    printf("found %.0lf pairs in %.3lf ms\n", count, difftime);

    start_time = clock();
    count = count_within_radius(&points, &qpoint, radius);
    end_time = clock();
    difftime = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000.0;
    printf("brute %.0lf pairs in %.3lf ms\n", count, difftime);

    balltree_free(tree);
    free(points.points);
    return 0;
}
