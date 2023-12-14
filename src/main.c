#include <stdio.h>

#include "point.h"
#include "balltree.h"


int main(int argc, char** argv)
{
    double xi, yi, zi;
    int n_records = 0;
    struct PointBuffer points = pointbuffer_create(100);
    if (points.size < 0) {
        return 1;
    }

    FILE *file = fopen("testing/points.txt", "r");
    if (file == NULL) {
        perror("Error opening file");
        free(points.points);
        return 1;
    }
    // load the data from file while dynamically growing the point buffer
    while (fscanf(file, "%lf %lf %lf", &xi, &yi, &zi) == 3) {
        if (n_records == points.size) {
            if (!pointbuffer_resize(&points, points.size * 2)) {
                free(points.points);
                fclose(file);
                return 1;
            }
        }
        struct Point point = {xi, yi, zi};
        points.points[n_records] = point;
        n_records++;
    }
    fclose(file);
    // truncate the buffer to the actual size
    if (!pointbuffer_resize(&points, n_records)) {
        free(points.points);
        return 1;
    }

    // build the tree and print it
    struct PointSlice slice = pointslice_from_buffer(points);
    struct BallTree *tree = balltree_build(&slice, 6);
    free(points.points);
    if (!tree) {
        printf("ERROR: tree building failed\n");
        return 1;
    }
    balltree_print(tree);
    balltree_free(tree);
    return 0;
}
