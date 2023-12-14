#include <stdio.h>

#include "point.h"
#include "balltree.h"


int main(int argc, char** argv)
{
    const int N = 9;
    struct Point points[N] = {
        {1.1, 2.1, 3.1},
        {4.1, 5.1, 6.1},
        {7.1, 8.1, 9.1},
        {10.1, 11.1, 12.1},

        {1.2, 2.2, 3.2},
        {4.2, 5.2, 6.2},
        {7.2, 8.2, 9.2},
        {10.2, 11.2, 12.2},

        {10.3, 11.3, 12.3},
    };
    struct PointSlice slice = {
        .start = 0,
        .end = N,
        .points = points,
    };
    struct BallTree *tree = balltree_build(&slice, 3);
    if (!tree) {
        printf("ERROR: tree building failed\n");
        balltree_print(tree);
        return 1;
    }
    balltree_print(tree);
    balltree_free(tree);
    return 0;
}
