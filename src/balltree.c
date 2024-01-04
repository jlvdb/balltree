#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "point.h"
#include "ballnode.h"
#include "balltree.h"

#define SUCCESS 1
#define FAILED  0

#define DEFAULT_LEAFSIZE 20

struct SectionHeader {
    int size;
    int bytes;
};

struct FileHeader {
    int leafsize;
    struct SectionHeader nodes;
    struct SectionHeader points;
};

BallTree* balltree_build_leafsize(const PointBuffer *buffer, int leafsize) {
    BallTree *tree = (BallTree*)malloc(sizeof(BallTree));
    if (!tree) {
        return NULL;
    }
    tree->leafsize = leafsize;

    // copy the input data buffer, since it will be partitioned while build the tree
    PointBuffer data = {buffer->size, NULL};
    size_t n_bytes = buffer->size * sizeof(Point);
    data.points = (Point*)malloc(n_bytes);
    if (!data.points) {
        free(tree);
        return NULL;
    }
    memcpy(data.points, buffer->points, n_bytes);
    tree->data = data;

    PointSlice *slice = pointslice_from_buffer(&data);
    tree->root = ballnode_build_recursive(slice, leafsize);
    if (!tree->root) {
        free(data.points);
        free(tree);
        return NULL;
    }
    return tree;
}

BallTree* balltree_build(const PointBuffer *buffer) {
    return balltree_build_leafsize(buffer, DEFAULT_LEAFSIZE);
}

void balltree_free(BallTree *tree) {
    if (tree->data.points) {
        free(tree->data.points);
    }
    if (tree->root) {
        free(tree->root);
    }
    free(tree);
    tree = NULL;
}

int balltree_count_nodes(const BallTree *tree) {
    return ballnode_count_nodes(tree->root);
}

int balltree_to_file(const BallTree *tree, const char *path) {
    size_t elements_written;
    FILE *file = fopen(path, "wb");
    if (!file) {
        fprintf(stderr, "ERROR: failed to open file\n");
        return FAILED;
    }

    // write the header
    int num_nodes = balltree_count_nodes(tree);
    int num_points = tree->data.size;
    struct FileHeader header = {
        .leafsize = tree->leafsize,
        .nodes = {
            .size = num_nodes,
            .bytes = num_nodes * sizeof(BallNodeSerialized)
        },
        .points = {
            .size = num_points,
            .bytes = num_points * sizeof(Point)
        },
    };
    elements_written = fwrite(&header, sizeof(struct FileHeader), 1, file);
    if (elements_written != 1) {
        fprintf(stderr, "ERROR: failed to write header\n");
        goto close_file;
    }

    // append serialised point buffer
    elements_written = fwrite(tree->data.points, sizeof(Point), header.points.size, file);
    if (elements_written != header.points.size) {
        fprintf(stderr, "ERROR: failed to write data points\n");
        goto close_file;
    }

    // serialise nodes
    int index_tracker = 1;  // 0 is already reserved for root node
    BallNodeBuffer node_buffer = {
        .size = header.nodes.size,
        .next_free = &index_tracker,
    };
    node_buffer.buffer = (BallNodeSerialized*)malloc(header.nodes.bytes);
    if (!node_buffer.buffer) {
        fprintf(stderr, "ERROR: failed to allocate memory for serialized node data\n");
        goto close_file;
    }
    int success = ballnode_serialise_recursive(node_buffer, tree->root, 0);
    if (!success) {
        goto dealloc_nodes;
    }

    // append serialised nodes
    elements_written = fwrite(node_buffer.buffer, sizeof(BallNodeSerialized), header.nodes.size, file);
    free(node_buffer.buffer);
    if (elements_written != header.nodes.size) {
        fprintf(stderr, "ERROR: failed to write node data\n");
        goto close_file;
    }

    if (fflush(file) == EOF) {
        fprintf(stderr, "ERROR: failed to flush file\n");
        goto close_file;
    }
    fclose(file);

    return SUCCESS;

    // alternative exit route which cleans up buffers
dealloc_nodes:
    free(node_buffer.buffer);
close_file:
    fclose(file);
    return FAILED;
}

BallTree* balltree_from_file(const char *path) {
    size_t elements_read;
    FILE *file = fopen(path, "rb");
    if (file == NULL) {
        fprintf(stderr, "ERROR: failed to open file\n");
        return NULL;
    }

    // read header
    struct FileHeader header;
    elements_read = fread(&header, sizeof(struct FileHeader), 1, file);
    if (elements_read != 1) {
        fprintf(stderr, "ERROR: failed to read header\n");
        goto close_file;
    }

    // allocate all memory required
    BallTree *tree = (BallTree*)malloc(sizeof(BallTree));
    if (!tree) {
        fprintf(stderr, "ERROR: failed to allocate tree\n");
        goto close_file;
    }
    Point *points = (Point*)malloc(header.points.bytes);
    if (!points) {
        fprintf(stderr, "ERROR: failed to allocate point buffer\n");
        goto dealloc_tree;
    }
    BallNodeSerialized *node_buffer = (BallNodeSerialized *)malloc(header.nodes.bytes);
    if (!node_buffer) {
        fprintf(stderr, "ERROR: failed to allocate node data\n");
        goto dealloc_points;
    }

    // read point buffer
    elements_read = fread(points, sizeof(Point), header.points.size, file);
    if (elements_read != header.points.size) {
        fprintf(stderr, "ERROR: failed to read data points\n");
        goto dealloc_nodes;
    }

    // read node data
    elements_read = fread(node_buffer, sizeof(BallNodeSerialized), header.nodes.size, file);
    if (elements_read != header.nodes.size) {
        fprintf(stderr, "ERROR: failed to read node data\n");
        goto dealloc_nodes;
    }

    // populate the tree
    tree->leafsize = header.leafsize;
    tree->data = (PointBuffer){
        .size = header.points.size,
        .points = points,
    };
    tree->root = ballnode_deserialise_recursive(node_buffer, header.nodes.size, &tree->data, 0);
    if (!tree->root) {
        fprintf(stderr, "ERROR: failed to reconstruct tree\n");
        goto dealloc_nodes;
    }

    fclose(file);
    return tree;

    // alternative exit route which cleans up buffers
dealloc_nodes:
    free(node_buffer);
dealloc_points:
    free(points);
dealloc_tree:
    free(tree);
close_file:
    fclose(file);
    return NULL;
}

double balltree_count_radius(
    const BallTree *tree,
    const Point *point,
    double radius
) {
    return ballnode_count_radius(tree->root, point, radius);
}

double balltree_count_range(
    const BallTree *tree,
    const Point *point,
    double rmin,
    double rmax
) {
    return ballnode_count_range(tree->root, point, rmin, rmax);
}

double balltree_dualcount_radius(
    const BallTree *tree1,
    const BallTree *tree2,
    double radius
) {
    return ballnode_dualcount_radius(tree1->root, tree2->root, radius);
}
