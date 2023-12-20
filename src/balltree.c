#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "point.h"
#include "ballnode.h"
#include "balltree.h"

#define SUCCESS 1
#define FAILED  0

#define DEFAULT_LEAFSIZE 20

struct FileHeader {
    int leafsize;
    int num_nodes;
    int num_points;
};

struct BallTree* balltree_build_leafsize(const struct PointBuffer *buffer, int leafsize)
{
    struct BallTree *tree = (struct BallTree*)malloc(sizeof(struct BallTree));
    if (!tree) {
        return NULL;
    }
    tree->leafsize = leafsize;

    // copy the input data buffer, since it will be partitioned while build the tree
    struct PointBuffer data = {buffer->size, NULL};
    size_t n_bytes = buffer->size * sizeof(struct Point);
    data.points = (struct Point*)malloc(n_bytes);
    if (!data.points) {
        free(tree);
        return NULL;
    }
    memcpy(data.points, buffer->points, n_bytes);
    tree->data = data;

    struct PointSlice *slice = pointslice_from_buffer(&data);
    tree->root = ballnode_build_recursive(slice, leafsize);
    if (!tree->root) {
        free(data.points);
        free(tree);
        return NULL;
    }
    return tree;
}

struct BallTree* balltree_build(const struct PointBuffer *buffer)
{
    return balltree_build_leafsize(buffer, DEFAULT_LEAFSIZE);
}

void balltree_free(struct BallTree *tree)
{
    if (tree->data.points) {
        free(tree->data.points);
    }
    if (tree->root) {
        free(tree->root);
    }
    free(tree);
    tree = NULL;
}

int balltree_count_nodes(const struct BallTree *tree)
{
    return ballnode_count_nodes(tree->root);
}

int balltree_to_file(const struct BallTree *tree, const char *path)
{
    size_t elements_written;
    FILE *file = fopen(path, "wb");
    if (!file) {
        fprintf(stderr, "ERROR: failed to open file\n");
        return FAILED;
    }

    // write the header
    struct FileHeader header = {
        .leafsize = tree->leafsize,
        .num_nodes = balltree_count_nodes(tree),
        .num_points = tree->data.size,
    };
    elements_written = fwrite(&header, sizeof(struct FileHeader), 1, file);
    if (elements_written != 1) {
        fprintf(stderr, "ERROR: failed to write header\n");
        goto close_file;
    }

    // append serialised point buffer
    elements_written = fwrite(tree->data.points, sizeof(struct Point), header.num_points, file);
    if (elements_written != header.num_points) {
        fprintf(stderr, "ERROR: failed to write data points\n");
        goto close_file;
    }

    // serialise nodes
    int index_tracker = 1;  // 0 is already reserved for root node
    struct BallNodeBuffer node_buffer = {
        .size = header.num_nodes,
        .next_free = &index_tracker,
    };
    size_t n_bytes = header.num_nodes * sizeof(struct BallNodeSerialized);
    node_buffer.buffer = (struct BallNodeSerialized*)malloc(n_bytes);
    if (!node_buffer.buffer) {
        fprintf(stderr, "ERROR: failed to allocate memory for serialized node data\n");
        goto close_file;
    }
    int success = ballnode_serialise_recursive(node_buffer, tree->root, 0);
    if (!success) {
        goto dealloc_nodes;
    }

    // append serialised nodes
    elements_written = fwrite(node_buffer.buffer, sizeof(struct BallNodeSerialized), header.num_nodes, file);
    free(node_buffer.buffer);
    if (elements_written != header.num_nodes) {
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

struct BallTree* balltree_from_file(const char *path)
{
    size_t elements_read;
    FILE *file = fopen(path, "rb");
    if (file == NULL) {
        fprintf(stderr, "ERROR: failed to open file\n");
        goto close_file;
    }

    // read header
    struct FileHeader header;
    elements_read = fread(&header, sizeof(struct FileHeader), 1, file);
    if (elements_read != 1) {
        fprintf(stderr, "ERROR: failed to read header\n");
        goto close_file;
    }

    struct BallTree *tree = (struct BallTree*)malloc(sizeof(struct BallTree));
    if (!tree) {
        fprintf(stderr, "ERROR: failed to allocate tree\n");
        goto dealloc_tree;
    }
    struct Point *points = (struct Point*)malloc(header.num_points * sizeof(struct Point));
    if (!points) {
        fprintf(stderr, "ERROR: failed to allocate point buffer\n");
        goto dealloc_points;
    }
    tree->data = (struct PointBuffer){
        .size = header.num_points,
        .points = points,
    };

    // read point buffer
    elements_read = fread(points, sizeof(struct Point), header.num_points, file);
    if (elements_read != header.num_points) {
        fprintf(stderr, "ERROR: failed to read data points\n");
        goto dealloc_points;
    }

    // read node data
    size_t n_bytes = header.num_nodes * sizeof(struct BallNodeSerialized);
    struct BallNodeSerialized *node_buffer = (struct BallNodeSerialized *)malloc(n_bytes);
    if (!node_buffer) {
        fprintf(stderr, "ERROR: failed to allocate node data\n");
        goto dealloc_nodes;
    }
    elements_read = fread(node_buffer, sizeof(struct BallNodeSerialized), header.num_nodes, file);
    if (elements_read != header.num_nodes) {
        fprintf(stderr, "ERROR: failed to read node data\n");
        goto dealloc_points;
    }

    tree->root = ballnode_deserialise_recursive(node_buffer, header.num_nodes, &tree->data, 0);
    if (!tree->root) {
        fprintf(stderr, "ERROR: failed to reconstruct tree\n");
        goto dealloc_nodes;
    }

    return tree;

    // alternative exit route which cleans up buffers
dealloc_nodes:
    free(node_buffer);
dealloc_points:
    free(points);
dealloc_tree:
    balltree_free(tree);
close_file:
    fclose(file);
    return NULL;
}

double balltree_count_radius(const struct BallTree *tree, const struct Point *point, double radius)
{
    return ballnode_count_radius(tree->root, point, radius);
}

double balltree_count_range(const struct BallTree *tree, const struct Point *point, double rmin, double rmax)
{
    return ballnode_count_range(tree->root, point, rmin, rmax);
}

double balltree_dualcount_radius(const struct BallTree *tree1, const struct BallTree *tree2, double radius)
{
    return ballnode_dualcount_radius(tree1->root, tree2->root, radius);
}
