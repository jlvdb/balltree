#include <stdio.h>
#include <stdlib.h>

#include "point.h"
#include "balltree.h"
#include "ballnode.h"

#define SUCCESS 0
#define FAILED  1

#define TRUE  1
#define FALSE 0

typedef struct {
    int start;
    int end;
} PointSliceSerialized;

typedef struct {
    Point center;
    double radius;
    double sum_weight;
    // the following members differ from BallNode
    int left;   // store the index of the childs in BallNodeBuffer
    int right;  // --
    PointSliceSerialized data;  // does not include the buffer pointer
} BallNodeSerialized;

typedef struct {
    int size;
    BallNodeSerialized *nodes;
} SerializedBuffer;

typedef struct {
    int size;
    int bytes;
} SectionHeader;

typedef struct {
    int leafsize;
    SectionHeader nodes;
    SectionHeader points;
} FileHeader;

static SerializedBuffer *serializedbuffer_new(int size);
static void serializedbuffer_free(SerializedBuffer *buffer);
static int bnode_serialise(SerializedBuffer *buffer, BallNode *node, int *next_index, int insertion_index);
static BallNode *bnode_deserialise(SerializedBuffer *buffer, const PointBuffer *points, int index);


static SerializedBuffer *serializedbuffer_new(int size) {
    if (size < 1) {
        fprintf(stderr, "ERROR: PointBuffer size must be positive\n");
        return NULL;
    }

    SerializedBuffer *buffer = (SerializedBuffer *)malloc(sizeof(SerializedBuffer));
    if (buffer == NULL) {
        fprintf(stderr, "ERROR: SerializedBuffer allocation failed\n");
        return NULL;
    }

    size_t n_bytes = size * sizeof(BallNodeSerialized);
    BallNodeSerialized *nodes = (BallNodeSerialized *)malloc(n_bytes);
    if (nodes == NULL) {
        fprintf(stderr, "ERROR: SerializedBuffer memory allocation failed\n");
        free(buffer);
        return NULL;
    }

    buffer->size = size;
    buffer->nodes = nodes;
    return buffer;
}

static void serializedbuffer_free(SerializedBuffer *buffer) {
    if (buffer->nodes != NULL) {
        free(buffer->nodes);
    }
    free(buffer);
}

static int bnode_serialise(SerializedBuffer *buffer, BallNode *node, int *next_index, int insertion_index) {
    if (*next_index > buffer->size) {
        return FAILED;
    }

    // use next free index in the buffer to reference left and right childs or
    // -1 to represent NULL pointer
    int index_left = (node->left != NULL) ? (*next_index)++ : -1;
    int index_right = (node->left != NULL) ? (*next_index)++ : -1;

    // serialize data and insert into output buffer
    buffer->nodes[insertion_index] = (BallNodeSerialized){
        .center = node->center,
        .radius = node->radius,
        .sum_weight = node->sum_weight,
        .left = index_left,
        .right = index_right,
        .data = (PointSliceSerialized){
            .start = node->data.start,
            .end = node->data.end,
        },
    };

    // handle the childs if any
    int state;
    if (index_left != -1) {
        state = bnode_serialise(buffer, node->left, next_index, index_left);
        if (state == FAILED) {
            return FAILED;
        }
    }
    if (index_right != -1) {
        state = bnode_serialise(buffer, node->right, next_index, index_right);
        if (state == FAILED) {
            return FAILED;
        }
    }
    return SUCCESS;
}

static BallNode *bnode_deserialise(SerializedBuffer *buffer, const PointBuffer *points, int index) {
    if (index >= buffer->size) {
        fprintf(stderr, "ERROR: node index exceeds node array size\n");
        return NULL;
    }
    BallNodeSerialized *serialized = buffer->nodes + index;
    if (serialized->data.end > points->size) {
        fprintf(stderr, "ERROR: point buffer does not match data slice expected by node\n");
        return NULL;
    }

    // reconstruct the original BallNode
    BallNode *node = (BallNode *)malloc(sizeof(BallNode));
    if (node == NULL) {
        fprintf(stderr, "ERROR: failed to allocate a new node\n");
        return NULL;
    }
    node->center = serialized->center;
    node->radius = serialized->radius;
    node->sum_weight = serialized->sum_weight;
    node->data = (PointSlice){
        .start = serialized->data.start,
        .end = serialized->data.end,
        .points = points->points,
    };
    node->left = NULL;
    node->right = NULL;

    // attach any missing childs
    if (serialized->left != -1) {
        node->left = bnode_deserialise(buffer, points, serialized->left);
        if (node->left == NULL) {
            bnode_free(node);
            return NULL;
        }
    }
    if (serialized->right != -1) {
        node->right = bnode_deserialise(buffer, points, serialized->right);
        if (node->right == NULL) {
            bnode_free(node);
            return NULL;
        }
    }
    return node;
}

int balltree_to_file(const BallTree *tree, const char *path) {
    size_t elements_written;
    FILE *file = fopen(path, "wb");
    if (file == NULL) {
        fprintf(stderr, "ERROR: failed to open file\n");
        return FAILED;
    }

    // write the header
    int num_nodes = balltree_count_nodes(tree);
    int num_points = tree->data.size;
    FileHeader header = {
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
    elements_written = fwrite(&header, sizeof(FileHeader), 1, file);
    if (elements_written != 1) {
        fprintf(stderr, "ERROR: failed to write header\n");
        fclose(file);
        return FAILED;
    }

    // append serialised point buffer
    elements_written = fwrite(tree->data.points, sizeof(Point), header.points.size, file);
    if (elements_written != header.points.size) {
        fprintf(stderr, "ERROR: failed to write data points\n");
        fclose(file);
        return FAILED;
    }

    // serialise nodes into buffer
    SerializedBuffer *node_buffer = serializedbuffer_new(header.nodes.size);
    if (node_buffer == NULL) {
        fclose(file);
        return FAILED;
    }
    int index_tracker = 1;  // 0 is already reserved for root node
    int success = bnode_serialise(node_buffer, tree->root, &index_tracker, 0);
    if (success == FAILED) {
        serializedbuffer_free(node_buffer);
        fclose(file);
        return FAILED;
    }

    // append serialised node buffer
    elements_written = fwrite(node_buffer->nodes, sizeof(BallNodeSerialized), header.nodes.size, file);
    serializedbuffer_free(node_buffer);
    if (elements_written != header.nodes.size) {
        fprintf(stderr, "ERROR: failed to write node data\n");
        fclose(file);
        return FAILED;
    }

    if (fflush(file) == EOF) {
        fprintf(stderr, "ERROR: failed to flush file\n");
        fclose(file);
        return FAILED;
    }

    fclose(file);
    return SUCCESS;
}

BallTree* balltree_from_file(const char *path) {
    size_t elements_read;
    FILE *file = fopen(path, "rb");
    if (file == NULL) {
        fprintf(stderr, "ERROR: failed to open file\n");
        return NULL;
    }

    // read header
    FileHeader header;
    elements_read = fread(&header, sizeof(FileHeader), 1, file);
    if (elements_read != 1) {
        fprintf(stderr, "ERROR: failed to read header\n");
        goto error_close_file;
    }

    // allocate all memory required
    BallTree *tree = (BallTree *)malloc(sizeof(BallTree));
    if (tree == NULL) {
        fprintf(stderr, "ERROR: failed to allocate tree\n");
        goto error_close_file;
    }
    Point *points = (Point *)malloc(header.points.bytes);
    if (points == NULL) {
        fprintf(stderr, "ERROR: failed to allocate point buffer\n");
        goto error_dealloc_tree;
    }
    SerializedBuffer *node_buffer = serializedbuffer_new(header.nodes.size);
    if (node_buffer == NULL) {
        fprintf(stderr, "ERROR: failed to allocate node data\n");
        goto error_dealloc_points;
    }

    // read point buffer
    elements_read = fread(points, sizeof(Point), header.points.size, file);
    if (elements_read != header.points.size) {
        fprintf(stderr, "ERROR: failed to read data points\n");
        goto error_dealloc_nodes;
    }

    // read node data
    elements_read = fread(node_buffer->nodes, sizeof(BallNodeSerialized), header.nodes.size, file);
    if (elements_read != header.nodes.size) {
        fprintf(stderr, "ERROR: failed to read node data\n");
        goto error_dealloc_nodes;
    }

    // populate the tree
    tree->leafsize = header.leafsize;
    tree->data = (PointBuffer){
        .size = header.points.size,
        .points = points,
    };
    tree->root = bnode_deserialise(node_buffer, &tree->data, 0);
    if (tree->root == NULL) {
        fprintf(stderr, "ERROR: failed to reconstruct tree\n");
        goto error_dealloc_nodes;
    }

    serializedbuffer_free(node_buffer);
    fclose(file);
    return tree;

    // alternative exit route which cleans up buffers
error_dealloc_nodes:
    serializedbuffer_free(node_buffer);
error_dealloc_points:
    free(points);
error_dealloc_tree:
    free(tree);
error_close_file:
    fclose(file);
    return NULL;
}
