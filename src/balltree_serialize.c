#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "point.h"
#include "balltree.h"
#include "ballnode.h"
#include "balltree_macros.h"

static int ptbuf_write(const PointBuffer *buffer, FILE *file);
static PointBuffer *ptbuf_read(int n_items, FILE *file);

typedef struct {
    int size;
    BallNode *nodes;
    int next_free;
} BNodeBuffer;

static BNodeBuffer *bnodebuffer_new(int size);
static void bnodebuffer_free(BNodeBuffer *buffer);
static int bnodebuffer_get_next_free(BNodeBuffer *buffer);
static int bnodebuffer_write(const BNodeBuffer *buffer, FILE *file);
static BNodeBuffer *bnodebuffer_read(int n_items, FILE *file);

typedef struct {
    int n_items;
    int itemsize;
} SectionHeader;

typedef struct {
    int leafsize;
    SectionHeader points;
    SectionHeader nodes;
} FileHeader;

static FileHeader *fileheader_new(const BallTree *tree);
static int fileheader_write(const FileHeader *header, FILE *file);
static FileHeader *fileheader_read(FILE *file);

static inline intptr_t child_ptr_substitute(const BallNode *child, BNodeBuffer *buffer);
static int bnode_serialise(const BallNode *node, BNodeBuffer *buffer, int buf_idx);
static BallNode *bnode_deserialise(const BNodeBuffer *buffer, Point *points, int buf_idx);


static int ptbuf_write(const PointBuffer *buffer, FILE *file) {
    size_t n_items = (size_t)buffer->size;
    size_t n_written = fwrite(buffer->points, sizeof(Point), n_items, file);
    if (n_written != n_items) {
        EMIT_ERR_MSG(IOError, "failed to write %zu data points", n_items);
        return BTR_FAILED;
    }
    return BTR_SUCCESS;
}

static PointBuffer *ptbuf_read(int n_items, FILE *file) {
    PointBuffer *buffer = ptbuf_new(n_items);
    if (buffer == NULL) {
        return NULL;
    }

    size_t n_read = fread(buffer->points, sizeof(Point), n_items, file);
    if (n_read != (size_t)n_items) {
        ptbuf_free(buffer);
        EMIT_ERR_MSG(IOError, "failed to read %d data points", n_items);
        return NULL;
    }
    return buffer;
}

static BNodeBuffer *bnodebuffer_new(int size) {
    BNodeBuffer *nodebuffer = (BNodeBuffer *)malloc(sizeof(BNodeBuffer));
    if (nodebuffer == NULL) {
        EMIT_ERR_MSG(MemoryError, "failed to allocate BNodeBuffer");
        return NULL;
    }

    nodebuffer->size = size;
    nodebuffer->next_free = 0;
    nodebuffer->nodes = (BallNode *)malloc(size * sizeof(BallNode));
    if (nodebuffer->nodes == NULL) {
        EMIT_ERR_MSG(MemoryError, "failed to allocate BNodeBuffer buffer");
        bnodebuffer_free(nodebuffer);
        return NULL;
    }
    return nodebuffer;
}

static void bnodebuffer_free(BNodeBuffer *buffer) {
    if (buffer->nodes != NULL) {
        free(buffer->nodes);
    }
}

static int bnodebuffer_get_next_free(BNodeBuffer *buffer) {
    return (buffer->next_free)++;
}

static int bnodebuffer_write(const BNodeBuffer *buffer, FILE *file) {
    size_t n_items = (size_t)buffer->size;
    size_t n_written = fwrite(buffer->nodes, sizeof(BallNode), n_items, file);
    if (n_written != n_items) {
        EMIT_ERR_MSG(IOError, "failed to write %zu nodes", n_items);
        return BTR_FAILED;
    }
    return BTR_SUCCESS;
}

static BNodeBuffer *bnodebuffer_read(int n_items, FILE *file) {
    BNodeBuffer *buffer = bnodebuffer_new(n_items);
    if (buffer == NULL) {
        return NULL;
    }

    size_t n_read = fread(buffer->nodes, sizeof(BallNode), n_items, file);
    if (n_read != (size_t)n_items) {
        EMIT_ERR_MSG(IOError, "failed to read %d nodes", n_items);
        bnodebuffer_free(buffer);
        return NULL;
    }
    return buffer;
}

static FileHeader *fileheader_new(const BallTree *tree) {
    FileHeader *header = (FileHeader *)malloc(sizeof(FileHeader));
    if (header == NULL) {
        EMIT_ERR_MSG(MemoryError, "failed to allocate FileHeader");
        return NULL;
    }

    header->leafsize = tree->leafsize;
    header->points = (SectionHeader){
        .n_items = tree->data.size,
        .itemsize = sizeof(*tree->data.points)
    };
    header->nodes = (SectionHeader){
        .n_items = balltree_count_nodes(tree),
        .itemsize = sizeof(*tree->root)
    };
    return header;
}

static int fileheader_write(const FileHeader *header, FILE *file) {
    const size_t n_items = 1;
    size_t n_written = fwrite(header, sizeof(*header), n_items, file);
    if (n_written != n_items) {
        EMIT_ERR_MSG(IOError, "failed to write file header");
        return BTR_FAILED;
    }
    return BTR_SUCCESS;
}

static FileHeader *fileheader_read(FILE *file) {
    FileHeader *header = (FileHeader *)malloc(sizeof(FileHeader));
    if (header == NULL) {
        EMIT_ERR_MSG(MemoryError, "failed to allocate FileHeader");
        return NULL;
    }

    const size_t n_items = 1;
    size_t n_read = fread(header, sizeof(FileHeader), n_items, file);
    if (n_read != n_items) {
        EMIT_ERR_MSG(IOError, "failed to read file header");
        return NULL;
    }
    return header;
}

static inline intptr_t child_ptr_substitute(const BallNode *child, BNodeBuffer *buffer) {
    return (child == NULL) ? 0 : bnodebuffer_get_next_free(buffer);
}

static int bnode_serialise(const BallNode *node, BNodeBuffer *buffer, int buf_idx) {
    if (buffer->next_free > buffer->size) {
        EMIT_ERR_MSG(IndexError, "buffer is too small to store further nodes");
        return BTR_FAILED;
    }

    // insert into buffer
    BallNode *stored = buffer->nodes + buf_idx;
    *stored = *node;
    stored->data.points = NULL;  // pointer not valid after deserialisation

    // replace pointers to childs by index in buffer they will be stored at
    intptr_t left_idx = child_ptr_substitute(node->left, buffer);
    intptr_t right_idx = child_ptr_substitute(node->right, buffer);
    stored->left = (BallNode *)left_idx;
    stored->right = (BallNode *)right_idx;

    // serialise childs recursively
    if (left_idx != 0) {
        if (bnode_serialise(node->left, buffer, left_idx) == BTR_FAILED) {
            return BTR_FAILED;
        }
    }
    if (right_idx != 0) {
        if (bnode_serialise(node->right, buffer, right_idx) == BTR_FAILED) {
            return BTR_FAILED;
        }
    }
    return BTR_SUCCESS;
}

int balltree_to_file(const BallTree *tree, const char *path) {
    FILE *file = fopen(path, "wb");
    if (file == NULL) {
        EMIT_ERR_MSG(OSError, "failed to open file: %s", path);
        return BTR_FAILED;
    }

    // create the file header and write it to the file
    FileHeader *header = fileheader_new(tree);
    if (header == NULL) {
        goto err_close_file;
    }
    if (fileheader_write(header, file) == BTR_FAILED) {
        goto err_dealloc_header;
    }

    // append the tree's data points to the file
    if (ptbuf_write(&tree->data, file) == BTR_FAILED) {
        goto err_dealloc_header;
    }

    // serialise the tree nodes and append them to the file
    BNodeBuffer *nodebuffer = bnodebuffer_new(header->nodes.n_items);
    if (nodebuffer == NULL) {
        goto err_dealloc_header;
    }
    free(header);
    int root_index = bnodebuffer_get_next_free(nodebuffer);
    if (bnode_serialise(tree->root, nodebuffer, root_index) == BTR_FAILED) {
        bnodebuffer_free(nodebuffer);
        goto err_close_file;
    }
    if (bnodebuffer_write(nodebuffer, file) == BTR_FAILED) {
        bnodebuffer_free(nodebuffer);
        goto err_close_file;
    }

    bnodebuffer_free(nodebuffer);
    if (fflush(file) == EOF) {
        EMIT_ERR_MSG(IOError, "failed to flush file");
        goto err_close_file;
    }
    fclose(file);
    return BTR_SUCCESS;

err_dealloc_header:
    free(header);
err_close_file:
    fclose(file);
    return BTR_FAILED;
}

static BallNode *bnode_deserialise(
    const BNodeBuffer *buffer,
    Point *points,
    int buf_idx
) {
    if (buf_idx >= buffer->size) {
        EMIT_ERR_MSG(IndexError, "node index exceeds node buffer size");
        return NULL;
    }
    BallNode *stored = buffer->nodes + buf_idx;

    // create a new node instance
    BallNode *node = (BallNode *)malloc(sizeof(BallNode));
    if (node == NULL) {
        EMIT_ERR_MSG(MemoryError, "failed to allocate BallNode");
        return NULL;
    }
    *node = *stored;

    // restore the pointers with valid addresses
    node->data.points = points;
    int left_idx = (intptr_t)node->left;
    int right_idx = (intptr_t)node->right;
    if (left_idx != 0) {
        node->left = bnode_deserialise(buffer, points, left_idx);
        if (node->left == NULL) {
            free(node);
            return NULL;
        }
    }
    if (right_idx != 0) {
        node->right = bnode_deserialise(buffer, points, right_idx);
        if (node->right == NULL) {
            bnode_free(node->left);
            free(node);
            return NULL;
        }
    }
    return node;
}

BallTree* balltree_from_file(const char *path) {
    BallTree *tree = (BallTree *)malloc(sizeof(BallTree));
    if (tree == NULL) {
        EMIT_ERR_MSG(MemoryError, "failed to allocate BallTree");
        return NULL;
    }

    FILE *file = fopen(path, "rb");
    if (file == NULL) {
        EMIT_ERR_MSG(OSError, "failed to open file: %s", path);
        goto err_dealloc_tree;
    }

    // read header from the file
    FileHeader *header = fileheader_read(file);
    if (header == NULL) {
        goto err_close_file;
    }
    tree->leafsize = header->leafsize;

    // read tree's data points from the file
    PointBuffer *points = ptbuf_read(header->points.n_items, file);
    if (points == NULL) {
        goto err_dealloc_header;
    }
    tree->data = *points;  // move ownership of underlying point buffer
    free(points);

    // read the tree nodes from the file and deserialise them
    BNodeBuffer *nodebuffer = bnodebuffer_read(header->nodes.n_items, file);
    if (nodebuffer == NULL) {
        goto err_dealloc_header;  // points deallocated at err_dealloc_tree
    }
    tree->root = bnode_deserialise(nodebuffer, tree->data.points, 0);
    bnodebuffer_free(nodebuffer);
    if (tree->root == NULL) {
        goto err_dealloc_header;  // nodes deallocated at err_dealloc_tree
    }

    free(header);
    fclose(file);
    return tree;

err_dealloc_header:
    free(header);
err_close_file:
    fclose(file);
err_dealloc_tree:
    balltree_free(tree);
    return NULL;
}
