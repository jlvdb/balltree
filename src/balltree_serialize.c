#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "point.h"
#include "balltree.h"
#include "ballnode.h"
#include "balltree_macros.h"

static int ptbuf_write(const PointBuffer *buffer, FILE *file);
static PointBuffer *ptbuf_read(int64_t n_items, FILE *file);

typedef struct {
    BallNode *nodes;
    size_t next_free;
    size_t size;
} BNodeBuffer;

static BNodeBuffer *bnodebuffer_new(int64_t size);
static void bnodebuffer_free(BNodeBuffer *buffer);
static size_t bnodebuffer_get_next_free(BNodeBuffer *buffer);
static int bnodebuffer_write(const BNodeBuffer *buffer, FILE *file);
static BNodeBuffer *bnodebuffer_read(int64_t n_items, FILE *file);

typedef struct {
    int64_t n_items;
    int itemsize;
} SectionHeader;

typedef struct {
    SectionHeader points;
    SectionHeader nodes;
    int leafsize;
    int data_owned;
} FileHeader;

static FileHeader *fileheader_new(const BallTree *tree);
static int fileheader_write(const FileHeader *header, FILE *file);
static FileHeader *fileheader_read(FILE *file);

static int bnode_serialise(const BallNode *node, BNodeBuffer *buffer, size_t buf_idx, PointBuffer *points);
static BallNode *bnode_deserialise(const BNodeBuffer *buffer, size_t buf_idx, PointBuffer *points);


static int ptbuf_write(const PointBuffer *buffer, FILE *file) {
    size_t n_items = (size_t)buffer->size;
    size_t n_written = fwrite(buffer->points, sizeof(Point), n_items, file);
    if (n_written != n_items) {
        EMIT_ERR_MSG(IOError, "failed to write %zu data points", n_items);
        return BTR_FAILED;
    }
    return BTR_SUCCESS;
}

static PointBuffer *ptbuf_read(int64_t n_items, FILE *file) {
    PointBuffer *buffer = ptbuf_new(n_items);
    if (buffer == NULL) {
        return NULL;
    }

    size_t n_read = fread(buffer->points, sizeof(Point), n_items, file);
    if (n_read != (size_t)n_items) {
        ptbuf_free(buffer);
        EMIT_ERR_MSG(IOError, "failed to read %lld data points", n_items);
        return NULL;
    }
    return buffer;
}

static BNodeBuffer *bnodebuffer_new(int64_t size) {
    BNodeBuffer *nodebuffer = malloc(sizeof(BNodeBuffer));
    if (nodebuffer == NULL) {
        EMIT_ERR_MSG(MemoryError, "failed to allocate BNodeBuffer");
        return NULL;
    }

    nodebuffer->nodes = malloc(size * sizeof(BallNode));
    if (nodebuffer->nodes == NULL) {
        EMIT_ERR_MSG(MemoryError, "failed to allocate BNodeBuffer buffer");
        bnodebuffer_free(nodebuffer);
        return NULL;
    }
    nodebuffer->size = size;
    nodebuffer->next_free = 0L;
    return nodebuffer;
}

static void bnodebuffer_free(BNodeBuffer *buffer) {
    if (buffer->nodes != NULL) {
        free(buffer->nodes);
    }
    free(buffer);
}

static size_t bnodebuffer_get_next_free(BNodeBuffer *buffer) {
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

static BNodeBuffer *bnodebuffer_read(int64_t n_items, FILE *file) {
    BNodeBuffer *buffer = bnodebuffer_new(n_items);
    if (buffer == NULL) {
        return NULL;
    }

    size_t n_read = fread(buffer->nodes, sizeof(BallNode), n_items, file);
    if (n_read != (size_t)n_items) {
        EMIT_ERR_MSG(IOError, "failed to read %lld nodes", n_items);
        bnodebuffer_free(buffer);
        return NULL;
    }
    return buffer;
}

static FileHeader *fileheader_new(const BallTree *tree) {
    FileHeader *header = malloc(sizeof(FileHeader));
    if (header == NULL) {
        EMIT_ERR_MSG(MemoryError, "failed to allocate FileHeader");
        return NULL;
    }

    header->points = (SectionHeader){
        .n_items = tree->data->size,
        .itemsize = sizeof(Point),
    };
    header->nodes = (SectionHeader){
        .n_items = balltree_count_nodes(tree),
        .itemsize = sizeof(BallNode),
    };
    header->leafsize = tree->leafsize;
    header->data_owned = tree->data_owned;
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
    FileHeader *header = malloc(sizeof(FileHeader));
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

static int bnode_serialise(
    const BallNode *node,
    BNodeBuffer *buffer,
    size_t buf_idx,
    PointBuffer *points
) {
    if (buffer->next_free > buffer->size) {
        EMIT_ERR_MSG(IndexError, "buffer is too small to store further nodes");
        return BTR_FAILED;
    }

    // insert copy of current node into buffer
    BallNode *stored = buffer->nodes + buf_idx;
    *stored = *node;

    if (BALLNODE_IS_LEAF(node)) {
        // replace pointers to slice start/end by index into point data buffer,
        // see bnode_deserialise()
        stored->data.start = (Point *)(node->data.start - points->points);
        stored->data.end = (Point *)(node->data.end - points->points);
    } else {
        // replace pointers to childs by index in buffer they will be stored at
        size_t left_idx = bnodebuffer_get_next_free(buffer);
        size_t right_idx = bnodebuffer_get_next_free(buffer);
        stored->childs.left = (BallNode *)left_idx;
        stored->childs.right = (BallNode *)right_idx;

        // serialise childs recursively
        if (bnode_serialise(node->childs.left, buffer, left_idx, points) == BTR_FAILED) {
            return BTR_FAILED;
        }
        if (bnode_serialise(node->childs.right, buffer, right_idx, points) == BTR_FAILED) {
            return BTR_FAILED;
        }
    }
    return BTR_SUCCESS;
}

int balltree_to_file(const BallTree *tree, const char *path) {
    int return_value = BTR_FAILED;

    FILE *file = fopen(path, "wb");
    if (file == NULL) {
        EMIT_ERR_MSG(OSError, "failed to open file: %s", path);
        goto error;
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
    if (ptbuf_write(tree->data, file) == BTR_FAILED) {
        goto err_dealloc_header;
    }

    // serialise the tree nodes and append them to the file
    BNodeBuffer *nodebuffer = bnodebuffer_new(header->nodes.n_items);
    if (nodebuffer == NULL) {
        goto err_dealloc_header;
    }
    size_t root_index = bnodebuffer_get_next_free(nodebuffer);
    if (bnode_serialise(tree->root, nodebuffer, root_index, tree->data) == BTR_FAILED) {
        goto err_dealloc_bnodebuffer;
    }
    if (bnodebuffer_write(nodebuffer, file) == BTR_FAILED) {
        goto err_dealloc_bnodebuffer;
    }

    if (fflush(file) == EOF) {
        EMIT_ERR_MSG(IOError, "failed to flush file");
        goto err_dealloc_bnodebuffer;
    }

    return_value = BTR_SUCCESS;

err_dealloc_bnodebuffer:
    bnodebuffer_free(nodebuffer);
err_dealloc_header:
    free(header);
err_close_file:
    fclose(file);
error:
    return return_value;
}

static BallNode *bnode_deserialise(
    const BNodeBuffer *buffer,
    size_t buf_idx,
    PointBuffer *points
) {
    if (buf_idx >= buffer->size) {
        EMIT_ERR_MSG(IndexError, "node index exceeds node buffer size");
        return NULL;
    }
    BallNode *stored = buffer->nodes + buf_idx;

    // create a new node instance
    BallNode *node = malloc(sizeof(BallNode));
    if (node == NULL) {
        EMIT_ERR_MSG(MemoryError, "failed to allocate BallNode");
        return NULL;
    }
    *node = *stored;

    if (BALLNODE_IS_LEAF(node)) {
        // restore pointers to slice start/end from their index into the point
        // data buffer, see bnode_serialise()
        node->data.start = points->points + (size_t)node->data.start;
        node->data.end = points->points + (size_t)node->data.end;
    } else {
        // restore child instances from their index into node buffer
        size_t left_idx = (size_t)node->childs.left;
        size_t right_idx = (size_t)node->childs.right;
        node->childs.left = bnode_deserialise(buffer, left_idx, points);
        if (node->childs.left == NULL) {
            free(node);
            return NULL;
        }
        node->childs.right = bnode_deserialise(buffer, right_idx, points);
        if (node->childs.right == NULL) {
            free(node);
            return NULL;
        }
    }

    return node;
}

BallTree* balltree_from_file(const char *path) {
    BallTree *tree = NULL;  // return value

    FILE *file = fopen(path, "rb");
    if (file == NULL) {
        EMIT_ERR_MSG(OSError, "failed to open file: %s", path);
        goto error;
    }

    // read header from the file
    FileHeader *header = fileheader_read(file);
    if (header == NULL) {
        goto err_close_file;
    }

    // read tree's data points from the file
    PointBuffer *points = ptbuf_read(header->points.n_items, file);
    if (points == NULL) {
        goto err_dealloc_header;
    }

    // read the tree nodes from the file and deserialise them
    BNodeBuffer *nodebuffer = bnodebuffer_read(header->nodes.n_items, file);
    if (nodebuffer == NULL) {
        goto err_dealloc_points;
    }
    BallNode *root = bnode_deserialise(nodebuffer, 0L, points);
    if (root == NULL) {
        goto err_dealloc_bnodebuffer;
    }
    
    tree = malloc(sizeof(BallTree));
    if (tree == NULL) {
        EMIT_ERR_MSG(MemoryError, "failed to allocate BallTree");
        goto err_dealloc_root;
    }
    tree->leafsize = header->leafsize;
    tree->data_owned = 1;  // after deserialising, the buffer from the file is used
    tree->data = points;  // move ownership of underlying point buffer
    tree->root = root;

err_dealloc_root:
    if (tree == NULL) {
        bnode_free(root);
    }
err_dealloc_bnodebuffer:
    bnodebuffer_free(nodebuffer);
err_dealloc_points:
    if (tree == NULL) {
        free(points);
    }
err_dealloc_header:
    free(header);
err_close_file:
    fclose(file);
error:
    return tree;
}
