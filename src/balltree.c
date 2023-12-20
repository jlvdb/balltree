#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "point.h"
#include "ballnode.h"
#include "balltree.h"

#define SUCCESS 1
#define FAILED  0

#define DEFAULT_LEAFSIZE 20

struct CompressionHeader {
    int size;
    int bytes;
    int compressed;
};


struct FileHeader {
    int leafsize;
    struct CompressionHeader nodes;
    struct CompressionHeader points;
};


void fileheader_print(struct FileHeader *header)
{
    printf("Header {\n");
    printf("    leafsize=%d,\n", header->leafsize);
    printf("    nodes={size=%d, bytes=%d, compressed=%d},\n", header->nodes.size, header->nodes.bytes, header->nodes.compressed);
    printf("    points={size=%d, bytes=%d, compressed=%d}\n", header->points.size, header->points.bytes, header->points.compressed);
    printf("}\n");
}


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

int compress_and_append(FILE *file, size_t n_bytes, void *buffer)
{
    // compress the buffer nodes
    uLongf compressed_size = compressBound(n_bytes);
    Bytef *compressed_buffer = (Bytef *)malloc(compressed_size);
    if (!compressed_buffer) {
        fprintf(stderr, "ERROR: failed to allocate memory for compressed node data\n");
        return -1;
    }
    if (compress(compressed_buffer, &compressed_size, (Bytef *)buffer, n_bytes) != Z_OK) {
        fprintf(stderr, "ERROR: failed to compress node data\n");
        free(compressed_buffer);
        return -1;
    }

    // write compressed buffer
    size_t elements_written = fwrite(compressed_buffer, sizeof(Bytef), compressed_size, file);
    free(compressed_buffer);
    if (elements_written != compressed_size) {
        fprintf(stderr, "ERROR: failed to write compressed data\n");
        return -1;
    }
    return compressed_size;
}

int balltree_to_file(const struct BallTree *tree, const char *path)
{
    int success;
    FILE *file = fopen(path, "wb");
    if (!file) {
        fprintf(stderr, "ERROR: failed to open file\n");
        return FAILED;
    }

    // initialise the header
    int num_nodes = balltree_count_nodes(tree);
    int num_points = tree->data.size;
    struct FileHeader header = {
        .leafsize = tree->leafsize,
        .nodes = (struct CompressionHeader){
            .size = num_nodes,
            .bytes = sizeof(struct BallNodeSerialized) * num_nodes,
            .compressed = -1,
        },
        .points = (struct CompressionHeader){
            .size = num_points,
            .bytes = sizeof(struct Point) * num_points,
            .compressed = -1,
        },
    };
    if (fwrite(&header, sizeof(struct FileHeader), 1, file) != 1) {
        fprintf(stderr, "ERROR: failed to write header\n");
        goto close_file;
    }

    // serialise nodes
    int index_tracker = 1;  // 0 is already reserved for root node
    struct BallNodeBuffer node_buffer = {
        .size = header.nodes.size,
        .next_free = &index_tracker,
    };
    node_buffer.buffer = (struct BallNodeSerialized*)malloc(header.nodes.bytes);
    if (!node_buffer.buffer) {
        fprintf(stderr, "ERROR: failed to allocate memory for serialized node data\n");
        goto close_file;
    }
    success = ballnode_serialise_recursive(node_buffer, tree->root, 0);
    if (!success) {
        goto dealloc_nodes;
    }

    // write compressed data
    header.points.compressed = compress_and_append(file, header.points.bytes, tree->data.points);
    if (header.points.compressed == -1) {
        goto dealloc_nodes;
    }
    header.nodes.compressed = compress_and_append(file, header.nodes.bytes, node_buffer.buffer);
    free(node_buffer.buffer);
    if (header.nodes.compressed == -1) {
        goto close_file;
    }

    // rewrite the header
    if (fseek(file, 0, SEEK_SET) != 0) {
        fprintf(stderr, "ERROR: failed to the beginning of the file\n");
        goto close_file;
    }
    if (fwrite(&header, sizeof(struct FileHeader), 1, file) != 1) {
        fprintf(stderr, "ERROR: failed to update header\n");
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

void* decompress_and_read(FILE *file, struct CompressionHeader header)
{
    Bytef *compressed_buffer = (Bytef *)malloc(header.compressed);
    if (!compressed_buffer) {
        fprintf(stderr, "ERROR: failed to allocate memory for compressed data\n");
        return NULL;
    }

    // read the compressed data
    size_t elements_read = fread(compressed_buffer, sizeof(Bytef), header.compressed, file);
    if (elements_read != (size_t)header.compressed) {
        fprintf(stderr, "ERROR: failed to read compressed data from file\n");
        goto dealloc_compressed;
    }

    // setup decompression
    z_stream stream;
    memset(&stream, 0, sizeof(stream));
    if (inflateInit(&stream) != Z_OK) {
        fprintf(stderr, "ERROR: failed to initialize zlib stream for decompression\n");
        goto dealloc_compressed;
    }
    stream.next_in = compressed_buffer;
    stream.avail_in = header.compressed;
    void *buffer = malloc(header.bytes);
    if (!buffer) {
        fprintf(stderr, "ERROR: failed to allocate memory for decompressed data\n");
        goto dealloc_inflate;
    }

    // decompress data
    stream.next_out = (Bytef *)buffer;
    stream.avail_out = header.bytes;
    if (inflate(&stream, Z_FINISH) != Z_STREAM_END) {
        fprintf(stderr, "ERROR: failed to decompress data\n");
        goto dealloc_inflate;
    }
    inflateEnd(&stream);

    if (stream.total_out != (size_t)header.bytes) {
        fprintf(stderr, "ERROR: decompressed size does not match expected size\n");
        goto dealloc_raw;
    }

    // Free the compressed buffer
    free(compressed_buffer);
    return buffer;

dealloc_inflate:
    inflateEnd(&stream);
dealloc_raw:
    free(buffer);
dealloc_compressed:
    free(compressed_buffer);
    return NULL;
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

    // read and decompress the point data
    struct Point *points = decompress_and_read(file, header.points);
    if (!points) {
        goto close_file;
    }
    tree->data = (struct PointBuffer){
        .size = header.points.size,
        .points = points,
    };

    // read and decompress node data
    struct BallNodeSerialized *node_buffer = decompress_and_read(file, header.nodes);
    if (!node_buffer) {
        goto dealloc_points;
    }
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
