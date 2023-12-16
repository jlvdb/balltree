#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_BUFFER_SIZE 4096 // Adjust this based on your requirements

struct Buffer {
    char *buffer;
    size_t total_size;
    size_t used_size;
    FILE *file_ptr;
};

void init_buffer(struct Buffer *buf, size_t size) {
    buf->buffer = (char *)malloc(size);
    if (buf->buffer == NULL) {
        perror("Failed to allocate memory for buffer");
        exit(EXIT_FAILURE);
    }
    buf->total_size = size;
    buf->used_size = 0;
    buf->file_ptr = NULL;
}

void add_bytes(struct Buffer *buf, const char *bytes, size_t size) {
    // Check if there is enough space in the buffer and flush to disk if needed
    if (buf->used_size + size > buf->total_size) {
        if (buf->file_ptr != NULL) {
            fwrite(buf->buffer, 1, buf->used_size, buf->file_ptr);
            fflush(buf->file_ptr);
        }
        buf->used_size = 0;
    }

    memcpy(buf->buffer + buf->used_size, bytes, size);
    buf->used_size += size;
}

void finalize_buffer(struct Buffer *buf) {
    // Flush any remaining bytes to disk
    if (buf->file_ptr != NULL) {
        fwrite(buf->buffer, 1, buf->used_size, buf->file_ptr);
        fflush(buf->file_ptr);
    }

    if (buf->file_ptr != NULL) {
        fclose(buf->file_ptr);
    }
    free(buf->buffer);
}
