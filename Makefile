CC = gcc
CFLAGS = -Iinclude -Wall -O3 -ffast-math
LDFLAGS = -shared
LDLIBS = -lz
SRCDIR = src
BUILDDIR = build
BINDIR = bin

SRCS = $(wildcard $(SRCDIR)/*.c)
OBJS = $(patsubst $(SRCDIR)/%.c, $(BUILDDIR)/%.o, $(SRCS))

LIBDIR = lib
LIBS = $(LIBDIR)/balltree.so

TARGET = $(BINDIR)/main.out

all: $(TARGET) $(LIBS)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(LDLIBS) $^ -o $@

$(BUILDDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(LIBS): $(OBJS)
	$(CC) $(LDFLAGS) $(LDLIBS) $^ -o $@

clean:
	rm -rf $(BUILDDIR)/*.o $(LIBDIR)/* $(BINDIR)/*

.PHONY: clean
