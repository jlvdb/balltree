CC = gcc
CFLAGS = -Iinclude -Wall -O3 -ffast-math -fPIC
LDFLAGS = -shared
SRCDIR = src
BUILDDIR = build
BINDIR = bin

SRCS = $(wildcard $(SRCDIR)/*.c)
OBJS = $(patsubst $(SRCDIR)/%.c, $(BUILDDIR)/%.o, $(SRCS))

LIBDIR = lib
LIBS = $(LIBDIR)/libballtree.so

TARGET = $(BINDIR)/main.out

all: $(TARGET) $(LIBS)

$(TARGET): $(OBJS)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $^ -o $@

$(BUILDDIR)/%.o: $(SRCDIR)/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

$(LIBS): $(OBJS)
	@mkdir -p $(@D)
	$(CC) $(LDFLAGS) $^ -o $@

clean:
	rm -rf $(BUILDDIR)/*.o $(BINDIR)/* $(LIBDIR)/*.so

.PHONY: clean
