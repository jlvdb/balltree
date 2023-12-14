CC = gcc
CFLAGS = -Iinclude -Wall -O
SRCDIR = src
BUILDDIR = build
BINDIR = bin

SRCS = $(wildcard $(SRCDIR)/*.c)
OBJS = $(patsubst $(SRCDIR)/%.c, $(BUILDDIR)/%.o, $(SRCS))

TARGET = $(BINDIR)/main.out

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@

$(BUILDDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILDDIR)/*.o $(BINDIR)/*

.PHONY: clean
