# Compiler settings
CC = gcc
CFLAGS = -Wall -O3

# Paths relative to the tests/ directory
INCLUDES = -I../src/cpyrcolate
SRCS = ../src/cpyrcolate/percolate_core.c test_core.c
TARGET = ../build/test_percolation

# Default target when you just type 'make'
.PHONY: all build run clean

all: run

# Compile the executable
build:
	mkdir -p ../build
	$(CC) $(CFLAGS) $(INCLUDES) $(SRCS) -o $(TARGET)

# Run the executable (depends on 'build' so it always compiles first if needed)
run: build
	./$(TARGET)

# Clean up the compiled binary
clean:
	rm -f $(TARGET)
