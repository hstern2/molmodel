CC = gcc
CXX = g++
CFLAGS = -Iinclude -Wall
CXXFLAGS = -Iinclude -Wall
SRC = $(wildcard src/*.c src/*.cpp)
OBJS = $(patsubst src/%.c, obj/%.o, $(SRC)) $(patsubst src/%.cpp, obj/%.o, $(SRC))
EXEC = bin/geo bin/cgmin-example bin/pppm-example bin/rmsd

# Platform-specific libraries
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S), Darwin) # macOS
    BLAS_LIB = -L/opt/homebrew/opt/openblas/lib -lopenblas
else # Assuming Linux (Ubuntu)
    BLAS_LIB = -lopenblas
endif

all: test

test: $(EXEC) $(OBJS)

bin/geo: obj/geoprog.o obj/geograd.o obj/fns.o
	$(CXX) -o $@ $^ $(BLAS_LIB) -lm

bin/cgmin-example: src/cgmin.c
	$(CC) $(CFLAGS) -DSIMPLE_EXAMPLE -o $@ $^ $(BLAS_LIB) -lm

bin/pppm-example: obj/pppm.o obj/array.o obj/util.o obj/pppm-example.o
	$(CC) $(CFLAGS) -DSIMPLE_EXAMPLE -o $@ $^ -lgsl -lfftw3 $(BLAS_LIB) -lm

bin/rmsd: obj/rmsdprog.o obj/rmsd.o obj/out.o obj/lapack.o obj/fns.o
	$(CXX) -o $@ $^ $(BLAS_LIB) -lm

obj/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

obj/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f obj/*.o $(EXEC)

test:
	bin/cgmin-example
	@echo
	bin/pppm-example
	@echo
	bin/rmsd dat/v1 dat/v2

.PHONY: all clean
