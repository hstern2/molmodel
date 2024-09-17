CC = gcc
CXX = g++
CFLAGS = -Iinclude -Wall
CXXFLAGS = -Iinclude -Wall
OBJS = obj/geoprog.o obj/geograd.o obj/kspace.o obj/pppm.o obj/array.o obj/util.o obj/cgmin.o
EXEC = bin/geo bin/cgmin-example bin/pppm-example

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

bin/pppm-example: src/pppm.c obj/array.o obj/util.o
	$(CC) $(CFLAGS) -DSIMPLE_EXAMPLE -o $@ $^ -lgsl -lfftw3 $(BLAS_LIB) -lm

obj/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

obj/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f obj/*.o $(EXEC)

test:
	bin/cgmin-example
	bin/pppm-example

.PHONY: all clean
