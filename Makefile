CC = gcc
CXX = g++
CFLAGS = -Iinclude -Wall
CXXFLAGS = -Iinclude -Wall
SRC = $(wildcard src/*.c src/*.cpp)
OBJS = $(patsubst src/%.c, obj/%.o, $(SRC)) $(patsubst src/%.cpp, obj/%.o, $(SRC))
EXEC = bin/geo bin/cgmin-example bin/pppm-example bin/rmsd bin/msim

# Platform-specific libraries
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S), Darwin) # macOS
    BLAS_LIB = -L/opt/homebrew/opt/openblas/lib -lopenblas
else # Assuming Linux (Ubuntu)
    BLAS_LIB = -lopenblas
endif

all: test

test: $(EXEC) $(OBJS)

bin/msim: \
obj/atom.o obj/atomio.o obj/msys.o obj/msysio.o obj/msimprog.o obj/msimprogio.o \
obj/nbrlist.o obj/mdyn.o obj/mdynio.o obj/elecljio.o obj/checkgr.o obj/checkgrio.o \
obj/rspace.o obj/kspace.o obj/callback.o obj/callbackio.o obj/rdf.o obj/rdfio.o obj/intra.o \
obj/intraio.o obj/mmin.o obj/mminio.o obj/cgmin.o obj/eleclj.o obj/pppm.o obj/pdb.o \
obj/geohist.o obj/geohistio.o obj/apattern.o obj/acid.o obj/acidio.o obj/acidana.o \
obj/constraint.o obj/constraintio.o obj/covalent.o obj/covalentio.o obj/contact.o \
obj/contactio.o obj/rmsdcb.o obj/rmsdcbio.o obj/waterintra.o obj/waterintraio.o \
obj/cmmsd.o obj/cmmsdio.o obj/harmonic.o obj/harmonicio.o obj/zmatrix.o obj/zmatrixio.o \
obj/espfit.o obj/espfitio.o obj/bcifit.o obj/bcifitio.o obj/avgmf.o obj/avgmfio.o \
obj/cmpoleio.o obj/ppot.o obj/ppotio.o obj/fns.o obj/util.o obj/array.o obj/progio.o \
obj/prepro.o obj/histio.o obj/prog.o obj/random.o obj/cheby.o obj/pbin.o obj/out.o \
obj/bcond.o obj/linalg.o obj/geograd.o obj/rmsd.o obj/euler.o obj/hist.o obj/bcondtype.o
	$(CXX) -o $@ $^ -lgsl -lfftw3 $(BLAS_LIB) -lm

bin/geo: obj/geoprog.o obj/geograd.o obj/fns.o
	$(CXX) -o $@ $^ $(BLAS_LIB) -lm

bin/cgmin-example: src/cgmin.c
	$(CC) $(CFLAGS) -DSIMPLE_EXAMPLE -o $@ $^ $(BLAS_LIB) -lm

bin/pppm-example: obj/pppm.o obj/array.o obj/util.o obj/pppm-example.o
	$(CC) $(CFLAGS) -DSIMPLE_EXAMPLE -o $@ $^ -lgsl -lfftw3 $(BLAS_LIB) -lm

bin/rmsd: obj/rmsdprog.o obj/rmsd.o obj/out.o obj/linalg.o obj/fns.o obj/util.o obj/array.o
	$(CXX) -o $@ $^ $(BLAS_LIB) -lm

obj/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

obj/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

src/%io.cpp: include/%.hpp
	./io -o $< > $@


clean:
	rm -f obj/*.o $(EXEC) src/*io.cpp

test:
	bin/cgmin-example
	@echo
	bin/pppm-example
	@echo
	bin/rmsd dat/v1 dat/v2

.PHONY: all clean
