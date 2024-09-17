CC = gcc
CXX = g++
CFLAGS = -Iinclude -Wall
CXXFLAGS = -Iinclude -Wall
OBJS = obj/geoprog.o obj/geograd.o obj/kspace.o obj/pppm.o obj/array.o obj/util.o obj/cgmin.o
EXEC = bin/geo bin/cgmin-example bin/pppm-example

all: test

test: $(EXEC) $(OBJS)

bin/geo: obj/geoprog.o obj/geograd.o obj/fns.o
	$(CXX) -o $@ $^ -lm
bin/cgmin-example: src/cgmin.c
	$(CC) $(CFLAGS) -DSIMPLE_EXAMPLE -o $@ $^ -lm
bin/pppm-example: src/pppm.c obj/array.o obj/util.o
	$(CC) $(CFLAGS) -DSIMPLE_EXAMPLE -o $@ $^ -lgsl -lfftw3 -lm
obj/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@ -lm
obj/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ -lm
clean:
	rm -f obj/*.o $(EXEC)
test:
	bin/cgmin-example
	bin/pppm-example

.PHONY: all clean
