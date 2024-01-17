CC = gcc
CXX = g++
CFLAGS = -Iinclude -Wall
CXXFLAGS = -Iinclude -Wall
OBJS = obj/geoprog.o obj/geograd.o obj/kspace.o obj/pppm.o obj/array.o obj/util.o obj/cgmin.o
EXEC = bin/geo bin/cgmin-example

all: $(EXEC) $(OBJS)

bin/geo: obj/geoprog.o obj/geograd.o obj/fns.o
	$(CC) -o $@ $^
bin/cgmin-example: src/cgmin.c
	$(CC) $(CFLAGS) -DSIMPLE_EXAMPLE -o $@ $^
obj/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@
obj/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
clean:
	rm -f obj/*.o $(EXEC)

.PHONY: all clean
