CC = gcc
CXX = g++
CFLAGS = -Iinclude -Wall
CXXFLAGS = -Iinclude -Wall
OBJS = obj/geoprog.o obj/geograd.o
EXEC = bin/geo

all: $(EXEC)

bin/geo: obj/geoprog.o obj/geograd.o obj/fns.o
	$(CC) -o $@ $^
obj/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@
obj/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
clean:
	rm -f obj/*.o $(EXEC)

.PHONY: all clean
