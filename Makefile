CXX=gcc
CXXFLAGS=-Wall
LFLAGS=-lm

OBJS = argon.o main.o

all: program

program: $(OBJS)
		$(CXX) $(CXXFLAGS) $^ -o $@

clean:
		rm -f *.o main

.PHONY: all clean
