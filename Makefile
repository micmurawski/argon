CC=gcc
CFLAGS=-O2 -ffast-math -Wall -pedantic -std=c99
LFLAGS=-lm

OBJS = argon.o main.o

all: program

program: $(OBJS)
		$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

clean:
		rm -f *.o main

.PHONY: all clean
