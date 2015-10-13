CC=gcc
CFLAGS=-Wall -pedantic -std=c99
LFLAGS=-lm

OBJS = argon.o main.o

all: program

program: $(OBJS)
		$(CC) $(LFLAGS) $(CFLAGS) $^ -o $@

clean:
		rm -f *.o main

.PHONY: all clean
