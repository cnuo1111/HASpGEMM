CC = gcc
CFLAGS = -O3 -fopenmp -g -I include

all: hash spa esc

hash: ./src/hash.c
	$(CC) $(CFLAGS) ./src/hash.c -o hash

spa: ./src/spa.c
	$(CC) $(CFLAGS) ./src/spa.c -o spa

esc: ./src/esc.c
	$(CC) $(CFLAGS) ./src/esc.c -o esc

clean:
	rm -f hash spa esc
