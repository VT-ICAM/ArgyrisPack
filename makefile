CFLAGS=-DLAPACKINDEX=int -DUSE_ROW_MAJOR -Wall -std=c99 -fPIC
CC=/usr/bin/cc

all : *.c *.h
	$(CC) $(CFLAGS) -O3 -c argyris_pack.c

so : all
	$(CC) -shared -o libargyris_pack.so argyris_pack.o -lblas -lm

test : so
	$(CC) -L./ -std=c99 -Wall test.c -o test -largyris_pack