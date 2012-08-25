CFLAGS=-DLAPACKINDEX=int -DUSE_ROW_MAJOR -std=c99
CC=gcc

all : *.c *.h
	$(CC) $(CFLAGS) -c argyris_pack.c -lblas