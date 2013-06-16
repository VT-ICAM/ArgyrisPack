NUMERIC_PATH=./ap/numeric
CFLAGS=-DLAPACKINDEX=int -D$(STORAGE_ORDER) -Wall -Wextra -pedantic
CFLAGS+=-std=c99 -fPIC -O3 -I$(NUMERIC_PATH)

all : $(NUMERIC_PATH)/*.c $(NUMERIC_PATH)/*.h
	$(CC) $(CFLAGS) -c $(NUMERIC_PATH)/argyris_pack.c

so : all
	$(CC) -shared -o libargyris_pack.so argyris_pack.o -lblas -lm
