CFLAGS=-DLAPACKINDEX=int -D$(STORAGE_ORDER) -Wall -std=c99 -fPIC -O3 -I./ap/numerical/
NUMERICAL_PATH=./ap/numerical/
all : $(NUMERICAL_PATH)*.c $(NUMERICAL_PATH)*.h
	$(CC) $(CFLAGS) -c ./ap/numerical/argyris_pack.c

so : all
	$(CC) -shared -o libargyris_pack.so argyris_pack.o -lblas -lm

# as a quick test, link the test.c file to the shared object built in the
# current directory.
test : so
	$(CC) -L./ -std=c99 -Wall -D$(STORAGE_ORDER) test.c -o test -largyris_pack
