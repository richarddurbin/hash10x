# makefile for hash10x, developed on Richard's Mac.

CFLAGS= -g -O3

all: fq2b hash10x

clean:
	$(RM) *.o *~ fq2b hash10x

### object files

UTILS_OBJS=hash.o dict.o array.o utils.o
UTILS_HEADERS=utils.h array.h dict.h hash.h
$(UTILS_OBJS): utils.h $(UTILS_HEADERS)

readseq.o: readseq.h

minhash.o: minhash.h

### programs

fq2b: fq2b.c $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ -lz

hash10x: hash10x.c readseq.o minhash.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ -lm

### end of file
