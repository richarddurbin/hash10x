# makefile for hash10x, developed on Richard's Mac.

#CFLAGS= -O3
CFLAGS= -g				# for debugging
#CFLAGS= -03 -DOMP -fopenmp		# for OMP parallelisation - doesn't compile on Mac

all: fq2b hash10x

clean:
	$(RM) *.o *~ fq2b hash10x moshmap moshasm

### object files

UTILS_OBJS=hash.o dict.o array.o utils.o
UTILS_HEADERS=utils.h array.h dict.h hash.h
$(UTILS_OBJS): utils.h $(UTILS_HEADERS)

readseq.o: readseq.h

seqhash.o: seqhash.h

moshset.c: moshset.h

### programs

fq2b: fq2b.c $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ -lz

hash10x: hash10x.c readseq.o seqhash.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ -lm

moshmap: moshmap.c readseq.o seqhash.o moshset.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ -lm

moshasm: moshasm.c readseq.o seqhash.o moshset.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ -lm

### end of file
