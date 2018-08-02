/*  File: hash10x.c
 *  Author: Richard Durbin (rd109@gen.cam.ac.uk)
 *  Copyright (C) Richard Durbin, 2018
 *-------------------------------------------------------------------
 * Description: minHash-based analysis of 10X Genomics linked read data sets
    reads fqb binary output, and writes as custom binary, by convention .hash 
    look to identify read clusters, and het kmers
    can evaluate with crib built from external fasta files
 * Exported functions:
 * HISTORY:
 * Last edited: Aug  2 09:50 2018 (rd)
 * Created: Mon Mar  5 11:38:47 2018 (rd)
 *-------------------------------------------------------------------
 */

#include "minhash.h"
#include "array.h"
#include <ctype.h>
#ifdef OMP
#include <omp.h>
#endif

struct {
  int k ;
  int w ;
  int r ;
  int B ;
  int N ;
  int chunkSize ;
  int clusterThreshold ;
} params ;

typedef struct {
  U32 hash ;			/* actually the hash index, not the original hash */
  U16 read ;			/* 0..(nRead-1) */
  U8 cluster ;		        /* 0 means not in a cluster, 1..nCluster assigns to a cluster */
  U8 isHet:1 ;			/* is heterozygous - single copy */
  U8 isHom:1 ;			/* is homozygous - two homologous copies */
  U8 isErr:1 ;			/* is an error */
  U8 isSure:1 ;			/* confident about Het/Hom/Err status */
} CodeHash ;

int compareCodeHash (const void *a, const void *b)
{ U64 ua = ((CodeHash*)a)->hash, ub = ((CodeHash*)b)->hash ;
  if (ua < ub) return -1 ; else if (ua == ub) return 0 ; else return 1 ;
}

int compareCodeRead (const void *a, const void *b)
{ U16 ra = ((CodeHash*)a)->read, rb = ((CodeHash*)b)->read ;
  if (ra < rb) return -1 ; else if (ra == rb) return 0 ; else return 1 ;
}

typedef struct {
  U32 nRead ;			/* actually number of read pairs */
  U32 nHash ;			/* number of unique hashes */
  U32 nCluster ;
  U32 clusterParent:24 ;        /* 1 + offset; 0 if these are original barcodes */
  CodeHash *codeHash ;		/* the unique hashes for this barcode  */
  double pointToMin ;		/* number of good hashes that point to the minimum hash for the cluster - as a fraction of nGoodHashes an indication of hash and so read density */
} BarcodeBlock ;

/* strategy for splitting out clusters into their own barcodes is to leave the original 
   barcode entries in place containing any unclustered hashes, then to follow with new
   clusters, grouped in order of their clusterParent */

/* globals */

MinHash *minHash ;

int hashTableBits = 28 ;	/* default size */
U32 hashTableSize ; 		/* will be 1 << hashTableBits */
U64 hashTableMask ;		/* hashTableSize - 1 */
U32 *hashIndex ;		/* index to be used into arrays */
Array hashCount ;		/* of U32 */
Array hashValue ;		/* of U64 */
Array hashCodes ;		/* of U32*, with hashCount[i] U32 entries */

Array barcodeBlocks ;		/* one BarcodeBlock per barcode */

BOOL isPrintTables = FALSE ;
BOOL isVerbose = FALSE ;
BOOL isInteractive = FALSE ;

int numThreads = 1 ;		/* default to serial - reset  */

FILE *outFile ;			/* initialise to stdout at start of main() */

/********************************************/

void unpackFQB (U32 *u, char *s1, char *q1, char *s2, char *q2)
{
  char *s ; int i,j ;
  if ((s = s1))
    for (i =  0 ; i < 10 ; ++i) for (j = 16 ; j-- ; ) *s++ = (u[i] >> (2*j)) & 3 ;
  if ((s = q1))
    for (i = 10 ; i < 15 ; ++i) for (j = 32 ; j-- ; ) *s++ = (u[i] >> j) & 3 ;
  if ((s = s2))
    for (i = 15 ; i < 25 ; ++i) for (j = 16 ; j-- ; ) *s++ = (u[i] >> (2*j)) & 3 ;
  if ((s = q2))
    for (i = 25 ; i < 30 ; ++i) for (j = 32 ; j-- ; ) *s++ = (u[i] >> j) & 3 ;
}

typedef struct { U64 hash ; int read, pos ; } HashTemp ; /* just for processing */

void seqAddHashes (char *s, int len, Array a, int readCount)
{
  MinHashRCiterator *mi = ranHashRCiterator (minHash, s, len) ;
  U64 hash ; int pos ;
  while (ranHashRCnext (mi, &hash, &pos))
    { HashTemp *h = arrayp(a,arrayMax(a),HashTemp) ;
      h->hash = hash ; h->read = readCount ; h->pos = pos ;
    }
  minHashRCiteratorDestroy (mi) ;
}

int compareHashTemp (const void *a, const void *b)
{ U64 ua = ((HashTemp*)a)->hash, ub = ((HashTemp*)b)->hash ;
  if (ua < ub) return -1 ; else if (ua == ub) return 0 ; else return 1 ;
}

static inline U32 hashIndexFind (U64 hash, int isAdd)
{
  static int nIndex = 0 ;
  U32 index ;
  U64 offset = hash & hashTableMask ;
  U64 diff = ((hash >> hashTableBits) & hashTableMask) | 1 ; /* odd number so coprime to size */
  while ((index = hashIndex[offset]) && (arr(hashValue, index, U64) != hash))
    offset = (offset + diff) & hashTableMask ;
  if (!index && isAdd)
    { index = hashIndex[offset] = arrayMax(hashValue) ;
      array(hashValue, index, U64) = hash ;
      if ((++nIndex >> 2) > hashTableSize) die ("hashTableSize is too small") ;
    }
  return index ;
}

void processBlock (U32 *readData, BarcodeBlock *b)
{
  Array a = arrayCreate (4096, HashTemp) ;
  int i, j, n ;
  char s1[160], s2[160] ;

  for (i = 0 ; i < b->nRead ; ++i, readData += 30)
    { unpackFQB (readData, s1, 0, s2, 0) ;
      seqAddHashes (&s1[23], 127, a, i) ; /* ignore 16bp of barcode and 7bp of spacer */
      seqAddHashes (s2, 150, a, i) ;
    }
  
  arraySort (a, compareHashTemp) ; /* sort then remove duplicate hash values */
  U64 hash, lastHash = arrp(a,0,HashTemp)->hash ;
  for (n = 1, j = 1 ; j < arrayMax(a) ; ++j)
    { hash = arrp(a,j,HashTemp)->hash ;
      if (hash != lastHash)
	arr(a,n++,HashTemp) = arr(a,j,HashTemp) ; lastHash = hash ;
    }

  b->nHash = n ;
  b->codeHash = new (n, CodeHash) ;
  for (i = 0 ; i < n ; ++i)
    { int index = hashIndexFind (arrp(a,i,HashTemp)->hash, TRUE) ;
      ++array(hashCount, index, U32) ;
      b->codeHash[i].hash = index ;
      b->codeHash[i].read = arrp(a,i,HashTemp)->read ;
    }

  qsort (b->codeHash, n, sizeof(CodeHash), compareCodeHash) ; /* sort so hashes are in order */
  
  arrayDestroy (a) ;
}

void readFQB (FILE *f, int N, int chunkSize)
{
  int nReads = 0 ;	/* this group for statistics accumulation */
  U64 nHashes = 0 ;
  
  printf ("  reading and processing sorted fqb file with chunkSize %d", chunkSize) ;
  if (N) printf (", first %d records", N) ;
  printf ("\n") ;
  
  U32 *uBuf = new (chunkSize*30, U32) ; /* read and process readPairs in chunks */
  U32 *uStart = uBuf ;			/* beginning of records for current block */
  U32 barcode = 0 ;
  BarcodeBlock *b = arrayp(barcodeBlocks, 1, BarcodeBlock) ; /* start at 1 */
  b->nRead = 0 ;
  while (!N || nReads < N)
    { if (b->nRead) { memmove (uBuf, uStart, 120*b->nRead) ; uStart = uBuf ; }

      U32 *u = &uBuf[30*b->nRead] ;
      int thisChunk = chunkSize - b->nRead ; if (thisChunk <= 0) die ("chunkSize too small") ;
      if (N && (nReads + thisChunk > N)) { thisChunk = N - nReads ; }
      int nRec = fread (u, 120, thisChunk, f) ;
      if (!nRec) { if (feof (f)) break ; else die ("file read problem") ; } /* loop exit here */
      
      int i ;
      if (!barcode) barcode = u[0] ; /* needed to start barcode matching */
      for (i = 0 ; i < nRec ; ++i)
	if (u[30*i] == barcode) ++b->nRead ;
	else
	  { processBlock (uStart, b) ;
	    nHashes += b->nHash ;
	    b = arrayp(barcodeBlocks, arrayMax(barcodeBlocks), BarcodeBlock) ; b->nRead = 1 ;
	    barcode = u[30*i] ; uStart = &u[30*i] ;
	  }

      nReads += nRec ;
    }
  fclose (f) ;

  fprintf (outFile, "  read %d read pair records for %d barcodes, mean %.2f read pairs per barcode\n",
	  nReads, arrayMax(barcodeBlocks)-1, nReads/(double)(arrayMax(barcodeBlocks)-1)) ;
  fprintf (outFile, "  created %ld hashes, mean %.2f hashes per read pair, %.2f per barcode\n",
	  (long)nHashes, nHashes/(double)nReads, nHashes/(double)(arrayMax(barcodeBlocks)-1)) ;
  if (outFile != stdin)
    { printf ("  read %d read pair records for %d barcodes, mean %.2f read pairs per barcode\n",
	  nReads, arrayMax(barcodeBlocks)-1, nReads/(double)(arrayMax(barcodeBlocks)-1)) ;
      printf ("  created %ld hashes, mean %.2f hashes per read pair, %.2f per barcode\n",
	  (long)nHashes, nHashes/(double)nReads, nHashes/(double)(arrayMax(barcodeBlocks)-1)) ;
    }
}

/*****************************************************************/

static U32 HASHFILE_VERSION = 1 ;
static U16 CODEHASH_SIZE = sizeof(CodeHash) ;
static U16 BCBLOCK_SIZE = sizeof(BarcodeBlock) ;

void writeHashFile (FILE *f)
{
  if ((fwrite ("10XH", 4, 1, f) != 1) || (fwrite (&HASHFILE_VERSION, 4, 1, f) != 1) ||
      (fwrite (&CODEHASH_SIZE, 2, 1, f) != 1) || (fwrite (&BCBLOCK_SIZE, 2, 1, f) != 1) ||
      (fwrite (&hashTableBits, 4, 1, f) != 1)) die ("write fail 1") ;
  if (fwrite (hashIndex, sizeof(U32), hashTableSize, f) != hashTableSize) die ("write fail 2") ;
  if (!arrayWrite (hashValue, f)) die ("failed to write hashValue array") ;
  if (!arrayWrite (hashCount, f)) die ("failed to write hashCount array") ;
  if (!arrayWrite (barcodeBlocks, f)) die ("failed to write barcodeBlocks array") ;
  int i ;
  for (i = 1 ; i < arrayMax(barcodeBlocks) ; ++i)
    { BarcodeBlock *b = arrp(barcodeBlocks, i, BarcodeBlock) ;
      if (fwrite (b->codeHash, sizeof(CodeHash), b->nHash, f) != b->nHash) die ("write fail 3") ;
    }
  fclose (f) ;

  fprintf (outFile, "  wrote %d hash table entries and %d barcode blocks\n",
	   hashTableSize, arrayMax(barcodeBlocks)) ;
  if (outFile != stdout)
    printf ("  wrote %d hash table entries and %d barcode blocks\n",
	    hashTableSize, arrayMax(barcodeBlocks)) ;
}

void readHashFile (FILE *f)
{
  int B ; U32 version ; U16 codeHashSize, bcBlockSize ;
  char name[5] = "abcd" ;
  if ((fread (name, 4, 1, f) != 1) || (fread (&version, 4, 1, f) != 1) ||
      (fread (&codeHashSize, 2, 1, f) != 1) || (fread (&bcBlockSize, 2, 1, f) != 1))
    die ("read fail 0") ;
  if (strcmp (name, "10XH")) die ("not a 10X hash file") ;
  if (version != HASHFILE_VERSION)
    die ("hash file version mismatch: file %d != code %d", version, HASHFILE_VERSION) ;
  if (codeHashSize != CODEHASH_SIZE)
    die ("CodeHash structure size mismatch: file %d != code %d", codeHashSize, CODEHASH_SIZE) ;
  if (bcBlockSize != BCBLOCK_SIZE)
    die ("BarcodeBlock structure size mismatch: file %d != code %d", bcBlockSize, BCBLOCK_SIZE) ;
  if (fread (&B, 4, 1, f) != 1) die ("read fail 1") ;
  if (B != hashTableBits) die ("incompatible hash table size: rerun with -B %d", B) ;
  if (fread (hashIndex, sizeof(U32), hashTableSize, f) != hashTableSize) die ("read fail 2") ;
  if (!(hashValue = arrayRead (f))) die ("failed to read hashValue array") ;
  if (!(hashCount = arrayRead (f))) die ("failed to read hashCount array") ;
  if (!(barcodeBlocks = arrayRead (f))) die ("failed to read barcodeBlocks array") ;
  int i ;
  U64 nReads = 0, nHashes = 0 ;
  for (i = 1 ; i < arrayMax(barcodeBlocks) ; ++i)
    { BarcodeBlock *b = arrayp(barcodeBlocks, i, BarcodeBlock) ;
      nReads += b->nRead ; nHashes += b->nHash ;
      b->codeHash = new(b->nHash, CodeHash) ;
      if (fread (b->codeHash, sizeof(CodeHash), b->nHash, f) != b->nHash) die ("read fail 3") ;
    }
  fclose (f) ;

  fprintf (outFile, "  read %ld hashes for %ld reads in %d barcode blocks\n",
	  (long)nHashes, (long)nReads, arrayMax(barcodeBlocks)) ;
  if (outFile != stdout)
    printf ("  read %ld hashes for %ld reads in %d barcode blocks\n",
	    (long)nHashes, (long)nReads, arrayMax(barcodeBlocks)) ;
}

void fillHashTable (void)
{
  int i, j ;
  long nHashes = 0 ;
  int nHashCount = arrayMax(hashCount) ;

  hashCodes = arrayCreate(nHashCount,U32*) ;
  arrayMax(hashCodes) = nHashCount ;
  
  U32 *count = arrp(hashCount,1,U32) ;
  for (i = 1 ; i < nHashCount ; ++i, ++count)
    if (*count)
      { arr(hashCodes,i,U32*) = new(*count, U32) ;
	nHashes += *count ;
      }

  U32 *hashN = new0 (nHashCount, U32) ;
  for (i = 1 ; i < arrayMax(barcodeBlocks) ; ++i)
    { BarcodeBlock *b = arrp(barcodeBlocks, i, BarcodeBlock) ;
      CodeHash *c = b->codeHash ;
      for (j = 0 ; j < b->nHash ; ++j, ++c)
	arr(hashCodes,c->hash,U32*)[hashN[c->hash]++] = i ;
    }
  fprintf (outFile, "  filled hash table: %ld hashes from %d barcodes in %d bins\n",
	  nHashes, arrayMax(barcodeBlocks), nHashCount) ;
#ifdef CHECK
  for (i = 0 ; i < nHashCount ; ++i)
    if (hashN[i] != arr(hashCount,i,U32))
      die ("hashN[%d] = %d != hashCount[%d] = %d", i, hashN[i], i, arr(hashCount,i,U32)) ;
#endif
  free (hashN) ;
}

/********************************************/

void histogramReport (char* prefix, Array a)
{
  U64 sum = 0, total = 0 ;
  int i ;
  for (i = 0 ; i < arrayMax(a) ; ++i) { sum += arr(a,i,int) ; total += i*arr(a,i,int) ; }
  U64 partSum = 0, partTotal = 0, max = 0, massMax = 0 ;
  int median = 0, massMedian = 0, n99 = 0, nMass99 = 0, mode, massMode ;
  int t50 = sum*0.5, tMass50 = total*0.5, t99 = sum*0.99, tMass99 = total*0.99 ;
  for (i = 0 ; i < arrayMax(a) ; ++i)
    { int n = arr(a,i,int) ;
      partSum += n ;
      partTotal += i*n ;
      fprintf (outFile, "%s_HIST %6d %d %.4f %.4f\n",
	       prefix, i, n, partSum / (double)sum, partTotal / (double)total) ;
      if (n > max) { mode = i ; max = n ; }
      if (i*n > massMax) { massMode = i ; massMax = i*n ; }
      if (partSum > t50 && !median) median = i ;
      if (partTotal > tMass50 && !massMedian) massMedian = i ;
      if (partSum > t99 && !n99) n99 = i ;
      if (partTotal > tMass99 && !nMass99) nMass99 = i ; 
    }
  fprintf (outFile, "%s_STATS MEAN %.1f", prefix, total/(double)sum) ;
  fprintf (outFile, "  MODE %d  MEDIAN %d  PERCENT99 %d", mode, median, n99) ;
  fprintf (outFile, "  MASS_MODE %d  N50 %d  N99 %d\n", massMode, massMedian, nMass99) ;
}

void hashCountHist (void)
{
  if (!arrayMax(hashCount)) { fprintf (stderr, "  no hash list to print stats for\n") ; return ; }
  Array a = arrayCreate (1024, int) ;
  int i ;
  for (i = 0 ; i < arrayMax(hashCount) ; ++i)
    ++array(a, arr(hashCount, i, U32), int) ;
  histogramReport ("HASH_COUNT", a) ;
  arrayDestroy (a) ;
}

void codeSizeHist (void)
{
  if (!barcodeBlocks || !arrayMax (barcodeBlocks))
    { fprintf (stderr, "  no barcodes to print stats for\n") ; return ; }
  Array hashHist = arrayCreate (1024, int) ;
  Array clusterHist = arrayCreate (1024, int) ;
  int i ;
  for (i = 0 ; i < arrayMax(barcodeBlocks) ; ++i)
    { BarcodeBlock *b = arrp (barcodeBlocks, i, BarcodeBlock) ;
      ++array(hashHist, b->nHash, int) ;
      ++array(clusterHist, b->nCluster, int) ;
    }
  histogramReport ("CODE_SIZE", hashHist) ;
  if (arrayMax(clusterHist) > 1) histogramReport ("CODE_CLUSTER", clusterHist) ;
}

/************************************************************/

typedef struct {
  I16 chr ;		/* negative numbers for multiple occurrences in one genome */
  U16 pos ;		/* in units of 10k */
} CribInfo ;

#define CRIB_ERR 0
#define CRIB_HTA 1
#define CRIB_HTB 2
#define CRIB_HOM 3
#define CRIB_MUL 4
#define CRIBTYPE_MAX 5
static char* cribTypeName[] = {"err","htA","htB","hom","mul"} ;
static char cribTypeChar[] = "eabhm" ;

CribInfo *crib = 0 ;	/* of CribInfo over hash index */
U8 *cribType = 0 ;	/* over hash index with values CRIB_* */
int cribChrMax = 0 ;

#include "readseq.h"

static void cribAddGenome (FILE *f, CribInfo *crib)
{
  char *seq ;
  int len, chr = 0, pos ;
  int nAbsent = 0, nPresent = 0 ;
  U32 index ;
  dna2indexConv['N'] = dna2indexConv['n'] = 0 ; /* to get 2-bit encoding */
  while (readSequence (f, dna2indexConv, &seq, 0, 0, &len))
    { MinHashRCiterator *mi = ranHashRCiterator (minHash, seq, len) ;
      U64 hash ; int pos ;
      ++chr ; if (chr >= cribChrMax) cribChrMax = chr+1 ;
      while (ranHashRCnext (mi, &hash, &pos))
	if ((index = hashIndexFind (hash, FALSE)))
	  { CribInfo *c = &crib[index] ;
	    if (!c->chr) { c->chr = chr ; c->pos = pos >> 10 ; }
	    else if (c->chr > 0) c->chr = -1 ;
	    else --c->chr ;
	    ++nPresent ;
	  }
	else
	  ++nAbsent ;
      minHashRCiteratorDestroy (mi) ;
    }
  fclose (f) ;

  fprintf (outFile, "  read %d known and %d unknown hashes from %d sequences in crib genome\n",
	  nPresent, nAbsent, chr) ;
  if (outFile != stdout)
    printf ("  read %d known and %d unknown hashes from %d sequences in crib genome\n",
	    nPresent, nAbsent, chr) ;
}

void printArrayStats (Array a)
{
  int i, *ai = arrp(a,0,int), sum = 0, min = -1, max = arrayMax(a)-1 ;
  double total = 0 ;
  for (i = 0 ; i < arrayMax(a) ; ++i, ++ai)
    if (*ai)
      { sum += *ai ; total += *ai * i ;
	if (min == -1) min = i ;
      }
  fprintf (outFile, "  %d mean %.1f min %d max %d\n", sum, total/sum, min, max) ;
}

void cribBuild (FILE *f1, FILE *f2)
{
  crib = new0 (arrayMax(hashValue), CribInfo) ; cribAddGenome (f1, crib) ;
  CribInfo *crib2 = new0 (arrayMax(hashValue), CribInfo) ; cribAddGenome (f2, crib2) ;

  cribType = new0 (arrayMax(hashValue), U8) ;
  Array aHom = arrayCreate (256, int), aHet = arrayCreate (256, int) ;
  Array aErr = arrayCreate (256, int), aMul = arrayCreate (256, int) ;
  int i ;
  for (i = 1 ; i < arrayMax(hashValue) ; ++i)
    if (crib[i].chr == 0 && crib2[i].chr == 0)
      { cribType[i] = CRIB_ERR ; ++array(aErr, arr(hashCount, i, int), int) ; }
    else if (crib[i].chr > 0 && crib2[i].chr > 0)
      { cribType[i] = CRIB_HOM ; ++array(aHom, arr(hashCount, i, int), int) ; }
    else if (crib[i].chr < 0 || crib2[i].chr < 0)
      { cribType[i] = CRIB_MUL ; ++array(aMul, arr(hashCount, i, int), int) ;
	if (crib2[i].chr < crib[i].chr) crib[i].chr = crib2[i].chr ;
      }
    else
      { cribType[i] = crib[i].chr ? CRIB_HTA : CRIB_HTB ;
	++array(aHet, arr(hashCount, i, int), int) ;
	if (!crib[i].chr) { crib[i].chr = crib2[i].chr ; crib[i].pos = crib2[i].pos ; }
      }

  fprintf (outFile, "  crib matches\n") ;
  fprintf (outFile, "    hom  ") ; printArrayStats (aHom) ;
  fprintf (outFile, "    het  ") ; printArrayStats (aHet) ;
  fprintf (outFile, "    mul ") ; printArrayStats (aMul) ;
  fprintf (outFile, "    err ") ; printArrayStats (aErr) ;

  if (isPrintTables)
    { fprintf (outFile, "CRIB_TABLE      i   err         het          hom         mul\n") ;
      for (i = 1 ; i < 256 ; ++i)
	fprintf (outFile, "CRIB_TABLE      %4d%12d%12d%12d%12d\n", i,
		 arrayMax(aErr)>i?arr(aErr,i,int):0, arrayMax(aHet)>i?arr(aHet,i,int):0,
		 arrayMax(aHom)>i?arr(aHom,i,int):0, arrayMax(aMul)>i?arr(aMul,i,int):0) ;
    }

  arrayDestroy (aHom) ; arrayDestroy (aHet) ; arrayDestroy (aMul) ; arrayDestroy (aErr) ;
  free (crib2) ;
}

char* cribText (int x)
{ static char text[32] ;
  char *s = text ;
  s += sprintf (s, "%d", x) ;
  if (crib) s += sprintf (s, ":%s", cribTypeName[cribType[x]]) ;
  if (crib && cribType[x] > CRIB_ERR && cribType[x] < CRIB_MUL)
    s += sprintf (s, "_%d.%d", crib[x].chr, crib[x].pos) ;
  sprintf (s, "-%d", arr(hashCount, x, U32)) ;
  return text ;
}

/************************************************************/

static BOOL *hashWithinRange = 0 ;
static int hashRangeMin = 0, hashRangeMax = 0 ;

BOOL hashWithinRangeBuild (int min, int max)
{ 
  if (hashWithinRange && min == hashRangeMin && max == hashRangeMax) return FALSE ;
  if (!hashWithinRange) hashWithinRange = new0 (arrayMax(hashCount), BOOL) ;
  int i, n ;
  for (i = 0 ; i < arrayMax(hashCount) ; ++i)
    { n = arr(hashCount, i, U32) ;
      if (n >= min && n < max) hashWithinRange[i] = TRUE ;
    }
  hashRangeMin = min ; hashRangeMax = max ;
  return TRUE ;
}

Array hashNeighbours (int x)
{ /* return array of hashes that share between min and max codes with x */
  /* abuse CodeHash structure for this array, storing the number of shared codes in ->read */
  if (!hashWithinRange) die ("hashNeighbours() called without setting hashCountRange") ;
  Array ch = arrayCreate (1000, CodeHash) ;
  int i,j ;
  U32 xCount = arr(hashCount,x,U32) ; if (!xCount) return 0 ;
  for (i = 0 ; i < xCount ; ++i)
    { BarcodeBlock *b = arrp(barcodeBlocks, arr(hashCodes,x,U32*)[i], BarcodeBlock) ;
      CodeHash *c = b->codeHash ;
      for (j = 0 ; j < b->nHash ; ++j, ++c)
	if (c->hash != x && hashWithinRange[c->hash])
	  array(ch, arrayMax(ch), CodeHash) = *c ;
    }
  arraySort (ch, compareCodeHash) ;
  Array shared = arrayCreate (256, CodeHash) ; /* abuse structure esp. read field */
  CodeHash *last = arrayp (shared, 0, CodeHash) ;
  last->hash = arrp(ch, 0, CodeHash)->hash ; last->read = 1 ;
  for (i = 1 ; i < arrayMax(ch) ; ++i)
    if (arrp(ch, i, CodeHash)->hash == last->hash) ++last->read ;
    else
      { last = arrayp (shared, arrayMax(shared), CodeHash) ;
	last->hash = arrp(ch, i, CodeHash)->hash ; last->read = 1 ;
      }
  arrayDestroy (ch) ;
  return shared ;
}

Array countHashNeighbours (int x)
{
  U32 xCount = arr(hashCount,x,U32) ; if (!xCount) return 0 ;
  HASH h = hashCreate (500000) ;
  Array a = arrayCreate (500000, int) ;
  int i,j ;
  for (i = 0 ; i < xCount ; ++i)
    { BarcodeBlock *b = arrp(barcodeBlocks, arr(hashCodes,x,U32*)[i], BarcodeBlock) ;
      CodeHash *c = b->codeHash ;
      for (j = 0 ; j < b->nHash ; ++j, ++c)
	if (c->hash != x && hashWithinRange[c->hash])
	  ++array (a, hashAdd(h,HASH_INT(c->hash)), int) ;
    }
  Array b = arrayCreate (256, int) ;
  for (i = 1 ; i < arrayMax(a) ; ++i) ++array (b, arr(a,i,int), int) ;
  arrayDestroy (a) ; hashDestroy (h) ;
  return b ;
}

void exploreHash (int x)
{
  if (!hashWithinRange)
    { fprintf (stderr, "exploreHash called without hashCountRange\n") ; return ; }
  Array shared = hashNeighbours (x) ; if (!shared) return ;
  fprintf (outFile, "  %d hashes sharing codes with %s\n", arrayMax(shared), cribText(x)) ;
  arraySort (shared, compareCodeRead) ;
  int i, count = 1, n = 0 ;
  CodeHash *d = arrp(shared, 0, CodeHash) ;
  for (i = 0 ; i <= arrayMax(shared) ; ++i, ++d)
    if (d->read != count || i == arrayMax(shared))
      { fprintf (outFile, "    %d sharing %d codes", n, count) ;
	if (n > 7) { fprintf (outFile, "  ...") ; n = 7 ; }
	while (n) { fprintf (outFile, " %s", cribText(d[-n].hash)) ; --n ; }
	fprintf (outFile, "\n") ;
	count = d->read ; n = 1 ;
      }
    else ++n ;
  
  arrayDestroy (shared) ;
}

void doubleShared (int x1, int x2)
{
  if (!hashWithinRange)
    { fprintf (stderr, "doubleShared called without hashCountRange\n") ; return ; }
  Array sh1 = hashNeighbours (x1) ; if (!sh1) return ;
  Array sh2 = hashNeighbours (x2) ; if (!sh2) { arrayDestroy (sh1) ; return ; }

  int i1 = 0, i2 = 0 ;
  CodeHash *c1 = arrp(sh1, 0, CodeHash), *c2 = arrp(sh2, 0, CodeHash) ;
  while (i1 < arrayMax(sh1) && i2 < arrayMax(sh2))
    if (c1->hash == c2->hash && c1->read > 2 && c2->read > 2)
      { fprintf (outFile, "  hash %-8d code %-6d n1 %-3d n2 %-3d crib %s\n", c1->hash,
		*arr(hashCodes,c1->hash,U32*), c1->read, c2->read, cribText(c1->hash)) ;
	++c1 ; ++i1 ; ++c2 ; ++i2 ;
      }
    else if (c1->hash < c2->hash) { ++c1 ; ++i1 ; }
    else { ++c2 ; ++i2 ; }

  arrayDestroy (sh1) ; arrayDestroy (sh2) ;
}

void hashInfo (int hMin, int hMax, int hIncrement)
{
  if (!hashWithinRange)
    { fprintf (stderr, "hashInfo called without hashCountRange\n") ; return ; }
  int x ;
  for (x = hMin ; x < hMax ; x += hIncrement)
    { fprintf (outFile, "HASH_INFO  %s", cribText(x)) ;
      if (hashWithinRange[x])
	{ Array shared = hashNeighbours (x) ;
	  arraySort (shared, compareCodeRead) ;
	  CodeHash *ch = arrp(shared, arrayMax(shared)-1, CodeHash) ;
	  fprintf (outFile, " max share %d with %s", ch->read, cribText (ch->hash)) ; 
	  arrayDestroy (shared) ;
	}
      fputc ('\n', outFile) ;
    }
}

/***************************************/

void errorFix (int hashMin, int hashMax)
{
  if (!crib) die ("need to set crib") ;
  U32 x ; int i ;
  int low[CRIBTYPE_MAX], nt[CRIBTYPE_MAX] ;
  double min[CRIBTYPE_MAX], max[CRIBTYPE_MAX], sum[CRIBTYPE_MAX] ;
  for (i = 0 ; i < CRIBTYPE_MAX ; ++i)
    { low[i] = 0 ; nt[i] = 0 ; min[i] = 1.0 ; max[i] = 0.0 ; sum[i] = 0.0 ; }
  for (x = hashMin ; x < hashMax ; ++x)
    { int type = cribType[x] ;
      int xCount = arr(hashCount, x, U32) ; if (xCount < 5) { ++low[type] ; continue ; }
      Array a = countHashNeighbours (x) ;
      --*arrp(a,arrayMax(a)-1,int) ; /* remove top match - often wrong */
      while (arrayMax(a) && !arr(a,arrayMax(a)-1,int)) --arrayMax(a) ;
      if (arrayMax(a) > 4)
	{ int n = 0 ; U64 sumSq = 0 ; int *ai = arrp(a,0,int) ;
	  for (i = 0 ; i < 4 ; ++i) n += *ai++ ;
	  for ( ; i < arrayMax(a) ; ++i) { sumSq += *ai * (i-3)*(i-3) ; n += *ai++ ; }
	  double score = sumSq/((double)n*(xCount-3)) ;
	  ++nt[type] ; sum[type] += score ;
	  if (score < min[type]) min[type] = score ;
	  if (score > max[type]) max[type] = score ;
	  if ((type == CRIB_ERR && score > 0.00005) || (type != CRIB_ERR && score < 0.00005))
	    { printf ("  %s count %d max %d score %.6f ",
		      cribText(x), xCount, arrayMax(a)-1, score) ;
	      for (i = 5 ; i < arrayMax(a) ; ++i) printf (" %d", arr(a,i,int)) ;
	      putchar ('\n') ;
	    }
	}
      else
	++low[type] ;
      arrayDestroy (a) ;
    }
  for (i = 0 ; i < CRIBTYPE_MAX ; ++i)
    printf ("  %s low %d n %d min %.6f av %.6f max %.6f\n",
	    cribTypeName[i], low[i], nt[i], min[i], sum[i] / (nt[i]?nt[i]:1), max[i]) ;
  timeUpdate (stdout) ;
}


void shareScan (int countMin, int countMax)
{
  U32 x ;
  int i, j, k, nTot = 0 ; 
  int *nBin = new0 (countMax-countMin, int) ;
  U32 *xBin = new (10*(countMax-countMin), U32) ;
  Array *aBin = new (10*(countMax-countMin), Array) ;
  for (x = countMax ; x < arrayMax(hashCount) ; x += 100)
    { int xCount = arr(hashCount, x, U32) ;
      if (xCount < countMin || xCount >= countMax) continue ; /* out of range */
      int bin = xCount - countMin ;
      if (nBin[bin] == 10) continue ;                         /* already have enough */
      xBin[10*bin + nBin[bin]] = x ;
      Array a = countHashNeighbours (x) ;
      --*arrp(a,arrayMax(a)-1,int) ;                          /* remove top match; often wrong */
      while (arrayMax(a) && !arr(a,arrayMax(a)-1,int)) --arrayMax(a) ;
      aBin[10*bin + nBin[bin]++] = a ;
      if (++nTot == 10*(countMax-countMin)) break ;           /* done */
    }

  for (i = 0 ; i < countMax-countMin ; ++i)
    for (j = 0 ; j < nBin[i] ; ++j)
      { fprintf (outFile,"SHARE_SCAN %3d %24s ", i+countMin, cribText(xBin[10*i+j])) ;
	Array a = aBin[10*i+j] ;
	for (k = 1 ; k < arrayMax(a) ; ++k) fprintf (outFile," %d", arr(a,k,int)) ;
	fputc ('\n', outFile) ;
      }

  timeUpdate (stdout) ;
}

/***************************************/

U16 **goodHashes = 0 ;
int *nGoodHashes = 0 ;

void goodHashesBuild (void)
{
  if (!hashWithinRange) die ("cluster code called without setting hashCountRange") ;
  goodHashes = new (arrayMax(barcodeBlocks), U16*) ;
  nGoodHashes = new (arrayMax(barcodeBlocks), int) ;
  Array a = arrayCreate (4096, U16) ;
  int c, i ;
  
  for (c = 0 ; c < arrayMax(barcodeBlocks) ; ++c)
    { BarcodeBlock *cb = arrp(barcodeBlocks, c, BarcodeBlock) ;
      if (cb->nHash > U16MAX)
	{ nGoodHashes[c] = 0 ; goodHashes[c] = new(1,U16) ;
	  fprintf (stderr, "ignoring barcode %d - too many hashes %d > %d\n",
		   c, cb->nHash, U16MAX) ;
	  continue ;
	}
      CodeHash *ch = cb->codeHash ;
      arrayMax(a) = 0 ;
      for (i = 0 ; i < cb->nHash ; ++i, ++ch)
	if (hashWithinRange[ch->hash]) array(a, arrayMax(a), U16) = i ;
      goodHashes[c] = new(arrayMax(a),U16) ;
      memcpy (goodHashes[c], arrp(a,0,U16), arrayMax(a)*sizeof(U16)) ;
      nGoodHashes[c] = arrayMax(a) ;
    }

  printf ("  made goodHashes arrays for hash range %d to %d\n  ", hashRangeMin, hashRangeMax) ;
  timeUpdate (stdout) ; fflush (stdout) ;
}

/*************************/

void codeClusterFind (int code, int clusterThreshold) /* clusters this barcode's good hashes */
{
  BarcodeBlock *b = arrp(barcodeBlocks, code, BarcodeBlock), *c ;
  CodeHash *cb = b->codeHash, *cc ;
  int i, j ;
  U32 x, nc ;
  int *minShare = new0 (arrayMax(barcodeBlocks), int) ; /* lowest hash shared by code and j */
  U16 *g, *gi, *gc, *gcn ;
  Array clusterMin = arrayCreate (64, int) ;

  int n = nGoodHashes[code] ; if (!n) return ;
  int *minShareCount = new (n, int) ;

  for (i = 0 ; i < n ; ++i) b->codeHash[goodHashes[code][i]].cluster = 0 ; /* wipe existing */

  b->nCluster = 0 ;
  b->pointToMin = 0.0 ;
  int nClustered = 0 ;
  for (i = 1 ; i < n ; ++i)
    { x = b->codeHash[goodHashes[code][i]].hash ;
      nc = arr(hashCount, x, U32) ;
      memset (minShareCount, 0, n*sizeof(int)) ;
					       
      for (j = 0 ; j < nc ; ++j)
	{ int cj = arr(hashCodes, x, U32*)[j] ;
	  if (cj == code || !nGoodHashes[cj]) continue ;
	  if (!minShare[cj])
	    { minShare[cj] = i + 1 ; /* default value; +1 to distinguish from 0 */
	      g = goodHashes[code] ; gi = &g[i] ; gc = goodHashes[cj] ;
	      cc = arrp(barcodeBlocks, cj, BarcodeBlock)->codeHash ;
	      U32 hb = cb[*g].hash, hc = cc[*gc].hash ;
	      while (g < gi)
		if (hb < hc) hb = cb[*++g].hash ;
		else if (hb == hc)
		  { minShare[cj] = g-goodHashes[code]+1 ;
		    break ;
		  }
		else hc = cc[*++gc].hash ;
	    }
	  ++minShareCount[minShare[cj]-1] ;
	}

      int msMax = 0, msTot = 0, msBest ;
      g = goodHashes[code] ;
      for (j = 0 ; j < i ; ++j)
	{ if (minShareCount[j] > msMax)  { msBest = j ; msMax = minShareCount[j] ; }
	  msTot += minShareCount[j] ;
	}
      if (msMax >= clusterThreshold)
	{ CodeHash *ch = &b->codeHash[g[msBest]] ;
	  if (!ch->cluster)	/* create a new cluster */
	    { if (++b->nCluster > U8MAX) /* abandon this clustering */
		{ b->nCluster = 0 ; nClustered = 0 ;
		  for (j = 0 ; j < i ; ++j) b->codeHash[g[j]].cluster = 0 ;
		  fprintf (stderr, "    code %d with %d good hashes has too many clusters\n",
			   code, nGoodHashes[code]) ;
		  break ;
		}
	      ch->cluster = b->nCluster ; ++nClustered ;
	      array(clusterMin,b->nCluster,int) = msBest ;
	    }
	  b->codeHash[g[i]].cluster = ch->cluster ; ++nClustered ;
	  b->pointToMin += minShareCount[arr(clusterMin,ch->cluster,int)] / (double)msTot ;
	}
    }
  free (minShare) ;
  free (minShareCount) ;
  arrayDestroy (clusterMin) ;
  if (isVerbose)
    { if (b->nCluster)
	fprintf (outFile, "  code %d with %d hashes, %d good hashes, of which %d cluster into %d raw",
		code, b->nHash, nGoodHashes[code], nClustered, b->nCluster) ;
      else
	fprintf (outFile, "  code %d with %d hashes, %d good hashes, 0 clusters\n",
		code, b->nHash, nGoodHashes[code]) ;
    }
}

void codeClusterReadMerge (int code)
{
  BarcodeBlock *b = arrp(barcodeBlocks, code, BarcodeBlock) ;
  if (!b->nCluster) return ;	/* nothing to do */
  int i, j ;
  int *readMap = new0 (b->nRead, int), *trueCluster = new (b->nCluster+1, int) ;
  int *deadCluster = new0 (b->nCluster+1, int) ;
  for (i = 0 ; i <= b->nCluster ; ++i) trueCluster[i] = i ; /* initialisation */
  for (i = 0 ; i < b->nHash ; ++i)
    { int hashCluster = trueCluster[b->codeHash[i].cluster] ; if (!hashCluster) continue ;
      int readCluster = trueCluster[readMap[b->codeHash[i].read]] ;
      if (hashCluster == readCluster) continue ;
      else if (!readCluster) readMap[b->codeHash[i].read] = hashCluster ;
      else			/* merge clusters hashCluster and readCluster */
	{ if (hashCluster > readCluster)
	    { int t = hashCluster ; hashCluster = readCluster ; readCluster = t ; }
	  for (j = 1 ; j <= b->nCluster ; ++j)
	    if (trueCluster[j] == readCluster) trueCluster[j] = hashCluster ;
	  deadCluster[readCluster] = 1 ;
	}
    }

  /* now compactify trueCluster, reusing deadCluster[] for the map removing dead indices */
  for (j = 1 ; j <= b->nCluster ; ++j) deadCluster[j] = deadCluster[j-1] + 1 - deadCluster[j] ;
  for (j = 1 ; j <= b->nCluster ; ++j) trueCluster[j] = deadCluster[trueCluster[j]] ;
  b->nCluster = deadCluster[b->nCluster] ;
  for (i = 0 ; i < b->nHash ; ++i) b->codeHash[i].cluster = trueCluster[b->codeHash[i].cluster] ;

  free (readMap) ; free(trueCluster) ; free(deadCluster) ;
  if (isVerbose) { fprintf (outFile, " then %d merged clusters\n", b->nCluster) ; fflush (stdout) ; }
}

void codeClusterReport (int codeMin, int codeMax)
{
  int code ;
  typedef struct {
    int i ;
    int n ;
    int nRead ;
    int nt[5] ;
    I16 chr ;
    U16 pMin, pMax ;
    int nBad ;
    int bad ;			/* point to bad hash */
  } ReportInfo ;
  Array info = arrayCreate (256, ReportInfo) ;
  Array readClus = arrayCreate (1024, int) ;
  Array badLink = arrayCreate (8192, int) ;
  ReportInfo *r ;
  int i, c ;
  U64 totalGoodHash = 0 ;
  double totalPointToMin = 0.0 ;

  for (code = codeMin ; code < codeMax ; ++code)
    { BarcodeBlock *b = arrp(barcodeBlocks, code, BarcodeBlock) ;
      arrayReCreate (info, b->nCluster+1, ReportInfo) ;
      arrayReCreate (readClus, b->nRead, int) ;
      arrayReCreate (badLink, b->nHash, int) ;
      int nClusHash = 0 ;
      for (i = 0 ; i < b->nHash ; ++i)
	{ if (!(c = b->codeHash[i].cluster)) continue ;
	  ++nClusHash ;
	  arr(readClus, b->codeHash[i].read, int) = c ;
	  r = arrp(info, c, ReportInfo) ;
	  ++r->n ;
	  U32 bh = b->codeHash[i].hash ;
	  if (crib)
	    { ++r->nt[cribType[bh]] ;
	      if (cribType[bh] > CRIB_ERR && cribType[bh] < CRIB_MUL)
		{ CribInfo ch = crib[bh] ;
		  if (!r->chr) { r->chr = ch.chr ; r->pMin = r->pMax = ch.pos ; }
		  else if (ch.chr == r->chr)
		    { if (ch.pos < r->pMin) r->pMin = ch.pos ;
		      if (ch.pos > r->pMax) r->pMax = ch.pos ;
		    }
		  else { ++r->nBad ; arr(badLink,i,int) = r->bad ; r->bad = i ; }
		}
	    }
	}

      int nClusRead = 0 ;
      for (i = 0 ; i < b->nRead ; ++i)
	if ((c = arr(readClus,i,int))) { ++arrp(info,c,ReportInfo)->nRead ; ++nClusRead ; }
      
      fprintf (outFile, "  CLUSTER_SUMMARY %d nRead %d nHash %d nGoodHash %d nClusHash %d nClusRead %d nCluster %d\n",
	      code, b->nRead, b->nHash, nGoodHashes[code], nClusHash, nClusRead, b->nCluster) ;
      for (i = 1 ; i <= b->nCluster ; ++i)
	{ r = arrp(info, i, ReportInfo) ; if (!r->n) continue ;
	  fprintf (outFile, "    CODE_CLUSTER %d %d : %d reads %d hashes",
		   code, i, r->nRead, r->n) ;
	  if (crib)
	    { fprintf (outFile, " %d hom, %d htA, %d htB, %d mul, %d err", r->nt[CRIB_HOM],
		      r->nt[CRIB_HTA], r->nt[CRIB_HTB], r->nt[CRIB_MUL], r->nt[CRIB_ERR]) ;
	      if (r->chr)
		fprintf (outFile, "  chr %d pos %d %d", r->chr, r->pMin, r->pMax-r->pMin+1) ;
	      if (r->nBad)
		{ fprintf (outFile, "  OTHER %d", r->nBad) ;
		  int x = r->bad, j = 10 ;
		  while (x && j--)
		    { fprintf (outFile, " %s", cribText(b->codeHash[x].hash)) ;
		      x = arr(badLink, x, int) ;
		    }
		}
	    }
	  fputc ('\n', outFile) ;
	}
      totalGoodHash += nGoodHashes[code] ;
      totalPointToMin += b->pointToMin ;
    }

  if (totalGoodHash)
    fprintf (outFile, "  MIN_POINT_DENSITY %.3f\n", totalPointToMin / totalGoodHash) ;
  
  arrayDestroy (info) ; arrayDestroy (readClus) ;
}

/************************************************************/

void clusterSplitCodes (void)
{
  int i, j, nClusters = 0 ;
  int nCodes = arrayMax (barcodeBlocks) ; 
  for (i = 0 ; i < nCodes ; ++i) nClusters += arrp(barcodeBlocks, i, BarcodeBlock)->nCluster ;
  Array newBlocks = arrayCreate (nCodes+nClusters, BarcodeBlock) ;
  arrayMax(newBlocks) = nCodes + nClusters ;
  BarcodeBlock *old = arrp(barcodeBlocks, 0, BarcodeBlock) ;
  BarcodeBlock *new1 = arrp(newBlocks, 0, BarcodeBlock) ;
  BarcodeBlock *new2 = arrp(newBlocks, nCodes-1, BarcodeBlock) ; //-1 allows for ->cluster offset
  Array clusterHashCount = arrayCreate (256, int) ; /* how many hashes in each cluster */
  Array readMap = arrayCreate (4096, int) ;
  for (i = 0 ; i < nCodes ; ++i, ++old)
    if (old->nCluster)
      { arrayReCreate (clusterHashCount, old->nCluster+1, int) ; /* zeroes array */
	CodeHash *c = old->codeHash ;
	for (j = 0 ; j < old->nHash ; ++j, ++c) ++arr(clusterHashCount,c->cluster,int) ;
	// new1->nCluster = 0 ; // not needed because already 0
	new1->nHash = arr(clusterHashCount,0,int) ;
	if (new1->nHash) new1->codeHash = new0(new1->nHash, CodeHash) ;
	for (j = 1 ; j < old->nCluster+1 ; ++j)
	  { new2[j].nHash = arr(clusterHashCount,j,int) ;
	    if (new2[j].nHash) new2[j].codeHash = new0(new2[j].nHash, CodeHash) ;
	    new2[j].clusterParent = i+1 ;
	    // new2[j].nCluster = 0 ;   // not needed because already 0
	  }
	arrayReCreate (clusterHashCount, old->nCluster+1, int) ;
	arrayReCreate (readMap, old->nRead+1, int) ;
	for (j = 0, c = old->codeHash ; j < old->nHash ; ++j, ++c)
	  { int clus = c->cluster ; c->cluster = 0 ;
	    if (!array(readMap,c->read,int)) arr(readMap,c->read,int) = ++new2[clus].nRead ;
	    // previous line requires that all hashes from same read are in same cluster
	    c->read = arr(readMap,c->read,int) - 1 ;
	    if (clus) new2[clus].codeHash[arr(clusterHashCount,clus,int)++] = *c ;
	    else new1->codeHash[arr(clusterHashCount,0,int)++] = *c ;
	  }
	free (old->codeHash) ;
	++new1 ;
	new2 += old->nCluster ;
      }
    else
      *new1++ = *old ; // also carries over existing codeHash
  
  fprintf (outFile, "  made %d additional new barcodes from clusters in %d original barcodes\n",
	   nClusters, arrayMax(barcodeBlocks)) ;
  if (outFile != stdout)
    printf ("  made %d additional new barcodes from clusters in %d original barcodes\n",
	    nClusters, arrayMax(barcodeBlocks)) ;
  arrayDestroy (barcodeBlocks) ;
  barcodeBlocks = newBlocks ;
  arrayDestroy (clusterHashCount) ; arrayDestroy (readMap) ;
  printf ("  cluster timepoint: ") ; timeUpdate (stdout) ;

  // now must rebuild the mapping from hashes to codes
  for (i = 1 ; i < arrayMax(hashCount) ; ++i)
    if (arr(hashCount,i,U32)) free (arr(hashCodes,i,U32*)) ;
  arrayDestroy (hashCodes) ;
  fillHashTable () ;
}

/************************************************************/

void usage (void)
{ fprintf (stderr, "Usage: hash10x <commands>\n") ;
  fprintf (stderr, "Commands can be parameter settings with -x, or operations:\n") ;
  fprintf (stderr, "Be sure to set relevant parameters before invoking an operation!\n") ;
  fprintf (stderr, "   -k <kmer size> [%d]\n", params.k) ;
  fprintf (stderr, "   -w <window> [%d]\n", params.w) ;
  fprintf (stderr, "   -r <random number seed> [%d]\n", params.r) ;
  fprintf (stderr, "   -B <hash index table bitcount> [%d]\n", params.B) ;
  fprintf (stderr, "   -N <num records to read: 0 for all> [%d]\n", params.N) ;
  fprintf (stderr, "   -c <file chunkSize in readPairs> [%d]\n", params.chunkSize) ;
  fprintf (stderr, "   -cT | --clusterThreshold <clusterThreshold> [%d]\n", params.clusterThreshold) ;
  fprintf (stderr, "   -t | --threads <number of threads for parallel ops> [%d]\n", numThreads) ;
  fprintf (stderr, "   -o | --output <output filename> : '-' for stdout\n") ;
  fprintf (stderr, "   --interactive: enter interactive mode: --commands only without --\n") ;
  fprintf (stderr, "   --tables: toggle to print tables for operations while this is set\n") ;
  fprintf (stderr, "   --verbose: toggle for verbose mode\n") ;
  fprintf (stderr, "   --readFQB <sorted fqb input file name>: must have this or readHash\n") ;
  fprintf (stderr, "   --readHash <hash input file name>\n") ;
  fprintf (stderr, "   --writeHash <hash output file name>\n") ;
  fprintf (stderr, "   --hashInfo <start> <end> <skip>: info for hashes in [start,end)\n") ;
  fprintf (stderr, "   --hashExplore <hash>: look for structure around hash\n") ;
  fprintf (stderr, "   --doubleShared <hash1> <hash2>: codes shared with both\n") ;
  fprintf (stderr, "   --hashCountRange <min> <max>: set limits for hash counts for Explore, cluster etc.\n") ;
  fprintf (stderr, "   --cribBuild <genome1.fa> <genome2.fa>: match to genomic hashes\n") ;
  fprintf (stderr, "   --cluster <codeMin> <codeMax>: cluster reads for range of barcodes (1, 0 for all)\n") ;
  fprintf (stderr, "   --clusterReport <codeMin> <codeMax>\n") ;
  fprintf (stderr, "   --clusterSplit\n") ;
  fprintf (stderr, "   --hashStats : distribution of hash counts and summary info\n") ;
  fprintf (stderr, "   --codeStats : distribution of barcode/cluster sizes and summary info\n") ;
  fprintf (stderr, "   --help : print this usage message\n") ;
  fprintf (stderr, "   --quit : end interactive input and exit program\n") ;
  fprintf (stderr, "   --exit : end interactive input and exit program\n") ;
}

void initialise (int k, int w, int r, int B)
{
  srandom (r) ;

  if (k <= 0 || w <= 0) die ("k %d, w %d must be > 0; run without args for usage", k, w) ;
  minHash = minHashCreate (k, w) ;

  hashTableBits = B ;
  if (hashTableBits < 20 || hashTableBits > 30)
    die ("hashTableBits %d out of range 20-30", hashTableBits) ;
  hashTableSize = (size_t)1 << hashTableBits ;
  hashTableMask = hashTableSize - 1 ;
  hashIndex = new0 (hashTableSize, U32) ;
  hashValue = arrayCreate (1 << 20, U64) ; array(hashValue,0,U64) = 0 ; // so true indexes != 0
  hashCount = arrayCreate (1 << 20, U32) ;

  fprintf (outFile, "hash10x initialised with k = %d, w = %d, random seed = %d, hashtable bits = %d\n",
	  k, w, r, B) ;
}

int main (int argc, char *argv[])
{
  --argc ; ++argv ;		/* eat program name */

  outFile = stdout ;
  
  timeUpdate (stdout) ;		/* initialise timer */

  /* command line parameters */
  params.k = 21 ;
  params.w = 31 ;
  params.r = 17 ;
  params.B = 28 ;
  params.N = 0 ;
  params.chunkSize = 100000 ;
  params.clusterThreshold = 5 ;
#ifdef OMP
  numThreads = omp_get_max_threads () ;
  omp_set_num_threads (numThreads) ;
#endif
  
  if (!argc) usage () ;

#ifdef SIZES
  printf ("sizeof BarcodeBlock is %d\n", sizeof(BarcodeBlock)) ;
  printf ("sizeof CodeHash is %d\n", sizeof(CodeHash)) ;
  exit (0) ;
#endif

  barcodeBlocks = arrayCreate (1200, BarcodeBlock) ;

  FILE *f ;
  char inputBuffer[1024] ;	/* use these for interactive mode */
  char *inputArgs[32] ;
  Array inputArray = arrayCreate (1024, char) ;
  
  while (argc) {
    if (**argv != '-')
      die ("option/command %s does not start with '-': run without arguments for usage", *argv) ;
    if (!isInteractive)
      { fprintf (outFile, "COMMAND %s", *argv) ;
	int i ;
	for (i = 1 ; i < argc && *argv[i] != '-' ; ++i) fprintf (outFile, " %s", argv[i]) ;
	fputc ('\n', outFile) ;
	if (outFile != stdout)
	  { printf ("COMMAND %s", *argv) ;
	    int i ;
	    for (i = 1 ; i < argc && *argv[i] != '-' ; ++i) fprintf (outFile, " %s", argv[i]) ;
	    putchar ('\n') ;
	  }
      }
#define ARGMATCH(x,n)	(!strcmp (*argv, x) && argc >= n && (argc -= n, argv += n))
    if (ARGMATCH("-k",2)) params.k = atoi(argv[-1]) ;
    else if (ARGMATCH("-w",2)) params.w = atoi(argv[-1]) ;
    else if (ARGMATCH("-r",2)) params.r = atoi(argv[-1]) ;
    else if (ARGMATCH("-B",2)) params.B = atoi(argv[-1]) ;
    else if (ARGMATCH("-N",2)) params.N = atoi(argv[-1]) ;
    else if (ARGMATCH("-c",2)) params.chunkSize = atoi(argv[-1]) ;
    else if (ARGMATCH("-t",2) || ARGMATCH ("--threads",2))
      {
#ifdef OMP
	numThreads = atoi(argv[-1]) ;
	if (numThreads > omp_get_max_threads ()) numThreads = omp_get_max_threads () ;
	omp_set_num_threads (numThreads) ;
	fprintf (outFile, "  setting number of threads to %d\n", numThreads) ;
#else
	fprintf (stderr, "  can't set thread number - not compiled with OMP\n") ;
#endif
      }
    else if (ARGMATCH("-o",2) || ARGMATCH("--output",2))
      { if (!strcmp (argv[-1], "-"))
	  outFile = stdout ;
	else if (!(outFile = fopen (argv[-1], "w")))
	  { fprintf (stderr, "can't open output file %s\n", argv[-1]) ; outFile = stdout ; }
      }
    else if (ARGMATCH("--interactive",1)) isInteractive = TRUE ;
    else if (ARGMATCH("--tables",1)) isPrintTables = !isPrintTables ;
    else if (ARGMATCH("--verbose",1)) isVerbose = !isVerbose ;
    else if (ARGMATCH("--readFQB",2))
      { if (!(f = fopen (argv[-1], "r"))) die ("failed to open fqb file %s", argv[-1]) ;
	initialise (params.k, params.w, params.r, params.B) ;
        readFQB (f, params.N, params.chunkSize) ;
	fillHashTable () ;
      }
    else if (ARGMATCH("--readHash",2))
      { if (!(f = fopen (argv[-1], "r"))) die ("failed to open hash file %s", argv[-1]) ;
	initialise (params.k, params.w, params.r, params.B) ;
	readHashFile (f) ;
	fillHashTable () ;
      }
    else if (ARGMATCH("--writeHash",2))
      { if (!(f = fopen (argv[-1], "w"))) die ("failed to open hash file %s", argv[-1]) ;
	writeHashFile (f) ;
      }
    else if (ARGMATCH("--hashCountRange",3))
      hashWithinRangeBuild (atoi(argv[-2]), atoi(argv[-1])) ;
    else if (ARGMATCH("--hashInfo",4))
      hashInfo (atoi(argv[-3]), atoi(argv[-2]), atoi(argv[-1])) ;
    else if (ARGMATCH("--hashExplore",2))
      exploreHash (atoi(argv[-1])) ;
    else if (ARGMATCH("--doubleShared",3))
      doubleShared (atoi(argv[-2]), atoi(argv[-1])) ;
    else if (ARGMATCH("--cribBuild",3))
      { FILE *f1, *f2 ;
	if (!(f1 = fopen (argv[-2], "r"))) die ("failed to open .fa file %s", argv[-2]) ;
	if (!(f2 = fopen (argv[-1], "r"))) die ("failed to open .fa file %s", argv[-1]) ;
	cribBuild (f1, f2) ;
      }
    else if (ARGMATCH("-ct",2) || ARGMATCH ("--clusterThreshold",2))
      params.clusterThreshold = atoi(argv[-1]) ;
    else if (ARGMATCH("--cluster",3))
      { int code, codeMin = atoi(argv[-2]), codeMax = atoi(argv[-1]) ;
	if (!codeMin) codeMin = 1 ;
	if (!codeMax) codeMax = arrayMax(barcodeBlocks) ;
	if (!goodHashes) goodHashesBuild () ;
#ifdef OMP
#pragma omp parallel for
#endif
	for (code = codeMin ; code < codeMax ; ++code)
	  { codeClusterFind (code, params.clusterThreshold) ;
	    codeClusterReadMerge (code) ;
	  }
	fprintf (outFile, "  clustered codes %d to %d\n", codeMin, codeMax) ;
	if (outFile != stdout)
	  printf ("  clustered codes %d to %d\n", codeMin, codeMax) ;
	timeUpdate (stdout) ;
      }
    else if (ARGMATCH("--clusterReport",3))
      { int code, codeMin = atoi(argv[-2]), codeMax = atoi(argv[-1]) ;
	if (!codeMax) codeMax = arrayMax(barcodeBlocks) ;
	codeClusterReport (codeMin, codeMax) ;
      }
    else if (ARGMATCH("--clusterSplit",1))
      clusterSplitCodes () ;
    else if (ARGMATCH("--hashStats",1))
      { if (!(f = fopen (argv[-1], "w"))) die ("failed to open hash stats file %s", argv[-1]) ;
	hashCountHist () ;
      }
    else if (ARGMATCH("--codeStats",1))
      { if (!(f = fopen (argv[-1], "w"))) die ("failed to open code stats file %s", argv[-1]) ;
	codeSizeHist () ;
      }
    else if (ARGMATCH("--errorFix",3))
      errorFix (atoi(argv[-2]), atoi(argv[-1])) ;
    else if (ARGMATCH("--shareScan",3))
      shareScan (atoi(argv[-2]), atoi(argv[-1])) ;
    else if (ARGMATCH("--help",1)) usage () ;
    else if (ARGMATCH("--quit",1) || ARGMATCH("--exit",1)) break ; /* end loop */
    else
      if (isInteractive) fprintf (stderr, "  unknown option/command %s", (*argv)+2) ;
      else die ("unknown option/command %s; run without arguments for usage", *argv) ;

    printf ("  ") ; timeUpdate (stdout) ; fflush (stdout) ;
    
    if (isInteractive)		/* get an input line from stdin and convert to *argv */
      { if (argc) fprintf (stderr, "INPUT ERROR - ignoring arguments starting %s\n", *argv) ;
	char c ;
	printf ("> ") ; fflush (stdout) ;
	while ((c = getchar())) if (!isspace (c) || c == '\n') break ; /* remove leading space */
	if (c == '\n') { argc = 1 ; argv = inputArgs ; *argv = "--help" ; continue ; }
	arr(inputArray,0,char) = '-' ; arr(inputArray,1,char) = '-' ; arrayMax(inputArray) = 2 ;
	while (c != '\n' && c != EOF)
	  { for ( ; isgraph(c) ; c = getchar ())
	      array(inputArray, arrayMax(inputArray), char) = c ;
	    array(inputArray, arrayMax(inputArray), char) = 0 ; /* string terminator */
	    while (isspace (c) && c != '\n') c = getchar() ;	/* eat whitespace separator */
	  }
	if (c == EOF) break ;	       /* end of file or error - exit program */
	argv = inputArgs ; argc = 0 ; int i = 0 ;
	while (i < arrayMax(inputArray))
	  { argv[argc++] = arrp(inputArray,i,char) ;
	    while (*arrp(inputArray,i++,char)) ;
	  }
      }
  }

  fprintf (outFile, "total resources used: ") ; timeTotal (outFile) ;
  if (outFile != stdout) { printf ("total resources used: ") ; timeTotal (stdout) ; }
}

/*************** end of file ***************/

#include "hash.h"

void exploreHash2 (int x, int clusterThreshold)
{
  Array hash1 = hashNeighbours (x) ; if (!hash1) return ;
  fprintf (outFile, "  %d hashes sharing codes with %s\n", arrayMax(hash1), cribText (x)) ;

  Array code1Count = arrayCreate (1024, U32) ;
  HASH code1Hash = hashCreate (1024) ;
  int i, j ;
  for (i = 0 ; i < arrayMax(hash1) ; ++i)
    if (arrp(hash1,i,CodeHash)->read >= clusterThreshold)
      { for (j = 0 ; j < 27 ; ++j) ; }
  
  arrayDestroy (hash1) ; arrayDestroy (code1Count) ; hashDestroy (hash1) ;
}

/************

basic situation

    hetA    -   hapA   -   hom    -  hapB   -   hetB
    hashes      codes     hashes     codes      hashes

Without loss of generality, if we start from a single putative het, then is it hetA.
All code1 codes on this are hapA codes.
Then neigh2 hashes are either hetA or hom, or elsewhere in the genome or errs - require > 3 hits.
Now from hets you get back the hapA codes (all of them this time, not just the code1 subset).

**************/
