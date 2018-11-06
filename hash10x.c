/*  File: hash10x.c
 *  Author: Richard Durbin (rd109@gen.cam.ac.uk)
 *  Copyright (C) Richard Durbin, 2018
 *-------------------------------------------------------------------
 * Description: seqhash-based analysis of 10X Genomics linked read data sets
    reads fqb binary output, and writes as custom binary, by convention .hash 
    look to identify read clusters, and het kmers
    can evaluate with crib built from external fasta files
 * Exported functions:
 * HISTORY:
 * Last edited: Nov  6 17:16 2018 (rd109)
 * Created: Mon Mar  5 11:38:47 2018 (rd)
 *-------------------------------------------------------------------
 */

// choose seed so that poly-A and poly-C don't make moshes

#include "seqhash.h"
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
  U8 subCluster ;	        /* 0 means not in a cluster, else 1..nSubCluster */
  U8 isHet:1 ;			/* is heterozygous - single copy */
  U8 isHom:1 ;			/* is homozygous - two homologous copies */
  U8 isErr:1 ;			/* is an error */
  U8 isSure:1 ;			/* confident about Het/Hom/Err status */
} ClusterHash ;

typedef struct {
  U32 parentCluster ;
  U8  parentSubCluster ;
  U8  isCis ;			/* on same haplotype as parentCluster */
  U8  isTrans ;			/* on alternative haplotype: unknown if neither cis nor trans */
} ClusterInfo ;

int compareClusterHash (const void *a, const void *b)
{ U64 ua = ((ClusterHash*)a)->hash, ub = ((ClusterHash*)b)->hash ;
  if (ua < ub) return -1 ; else if (ua == ub) return 0 ; else return 1 ;
}

int compareClusterRead (const void *a, const void *b)
{ U16 ra = ((ClusterHash*)a)->read, rb = ((ClusterHash*)b)->read ;
  if (ra < rb) return -1 ; else if (ra == rb) return 0 ; else return 1 ;
}

typedef struct {
  U32 nRead ;			/* actually number of read pairs */
  U32 nHash ;			/* number of unique hashes */
  U32 nSubCluster ;
  U32 clusterParent ; 		/* 1 + offset; 0 if these are original barcodes */
  ClusterHash *clusHash ;	/* the unique hashes for this cluster - sorted on hashN  */
  //  Array cluster ;		/* of ClusterInfo */
  double pointToMin ;		/* number of good hashes that point to the minimum hash for the cluster - as a fraction of nGoodHashes an indication of hash and so read density */
} ClusterBlock ;

typedef struct {
  int n ;			/* number of clusters in block */
  ClusterBlock *blocks ;	/* information per cluster block - size n */
  U32  *hCount ;		/* count of each hash in this cluster - size hashNumber */
  U32* *hClusters ;		/* clusters that each  */
} ClusterSet ;

/* strategy for splitting out clusters into their own barcodes is to leave the original 
   barcode entries in place containing any unclustered hashes, then to follow with new
   clusters, grouped in order of their clusterParent */

/* globals */

Seqhash *hasher ;

int hashTableBits = 28 ;	/* default 28, max 34 so hashNumber < 2^32 so index is 32bit */
U64 hashTableSize ; 		/* will be 1 << hashTableBits */
U64 hashTableMask ;		/*  = hashTableSize - 1 */
U32 hashNumber ;		/* number of entries; U32 since less than (hashTableSize >> 2) */
U32 *hashIndex ;		/* index to be used into arrays - hashTableSize */
U64 *hashValue ;		/* value indexed by index - hashTableSize >> 2 */

Array hashDepth ;		/* of U32 */
Array hashCodes ;		/* of U32*, with hashDepth[i] U32 entries */
Array clusterBlocks ;		/* one ClusterBlock per barcode */

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
    for (i = 10 ; i < 15 ; ++i) for (j = 32 ; j-- ; ) *s++ = (u[i] >> j) & 1 ;
  if ((s = s2))
    for (i = 15 ; i < 25 ; ++i) for (j = 16 ; j-- ; ) *s++ = (u[i] >> (2*j)) & 3 ;
  if ((s = q2))
    for (i = 25 ; i < 30 ; ++i) for (j = 32 ; j-- ; ) *s++ = (u[i] >> j) & 1 ;
}

typedef struct { U64 hash ; int read, pos ; } HashTemp ; /* just for processing */

void seqAddHashes (char *s, int len, Array a, int readCount)
{
  SeqhashRCiterator *mi = moshRCiterator (hasher, s, len) ;
  U64 hash ; int pos ;
  while (moshRCnext (mi, &hash, &pos))
    { HashTemp *h = arrayp(a,arrayMax(a),HashTemp) ;
      h->hash = hash ; h->read = readCount ; h->pos = pos ;
    }
  seqhashRCiteratorDestroy (mi) ;
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
  while ((index = hashIndex[offset]) && (hashValue[index] != hash))
    offset = (offset + diff) & hashTableMask ;
  if (!index && isAdd)
    { index = hashIndex[offset] = hashNumber++ ;
      hashValue[index] = hash ;
      if ((++nIndex << 2) > hashTableSize) die ("hashTableSize is too small") ;
      // Alex changed this to be ((++nIndex << 1) > hashTableSize) (and nIndex to be unsigned 32)
    }
  return index ;
}

void processBlock (U32 *readData, ClusterBlock *b)
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
  b->clusHash = new (n, ClusterHash) ;
  for (i = 0 ; i < n ; ++i)
    { int index = hashIndexFind (arrp(a,i,HashTemp)->hash, TRUE) ;
      ++array(hashDepth, index, U32) ;
      b->clusHash[i].hash = index ;
      b->clusHash[i].read = arrp(a,i,HashTemp)->read ;
    }

  qsort (b->clusHash, n, sizeof(ClusterHash), compareClusterHash) ; /* sort so hashes are in order */
  
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
  ClusterBlock *b = arrayp(clusterBlocks, 1, ClusterBlock) ; /* start at 1 */
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
	    b = arrayp(clusterBlocks, arrayMax(clusterBlocks), ClusterBlock) ; b->nRead = 1 ;
	    barcode = u[30*i] ; uStart = &u[30*i] ;
	  }

      nReads += nRec ;
    }
  fclose (f) ;

  fprintf (outFile, "  read %d read pair records for %d barcodes, mean %.2f read pairs per barcode\n",
	  nReads, arrayMax(clusterBlocks)-1, nReads/(double)(arrayMax(clusterBlocks)-1)) ;
  fprintf (outFile, "  created %ld hashes, mean %.2f hashes per read pair, %.2f per barcode\n",
	  (long)nHashes, nHashes/(double)nReads, nHashes/(double)(arrayMax(clusterBlocks)-1)) ;
  if (outFile != stdin)
    { printf ("  read %d read pair records for %d barcodes, mean %.2f read pairs per barcode\n",
	  nReads, arrayMax(clusterBlocks)-1, nReads/(double)(arrayMax(clusterBlocks)-1)) ;
      printf ("  created %ld hashes, mean %.2f hashes per read pair, %.2f per barcode\n",
	  (long)nHashes, nHashes/(double)nReads, nHashes/(double)(arrayMax(clusterBlocks)-1)) ;
    }
}

/*****************************************************************/

static U32 HASHFILE_VERSION = 2 ;
static U16 CLUSTERHASH_SIZE = sizeof(ClusterHash) ;
static U16 CLUSTERBLOCK_SIZE = sizeof(ClusterBlock) ;

void writeHashFile (FILE *f)
{
  if ((fwrite ("10XH", 4, 1, f) != 1) || (fwrite (&HASHFILE_VERSION, 4, 1, f) != 1) ||
      (fwrite (&CLUSTERHASH_SIZE, 2, 1, f) != 1) || (fwrite (&CLUSTERBLOCK_SIZE, 2, 1, f) != 1) ||
      (fwrite (&hashTableBits, 4, 1, f) != 1)) die ("write fail 1") ;
  if (fwrite (hashIndex, sizeof(U32), hashTableSize, f) != hashTableSize) die ("write fail 2") ;
  if (fwrite (&hashNumber, sizeof(U32), 1, f) != 1) die ("failed to write hashNumber") ;
  if (fwrite (hashValue, sizeof(U64), hashNumber, f)) die ("failed to write hashValue") ;
  
  if (!arrayWrite (hashDepth, f)) die ("failed to write hashDepth array") ;
  if (!arrayWrite (clusterBlocks, f)) die ("failed to write clusterBlocks array") ;
  int i ;
  for (i = 1 ; i < arrayMax(clusterBlocks) ; ++i)
    { ClusterBlock *b = arrp(clusterBlocks, i, ClusterBlock) ;
      if (fwrite (b->clusHash, sizeof(ClusterHash), b->nHash, f) != b->nHash) die ("write fail 3") ;
    }
  fclose (f) ;

  fprintf (outFile, "  wrote %lld hash table entries and %d barcode blocks\n",
	   hashTableSize, arrayMax(clusterBlocks)) ;
  if (outFile != stdout)
    printf ("  wrote %lld hash table entries and %d barcode blocks\n",
	    hashTableSize, arrayMax(clusterBlocks)) ;
}

void readHashFile (FILE *f)
{
  int B ; U32 version ; U16 clusterHashSize, clusterBlockSize ;
  char name[5] = "abcd" ;
  if ((fread (name, 4, 1, f) != 1) || (fread (&version, 4, 1, f) != 1) ||
      (fread (&clusterHashSize, 2, 1, f) != 1) || (fread (&clusterBlockSize, 2, 1, f) != 1))
    die ("read fail 0") ;
  if (strcmp (name, "10XH")) die ("not a 10X hash file") ;
  if (version > HASHFILE_VERSION)
    die ("hash file version mismatch: file %d > code %d", version, HASHFILE_VERSION) ;
  if (clusterHashSize != CLUSTERHASH_SIZE)
    die ("ClusterHash structure size mismatch: file %d != code %d", clusterHashSize, CLUSTERHASH_SIZE) ;
  if (clusterBlockSize != CLUSTERBLOCK_SIZE)
    die ("ClusterBlock structure size mismatch: file %d != code %d", clusterBlockSize, CLUSTERBLOCK_SIZE) ;
  if (fread (&B, 4, 1, f) != 1) die ("read fail 1") ;
  if (B != hashTableBits) die ("incompatible hash table size: rerun with -B %d", B) ;
  if (fread (hashIndex, sizeof(U32), hashTableSize, f) != hashTableSize) die ("read fail 2") ;
  if (version == 1)
    { Array a ;
      if (!(a = arrayRead (f))) die ("failed to read hashValue array") ;
      hashNumber = arrayMax(a) ; hashValue = new(hashNumber, U64) ;
      memcpy (hashValue, arrp(a,0,U64), hashNumber*sizeof(U64)) ;
      arrayDestroy (a) ;
    }
  else if (version == 2)
    { if (fread (&hashNumber, sizeof(U32), 1, f) != 1) die ("failed to read hashNumber") ;
      if (fread (hashValue, sizeof(U64), hashNumber, f)) die ("failed to read hashValue") ;
    }

  if (!(hashDepth = arrayRead (f))) die ("failed to read hashDepth array") ;
  if (!(clusterBlocks = arrayRead (f))) die ("failed to read clusterBlocks array") ;
  int i ;
  U64 nReads = 0, nHashes = 0 ;
  for (i = 1 ; i < arrayMax(clusterBlocks) ; ++i)
    { ClusterBlock *b = arrayp(clusterBlocks, i, ClusterBlock) ;
      nReads += b->nRead ; nHashes += b->nHash ;
      b->clusHash = new(b->nHash, ClusterHash) ;
      if (fread (b->clusHash, sizeof(ClusterHash), b->nHash, f) != b->nHash) die ("read fail 3") ;
    }
  fclose (f) ;

  fprintf (outFile, "  read %ld hashes for %ld reads in %d barcode blocks\n",
	  (long)nHashes, (long)nReads, arrayMax(clusterBlocks)) ;
  if (outFile != stdout)
    printf ("  read %ld hashes for %ld reads in %d barcode blocks\n",
	    (long)nHashes, (long)nReads, arrayMax(clusterBlocks)) ;
}

void fillHashTable (void)
{
  int i, j ;
  long nHashes = 0 ;

  hashCodes = arrayCreate(hashNumber,U32*) ;
  arrayMax(hashCodes) = hashNumber ;
  
  U32 *count = arrp(hashDepth,1,U32) ;
  for (i = 1 ; i < hashNumber ; ++i, ++count)
    if (*count)
      { arr(hashCodes,i,U32*) = new(*count, U32) ;
	nHashes += *count ;
      }

  U32 *hashN = new0 (hashNumber, U32) ;
  for (i = 1 ; i < arrayMax(clusterBlocks) ; ++i)
    { ClusterBlock *b = arrp(clusterBlocks, i, ClusterBlock) ;
      ClusterHash *c = b->clusHash ;
      for (j = 0 ; j < b->nHash ; ++j, ++c)
	arr(hashCodes,c->hash,U32*)[hashN[c->hash]++] = i ;
    }
  fprintf (outFile, "  filled hash table: %ld hashes from %d barcodes in %d bins\n",
	  nHashes, arrayMax(clusterBlocks), hashNumber) ;
#ifdef CHECK
  for (i = 0 ; i < hashNumber ; ++i)
    if (hashN[i] != arr(hashDepth,i,U32))
      die ("hashN[%d] = %d != hashDepth[%d] = %d", i, hashN[i], i, arr(hashDepth,i,U32)) ;
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

void hashDepthHist (void)
{
  if (!arrayMax(hashDepth)) { fprintf (stderr, "  no hash list to print stats for\n") ; return ; }
  Array a = arrayCreate (1024, int) ;
  int i ;
  for (i = 0 ; i < arrayMax(hashDepth) ; ++i)
    ++array(a, arr(hashDepth, i, U32), int) ;
  histogramReport ("HASH_COUNT", a) ;
  arrayDestroy (a) ;
}

void codeSizeHist (void)
{
  if (!clusterBlocks || !arrayMax (clusterBlocks))
    { fprintf (stderr, "  no barcodes to print stats for\n") ; return ; }
  Array hashHist = arrayCreate (1024, int) ;
  Array clusterHist = arrayCreate (1024, int) ;
  int i ;
  for (i = 0 ; i < arrayMax(clusterBlocks) ; ++i)
    { ClusterBlock *b = arrp (clusterBlocks, i, ClusterBlock) ;
      ++array(hashHist, b->nHash, int) ;
      ++array(clusterHist, b->nSubCluster, int) ;
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
    { SeqhashRCiterator *mi = moshRCiterator (hasher, seq, len) ;
      U64 hash ; int pos ;
      ++chr ; if (chr >= cribChrMax) cribChrMax = chr+1 ;
      while (moshRCnext (mi, &hash, &pos))
	if ((index = hashIndexFind (hash, FALSE)))
	  { CribInfo *c = &crib[index] ;
	    if (!c->chr) { c->chr = chr ; c->pos = pos >> 10 ; }
	    else if (c->chr > 0) c->chr = -1 ;
	    else --c->chr ;
	    ++nPresent ;
	  }
	else
	  ++nAbsent ;
      seqhashRCiteratorDestroy (mi) ;
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
  crib = new0 (hashNumber, CribInfo) ; cribAddGenome (f1, crib) ;
  CribInfo *crib2 = new0 (hashNumber, CribInfo) ; cribAddGenome (f2, crib2) ;

  cribType = new0 (hashNumber, U8) ;
  Array aHom = arrayCreate (256, int), aHet = arrayCreate (256, int) ;
  Array aErr = arrayCreate (256, int), aMul = arrayCreate (256, int) ;
  int i ;
  for (i = 1 ; i < hashNumber ; ++i)
    if (crib[i].chr == 0 && crib2[i].chr == 0)
      { cribType[i] = CRIB_ERR ; ++array(aErr, arr(hashDepth, i, int), int) ; }
    else if (crib[i].chr > 0 && crib2[i].chr > 0)
      { cribType[i] = CRIB_HOM ; ++array(aHom, arr(hashDepth, i, int), int) ; }
    else if (crib[i].chr < 0 || crib2[i].chr < 0)
      { cribType[i] = CRIB_MUL ; ++array(aMul, arr(hashDepth, i, int), int) ;
	if (crib2[i].chr < crib[i].chr) crib[i].chr = crib2[i].chr ;
      }
    else
      { cribType[i] = crib[i].chr ? CRIB_HTA : CRIB_HTB ;
	++array(aHet, arr(hashDepth, i, int), int) ;
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
  sprintf (s, "-%d", arr(hashDepth, x, U32)) ;
  return text ;
}

/************************************************************/

static BOOL *hashWithinRange = 0 ;
static int hashRangeMin = 0, hashRangeMax = 0 ;

BOOL hashWithinRangeBuild (int min, int max)
{ 
  if (hashWithinRange && min == hashRangeMin && max == hashRangeMax) return FALSE ;
  if (!hashWithinRange) hashWithinRange = new0 (hashNumber, BOOL) ;
  int i, n ;
  for (i = 0 ; i < hashNumber ; ++i)
    { n = arr(hashDepth, i, U32) ;
      if (n >= min && n < max) hashWithinRange[i] = TRUE ;
    }
  hashRangeMin = min ; hashRangeMax = max ;
  return TRUE ;
}

Array hashNeighbours (int x)
{ /* return array of hashes that share between min and max codes with x */
  /* abuse ClusterHash structure for this array, storing the number of shared codes in ->read */
  if (!hashWithinRange) die ("hashNeighbours() called without setting hashDepthRange") ;
  Array ch = arrayCreate (1000, ClusterHash) ;
  int i,j ;
  U32 xCount = arr(hashDepth,x,U32) ; if (!xCount) return 0 ;
  for (i = 0 ; i < xCount ; ++i)
    { ClusterBlock *b = arrp(clusterBlocks, arr(hashCodes,x,U32*)[i], ClusterBlock) ;
      ClusterHash *c = b->clusHash ;
      for (j = 0 ; j < b->nHash ; ++j, ++c)
	if (c->hash != x && hashWithinRange[c->hash])
	  array(ch, arrayMax(ch), ClusterHash) = *c ;
    }
  arraySort (ch, compareClusterHash) ;
  Array shared = arrayCreate (256, ClusterHash) ; /* abuse structure esp. read field */
  ClusterHash *last = arrayp (shared, 0, ClusterHash) ;
  last->hash = arrp(ch, 0, ClusterHash)->hash ; last->read = 1 ;
  for (i = 1 ; i < arrayMax(ch) ; ++i)
    if (arrp(ch, i, ClusterHash)->hash == last->hash) ++last->read ;
    else
      { last = arrayp (shared, arrayMax(shared), ClusterHash) ;
	last->hash = arrp(ch, i, ClusterHash)->hash ; last->read = 1 ;
      }
  arrayDestroy (ch) ;
  return shared ;
}

Array countHashNeighbours (int x)
{
  U32 xCount = arr(hashDepth,x,U32) ; if (!xCount) return 0 ;
  HASH h = hashCreate (500000) ;
  Array a = arrayCreate (500000, int) ;
  int i,j ;
  for (i = 0 ; i < xCount ; ++i)
    { ClusterBlock *b = arrp(clusterBlocks, arr(hashCodes,x,U32*)[i], ClusterBlock) ;
      ClusterHash *c = b->clusHash ;
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
    { fprintf (stderr, "exploreHash called without hashDepthRange\n") ; return ; }
  Array shared = hashNeighbours (x) ; if (!shared) return ;
  fprintf (outFile, "  %d hashes sharing codes with %s\n", arrayMax(shared), cribText(x)) ;
  arraySort (shared, compareClusterRead) ;
  int i, count = 1, n = 0 ;
  ClusterHash *d = arrp(shared, 0, ClusterHash) ;
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
    { fprintf (stderr, "doubleShared called without hashDepthRange\n") ; return ; }
  Array sh1 = hashNeighbours (x1) ; if (!sh1) return ;
  Array sh2 = hashNeighbours (x2) ; if (!sh2) { arrayDestroy (sh1) ; return ; }

  int i1 = 0, i2 = 0 ;
  ClusterHash *c1 = arrp(sh1, 0, ClusterHash), *c2 = arrp(sh2, 0, ClusterHash) ;
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
    { fprintf (stderr, "hashInfo called without hashDepthRange\n") ; return ; }
  int x ;
  for (x = hMin ; x < hMax ; x += hIncrement)
    { fprintf (outFile, "HASH_INFO  %s", cribText(x)) ;
      if (hashWithinRange[x])
	{ Array shared = hashNeighbours (x) ;
	  arraySort (shared, compareClusterRead) ;
	  ClusterHash *ch = arrp(shared, arrayMax(shared)-1, ClusterHash) ;
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
      int xCount = arr(hashDepth, x, U32) ; if (xCount < 5) { ++low[type] ; continue ; }
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
  for (x = countMax ; x < hashNumber ; x += 100)
    { int xCount = arr(hashDepth, x, U32) ;
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
}

/***************************************/

U16 **goodHashes = 0 ;
int *nGoodHashes = 0 ;
static ClusterHash *clusHashBase ;

int compareGoodHashes (const void *a, const void *b) /* sort by increasing hashDepth */
{ U32 ha = clusHashBase[*(U16*)a].hash, hb = clusHashBase[*(U16*)b].hash ;
  U32 ca = arr(hashDepth, ha, U32), cb = arr(hashDepth, hb, U32) ;
  if (ca < cb) return -1 ; else if (ca == cb) return 0 ; else return 1 ;
}

/* As currently written this is not threadsafe, because we set clusHashBase to access in the sort.
   We could create a data structure with the index and the value to sort on, and then build
   a pre-allocated array of size b->nHash of that and sort it, and then copy into goodHashes.
   Then you could parallelise this.
*/

void goodHashesBuild (void)
{
  if (!hashWithinRange) die ("cluster code called without setting hashDepthRange") ;
  goodHashes = new (arrayMax(clusterBlocks), U16*) ;
  nGoodHashes = new (arrayMax(clusterBlocks), int) ;
  Array a = arrayCreate (4096, U16) ;
  int c, i ;
  
  for (c = 0 ; c < arrayMax(clusterBlocks) ; ++c)
    { ClusterBlock *cb = arrp(clusterBlocks, c, ClusterBlock) ;
      if (cb->nHash > U16MAX)
	{ nGoodHashes[c] = 0 ; goodHashes[c] = new(1,U16) ;
	  fprintf (stderr, "ignoring barcode %d - too many hashes %d > %d\n",
		   c, cb->nHash, U16MAX) ;
	  continue ;
	}
      ClusterHash *ch = clusHashBase = cb->clusHash ;
      arrayMax(a) = 0 ;
      for (i = 0 ; i < cb->nHash ; ++i, ++ch)
	if (hashWithinRange[ch->hash]) array(a, arrayMax(a), U16) = i ;
      arraySort (a, compareGoodHashes) ; /* sort by increasing hashDepth */
      goodHashes[c] = new(arrayMax(a),U16) ;
      memcpy (goodHashes[c], arrp(a,0,U16), arrayMax(a)*sizeof(U16)) ;
      nGoodHashes[c] = arrayMax(a) ;
    }

  printf ("  made goodHashes arrays for hash range %d to %d\n  ", hashRangeMin, hashRangeMax) ;
  timeUpdate (outFile) ; fflush (outFile) ;
}

/*************************/

void codeClusterFind (int code, int clusterThreshold) /* clusters this barcode's good hashes */
{
  ClusterBlock *b = arrp(clusterBlocks, code, ClusterBlock), *c ;
  ClusterHash *cb = b->clusHash, *cc ;
  int i, j ;
  U32 x, nc ;
  int *minShare = new0 (arrayMax(clusterBlocks), int) ; /* lowest hash shared with code */
  U16 *g, *gi, *gc, *gcn ;
  Array clusterMin = arrayCreate (64, int) ;

  int n = nGoodHashes[code] ; if (!n) return ;
  int *minShareCount = new (n, int) ;

  for (i = 0 ; i < n ; ++i) b->clusHash[goodHashes[code][i]].subCluster = 0 ; /* wipe existing */

  b->nSubCluster = 0 ;
  //  if (!b->cluster) b->cluster = arrayCreate (32, ClusterInfo) ;
  b->pointToMin = 0.0 ;
  int nSubClustered = 0 ;
  for (i = 1 ; i < n ; ++i)	// shouldn't this be from i = 0 ?
    { x = b->clusHash[goodHashes[code][i]].hash ;
      nc = arr(hashDepth, x, U32) ;
      memset (minShareCount, 0, n*sizeof(int)) ;
					       
      for (j = 0 ; j < nc ; ++j)
	{ int cj = arr(hashCodes, x, U32*)[j] ;
	  if (cj == code) continue ;
	  if (!minShare[cj]) minShare[cj] = i+1 ;
	  ++minShareCount[minShare[cj]-1] ;
	}

      int msMax = 0, msTot = 0, msBest ;
      g = goodHashes[code] ;
      for (j = 0 ; j < i ; ++j)
	{ if (minShareCount[j] > msMax)  { msBest = j ; msMax = minShareCount[j] ; }
	  msTot += minShareCount[j] ;
	}
      if (msMax >= clusterThreshold)
	{ ClusterHash *ch = &b->clusHash[g[msBest]] ;
	  if (!ch->subCluster)	/* create a new cluster */
	    { if (++b->nSubCluster > U8MAX) /* abandon this clustering */
		{ b->nSubCluster = 0 ; nSubClustered = 0 ;
		  for (j = 0 ; j < i ; ++j) b->clusHash[g[j]].subCluster = 0 ;
		  fprintf (stderr, "    code %d with %d good hashes has too many clusters\n",
			   code, nGoodHashes[code]) ;
		  break ;
		}
	      ch->subCluster = b->nSubCluster ; ++nSubClustered ;
	      array(clusterMin,b->nSubCluster,int) = msBest ;
	    }
	  b->clusHash[g[i]].subCluster = ch->subCluster ; ++nSubClustered ;
	  b->pointToMin += minShareCount[arr(clusterMin,ch->subCluster,int)] / (double)msTot ;
	}
    }
  free (minShare) ;
  free (minShareCount) ;
  arrayDestroy (clusterMin) ;
  if (isVerbose)
    { if (b->nSubCluster)
	fprintf (outFile, "  code %d with %d reads %d hashes, %d good hashes, of which %d cluster into %d raw",
		 code, b->nRead, b->nHash, nGoodHashes[code], nSubClustered, b->nSubCluster) ;
      else
	fprintf (outFile, "  code %d with %d reads %d hashes, %d good hashes, 0 clusters\n",
		 code, b->nRead, b->nHash, nGoodHashes[code]) ;
    }
}

void codeClusterReadMerge (int code)
{
  ClusterBlock *b = arrp(clusterBlocks, code, ClusterBlock) ;
  if (!b->nSubCluster) return ;	/* nothing to do */
  int i, j ;
  int *readMap = new0 (b->nRead, int), *trueCluster = new (b->nSubCluster+1, int) ;
  int *deadCluster = new0 (b->nSubCluster+1, int) ;
  for (i = 0 ; i <= b->nSubCluster ; ++i) trueCluster[i] = i ; /* initialisation */
  for (i = 0 ; i < b->nHash ; ++i)
    { int hashCluster = trueCluster[b->clusHash[i].subCluster] ; if (!hashCluster) continue ;
      int readCluster = trueCluster[readMap[b->clusHash[i].read]] ;
      if (hashCluster == readCluster) continue ;
      else if (!readCluster)
	readMap[b->clusHash[i].read] = hashCluster ;
      else			/* merge clusters hashCluster and readCluster */
	{ if (hashCluster > readCluster)
	    { int t = hashCluster ; hashCluster = readCluster ; readCluster = t ; }
	  for (j = 1 ; j <= b->nSubCluster ; ++j)
	    if (trueCluster[j] == readCluster) trueCluster[j] = hashCluster ;
	  deadCluster[readCluster] = 1 ;
	}
    }

  /* now compactify trueCluster, reusing deadCluster[] for the map removing dead indices */
  for (j = 1 ; j <= b->nSubCluster ; ++j) deadCluster[j] = deadCluster[j-1] + 1 - deadCluster[j] ;
  for (j = 1 ; j <= b->nSubCluster ; ++j) trueCluster[j] = deadCluster[trueCluster[j]] ;
  b->nSubCluster = deadCluster[b->nSubCluster] ;
  for (i = 0 ; i < b->nHash ; ++i) b->clusHash[i].subCluster = trueCluster[b->clusHash[i].subCluster] ;

  free (readMap) ; free(trueCluster) ; free(deadCluster) ;
  if (isVerbose) { fprintf (outFile, " then %d merged clusters\n", b->nSubCluster) ; fflush (stdout) ; }
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
    { ClusterBlock *b = arrp(clusterBlocks, code, ClusterBlock) ;
      arrayReCreate (info, b->nSubCluster+1, ReportInfo) ;
      arrayReCreate (readClus, b->nRead, int) ;
      arrayReCreate (badLink, b->nHash, int) ;
      int nClusHash = 0 ;
      for (i = 0 ; i < b->nHash ; ++i)
	{ if (!(c = b->clusHash[i].subCluster)) continue ;
	  ++nClusHash ;
	  arr(readClus, b->clusHash[i].read, int) = c ;
	  r = arrp(info, c, ReportInfo) ;
	  ++r->n ;
	  U32 bh = b->clusHash[i].hash ;
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
      
      fprintf (outFile, "  CLUSTER_SUMMARY %d nRead %d nHash %d nGoodHash %d nClusHash %d nClusRead %d nSubCluster %d\n",
	      code, b->nRead, b->nHash, nGoodHashes[code], nClusHash, nClusRead, b->nSubCluster) ;
      for (i = 1 ; i <= b->nSubCluster ; ++i)
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
		    { fprintf (outFile, " %s", cribText(b->clusHash[x].hash)) ;
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
  int i, j, nSubClusters = 0 ;
  int nCodes = arrayMax (clusterBlocks) ; 
  for (i = 0 ; i < nCodes ; ++i) nSubClusters += arrp(clusterBlocks, i, ClusterBlock)->nSubCluster ;
  Array newBlocks = arrayCreate (nCodes+nSubClusters, ClusterBlock) ;
  arrayMax(newBlocks) = nCodes + nSubClusters ;
  ClusterBlock *old = arrp(clusterBlocks, 0, ClusterBlock) ;
  ClusterBlock *new1 = arrp(newBlocks, 0, ClusterBlock) ;
  ClusterBlock *new2 = arrp(newBlocks, nCodes-1, ClusterBlock) ; //-1 allows for ->subCluster offset
  Array clusterHashCount = arrayCreate (256, int) ; /* how many hashes in each cluster */
  Array readMap = arrayCreate (4096, int) ;
  for (i = 0 ; i < nCodes ; ++i, ++old)
    if (old->nSubCluster)
      { arrayReCreate (clusterHashCount, old->nSubCluster+1, int) ; /* zeroes array */
	ClusterHash *c = old->clusHash ;
	for (j = 0 ; j < old->nHash ; ++j, ++c) ++arr(clusterHashCount,c->subCluster,int) ;
	new1->clusHash = new0(arr(clusterHashCount,0,int), ClusterHash) ;
	for (j = 1 ; j < old->nSubCluster+1 ; ++j)
	  { new2[j].clusHash = new0(arr(clusterHashCount,j,int), ClusterHash) ;
	    new2[j].clusterParent = i+1 ;
	  }
	arrayReCreate (clusterHashCount, old->nSubCluster+1, int) ;
	arrayReCreate (readMap, old->nRead+1, int) ;
	for (j = 0, c = old->clusHash ; j < old->nHash ; ++j, ++c)
	  { int clus = c->subCluster ; c->subCluster = 0 ;
	    if (clus)
	      { if (!array(readMap,c->read,int)) arr(readMap,c->read,int) = ++new2[clus].nRead ;
	    // previous line requires that all hashes from same read are in same cluster
		c->read = arr(readMap,c->read,int) - 1 ;
		new2[clus].clusHash[new2[clus].nHash++] = *c ;
	      }
	    else new1->clusHash[new1->nHash++] = *c ;
	  }
	new1->nRead = old->nRead ; /* not strictly true any more, but safe */
	++new1 ;
	new2 += old->nSubCluster ;
	free (old->clusHash) ;
      }
    else
      *new1++ = *old ; // also carries over existing clusHash
  
  fprintf (outFile, "  made %d additional new barcodes from clusters in %d original barcodes\n",
	   nSubClusters, arrayMax(clusterBlocks)) ;
  if (outFile != stdout)
    printf ("  made %d additional new barcodes from clusters in %d original barcodes\n",
	    nSubClusters, arrayMax(clusterBlocks)) ;
  arrayDestroy (clusterBlocks) ;
  clusterBlocks = newBlocks ;
  arrayDestroy (clusterHashCount) ; arrayDestroy (readMap) ;
  printf ("  cluster timepoint: ") ; timeUpdate (stdout) ;

  // now must rebuild the mapping from hashes to codes
  for (i = 1 ; i < hashNumber ; ++i)
    if (arr(hashDepth,i,U32)) free (arr(hashCodes,i,U32*)) ;
  arrayDestroy (hashCodes) ;
  fillHashTable () ;
}

/********************/

void cribSummary (void)
{
  if (!crib) { fprintf (stderr, "cribSummary requires crib\n") ; return ; }
  int i, j, nBaseCode = 0, nSubClusterCode = 0 ;
  U64 countBase[CRIBTYPE_MAX], countCluster[CRIBTYPE_MAX] ;
  HASH clusterHashes[CRIBTYPE_MAX], baseHashes[CRIBTYPE_MAX] ;
  for (i = 0 ; i < CRIBTYPE_MAX ; ++i)
    { countBase[i] = countCluster[i] = 0 ;
      clusterHashes[i] = hashCreate (1<<20) ;
      baseHashes[i] = hashCreate (1<<20) ;
    }
  fprintf (stderr, "made hash objects\n") ;
  for (i = 0 ; i < arrayMax(clusterBlocks) ; ++i)
    { ClusterBlock *b = arrp(clusterBlocks, i, ClusterBlock) ;
      if (b->clusterParent)
	{ ++nSubClusterCode ;
	  for (j = 0 ; j < b->nHash ; ++j)
	    { U32 h = b->clusHash[j].hash ;
	      ++countCluster[cribType[h]] ;
	      hashAdd (clusterHashes[cribType[h]], HASH_INT(h)) ;
	    }
	}
      else
	{ ++nBaseCode ;
	  for (j = 0 ; j < b->nHash ; ++j)
	    { U32 h = b->clusHash[j].hash ;
	      ++countBase[cribType[h]] ;
	      hashAdd (baseHashes[cribType[h]], HASH_INT(h)) ;
	    }	      
	}
    }
  fprintf (outFile, "  %d base codes ", nBaseCode) ;
  for (i = 0 ; i < CRIBTYPE_MAX ; ++i)
    fprintf (outFile, " %s %llu %d %.1f",
	     cribTypeName[i], countBase[i], hashCount(baseHashes[i]),
	     countBase[i] / (double)hashCount(baseHashes[i])) ;
  fprintf (outFile, "\n  %d cluster codes ", nSubClusterCode) ;
  for (i = 0 ; i < CRIBTYPE_MAX ; ++i)
    fprintf (outFile, " %s %llu %d %.1f",
	     cribTypeName[i], countCluster[i], hashCount(clusterHashes[i]),
	     countCluster[i] / (double)hashCount(clusterHashes[i])) ;
  fputc ('\n', outFile) ;
for (i = 0 ; i < CRIBTYPE_MAX ; ++i)
    { hashDestroy (baseHashes[i]) ; hashDestroy (clusterHashes[i]) ; }
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
  fprintf (stderr, "   --hashDepthRange <min> <max>: set limits for hash counts for Explore, cluster etc.\n") ;
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
  hasher = seqhashCreate (k, w) ;

  hashTableBits = B ;
  if (hashTableBits < 20 || hashTableBits > 30)
    die ("hashTableBits %d out of range 20-30", hashTableBits) ;
  hashTableSize = (size_t)1 << hashTableBits ;
  hashTableMask = hashTableSize - 1 ;
  hashIndex = new0 (hashTableSize, U32) ;
  hashValue = new0 (hashTableSize>>2, U64) ;
  hashNumber = 1 ; // so true indexes != 0
  hashDepth = arrayCreate (1 << 20, U32) ;

  fprintf (outFile, "hash10x initialised with k = %d, w = %d, random seed = %d, hashtable bits = %d\n",
	  k, w, r, B) ;
}

void codeExplore (int code, int clusterThreshold) ; /* forward declaration - code below */

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
  printf ("sizeof ClusterBlock is %d\n", sizeof(ClusterBlock)) ;
  printf ("sizeof ClusterHash is %d\n", sizeof(ClusterHash)) ;
  exit (0) ;
#endif

  clusterBlocks = arrayCreate (1200, ClusterBlock) ;

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
    else if (ARGMATCH("--hashDepthRange",3))
      { hashWithinRangeBuild (atoi(argv[-2]), atoi(argv[-1])) ;
	goodHashesBuild () ;
      }
    else if (ARGMATCH("--hashInfo",4))
      hashInfo (atoi(argv[-3]), atoi(argv[-2]), atoi(argv[-1])) ;
    else if (ARGMATCH("--hashExplore",2))
      exploreHash (atoi(argv[-1])) ;
    else if (ARGMATCH("--doubleShared",3))
      doubleShared (atoi(argv[-2]), atoi(argv[-1])) ;
    else if (ARGMATCH("--codeExplore",2))
      if (goodHashes)
	codeExplore (atoi(argv[-1]), params.clusterThreshold) ;
      else
	{ fprintf (outFile, "!! you must set hashDepthRange before codeExplore\n") ;
	  if (outFile != stdout) fprintf (stderr, "!! you must set hashDepthRange before codeExplore\n") ;
	}
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
	if (!codeMax) codeMax = arrayMax(clusterBlocks) ;
	if (goodHashes)
	  {
#ifdef OMP
#pragma omp parallel for
#endif
	    for (code = codeMin ; code < codeMax ; ++code)
	      { codeClusterFind (code, params.clusterThreshold) ;
		codeClusterReadMerge (code) ;
	      }
	    fprintf (outFile, "  clustered codes %d to %d\n", codeMin, codeMax) ;
	    if (outFile != stdout) printf ("  clustered codes %d to %d\n", codeMin, codeMax) ;
	  }
      else
	{ fprintf (outFile, "!! you must set hashDepthRange before cluster\n") ;
	  if (outFile != stdout) fprintf (stderr, "!! you must set hashDepthRange before cluster\n") ;
	}
      }
    else if (ARGMATCH("--clusterReport",3))
      { int code, codeMin = atoi(argv[-2]), codeMax = atoi(argv[-1]) ;
	if (!codeMax) codeMax = arrayMax(clusterBlocks) ;
	codeClusterReport (codeMin, codeMax) ;
      }
    else if (ARGMATCH("--clusterSplit",1)) clusterSplitCodes () ;
    else if (ARGMATCH("--hashStats",1)) hashDepthHist () ;
    else if (ARGMATCH("--codeStats",1)) codeSizeHist () ;
    else if (ARGMATCH("--cribSummary",1)) cribSummary () ;
    else if (ARGMATCH("--errorFix",3)) errorFix (atoi(argv[-2]), atoi(argv[-1])) ;
    else if (ARGMATCH("--shareScan",3)) shareScan (atoi(argv[-2]), atoi(argv[-1])) ;
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

/************

basic situation

    hetA    -   hapA   -   hom    -  hapB   -   hetB
    hashes      codes     hashes     codes      hashes

Without loss of generality, if we start from a single putative het, then it is hetA.
All code1 codes on this are hapA codes.
Then neigh2 hashes are either hetA or hom, or elsewhere in the genome or errs - require > 3 hits.
Now from hets you get back the hapA codes (all of them this time, not just the code1 subset).

Alternatively, if we start from a hom, then the neigh2 hashes are either hetA or hetB or hom.
Connect only via clustered codes, and require > clusterThreshold connectivity to remove errs, and we hope mislocalised hashes. 
Single-linkage clustering of het hashes on connection via a shared code should give two clusters.
But a hom would link these clusters.  We could sort the hashes on depth, and process in order of increasing depth.
The neigh2 of a single het is quite big.  Either there will be no overlap, or quite a lot.


**************/

void exploreHash2 (int x, int clusterThreshold)
{
  Array hash1 = hashNeighbours (x) ; if (!hash1) return ;
  fprintf (outFile, "  %d hashes sharing codes with %s\n", arrayMax(hash1), cribText (x)) ;

  Array code1Count = arrayCreate (1024, U32) ;
  HASH code1Hash = hashCreate (1024) ;
  int i, j ;
  for (i = 0 ; i < arrayMax(hash1) ; ++i)
    if (arrp(hash1,i,ClusterHash)->read >= clusterThreshold)
      { for (j = 0 ; j < 27 ; ++j) ; }
  
  arrayDestroy (hash1) ; arrayDestroy (code1Count) ; hashDestroy (hash1) ;
}

/* this takes code from codeClusterFind() and codeClusterMerge() and extends to try to assign 
   het, hom, err and mul status to hashes
*/

typedef struct { int i ; int count ; } CountStruct ;
static int compareCount (const void *a, const void *b)
{ int ca = ((CountStruct*)a)->count, cb = ((CountStruct*)b)->count ; return ca-cb ; }

void codeExplore (int code, int clusterThreshold) /* clusters this barcode's good hashes */
{
  ClusterBlock *b = arrp(clusterBlocks, code, ClusterBlock), *c ;
  ClusterHash *cb = b->clusHash, *cc ;
  int i, j ;
  U32 x, nc ;
  int *countShare = new0 (arrayMax(clusterBlocks), int) ; /* number of hashes linking to code */
  int *minShare = new0 (arrayMax(clusterBlocks), int) ; /* lowest hash shared with code */
  Array clusterMin = arrayCreate (64, int) ;

  int n = nGoodHashes[code] ; if (!n) return ;
  int *minShareCount = new (n, int) ;
  U16 *g = goodHashes[code] ;
  
  for (i = 0 ; i < n ; ++i) b->clusHash[goodHashes[code][i]].subCluster = 0 ; /* wipe existing */

  b->nSubCluster = 0 ;
  //  if (!b->cluster) b->cluster = arrayCreate (32, ClusterInfo) ;
  b->pointToMin = 0.0 ;
  int nSubClustered = 0 ;
  for (i = 0 ; i < n ; ++i)	// shouldn't this be from i = 0 ?
    { x = b->clusHash[g[i]].hash ;
      nc = arr(hashDepth, x, U32) ;
      memset (minShareCount, 0, n*sizeof(int)) ;
					       
      for (j = 0 ; j < nc ; ++j)
	{ int cj = arr(hashCodes, x, U32*)[j] ;
	  if (cj == code) continue ;
	  if (!minShare[cj]) minShare[cj] = i + 1 ;
	  ++minShareCount[minShare[cj]-1] ;
	  ++countShare[cj] ;
	}

      int msMax = 0, msTot = 0, msBest ;
      for (j = 0 ; j < i ; ++j)
	{ if (minShareCount[j] > msMax)  { msBest = j ; msMax = minShareCount[j] ; }
	  msTot += minShareCount[j] ;
	}
      if (msMax >= clusterThreshold)
	{ ClusterHash *ch = &b->clusHash[g[msBest]] ;
	  if (!ch->subCluster)	/* create a new cluster */
	    { if (++b->nSubCluster > U8MAX) /* abandon this clustering */
		{ b->nSubCluster = 0 ; nSubClustered = 0 ;
		  for (j = 0 ; j < i ; ++j) b->clusHash[g[j]].subCluster = 0 ;
		  fprintf (stderr, "    code %d with %d good hashes has too many clusters\n",
			   code, nGoodHashes[code]) ;
		  break ;
		}
	      ch->subCluster = b->nSubCluster ; ++nSubClustered ;
	      array(clusterMin,b->nSubCluster,int) = msBest ;
	    }
	  b->clusHash[g[i]].subCluster = ch->subCluster ; ++nSubClustered ;
	  b->pointToMin += minShareCount[arr(clusterMin,ch->subCluster,int)] / (double)msTot ;
	}
    }

  if (b->nSubCluster)
    fprintf (outFile, "  code %d with %d hashes, %d good hashes, of which %d cluster into %d raw",
	     code, b->nHash, nGoodHashes[code], nSubClustered, b->nSubCluster) ;
  else
    fprintf (outFile, "  code %d with %d hashes, %d good hashes, 0 clusters\n",
	     code, b->nHash, nGoodHashes[code]) ;

  if (!b->nSubCluster) return ;	/* nothing to do */
  int *readMap = new0 (b->nRead, int), *trueCluster = new (b->nSubCluster+1, int) ;
  int *deadCluster = new0 (b->nSubCluster+1, int) ;
  for (i = 0 ; i <= b->nSubCluster ; ++i) trueCluster[i] = i ; /* initialisation */
  for (i = 0 ; i < b->nHash ; ++i)
    { int hashCluster = trueCluster[b->clusHash[i].subCluster] ; if (!hashCluster) continue ;
      int readCluster = trueCluster[readMap[b->clusHash[i].read]] ;
      if (hashCluster == readCluster) continue ;
      else if (!readCluster) readMap[b->clusHash[i].read] = hashCluster ;
      else			/* merge clusters hashCluster and readCluster */
	{ if (hashCluster > readCluster)
	    { int t = hashCluster ; hashCluster = readCluster ; readCluster = t ; }
	  for (j = 1 ; j <= b->nSubCluster ; ++j)
	    if (trueCluster[j] == readCluster) trueCluster[j] = hashCluster ;
	  deadCluster[readCluster] = 1 ;
	}
    }

  /* now compactify trueCluster, reusing deadCluster[] for the map removing dead indices */
  for (j = 1 ; j <= b->nSubCluster ; ++j) deadCluster[j] = deadCluster[j-1] + 1 - deadCluster[j] ;
  for (j = 1 ; j <= b->nSubCluster ; ++j) trueCluster[j] = deadCluster[trueCluster[j]] ;
  b->nSubCluster = deadCluster[b->nSubCluster] ;
  for (i = 0 ; i < b->nHash ; ++i) b->clusHash[i].subCluster = trueCluster[b->clusHash[i].subCluster] ;

  fprintf (outFile, " then %d merged clusters\n", b->nSubCluster) ; fflush (stdout) ;

  /* now go through countShare and identify clusters that share at least 3 hashes */

  Array h = arrayCreate (256, int) ;
  Array z = arrayCreate (1024, CountStruct) ;
  for (i = 0 ; i < arrayMax(clusterBlocks) ; ++i)
    { ++array(h, countShare[i], int) ;
      if (countShare[i] > 0)
	{ CountStruct *c = arrayp(z,arrayMax(z), CountStruct) ;
	  c->i = i ; c->count = countShare[i] ;
	}
    }
  histogramReport ("COUNT_SHARE", h) ; arrayDestroy (h) ;
  codeClusterReport (code, code+1) ;
  arraySort (z, compareCount) ;
  for (j = arrayMax(z) ; j-- ; ) /* run backwards through z */
    { CountStruct *cz = arrp(z, j, CountStruct) ;
      fprintf (outFile, "  SHARE code %d shared_hashes %d counts", cz->i, cz->count) ;
      int nA = 0, nB = 0 ;
      ClusterBlock *czb = arrp(clusterBlocks, cz->i, ClusterBlock) ;
      for (i = 0 ; i < czb->nHash ; ++i)
	if (cribType[czb->clusHash[i].hash] == CRIB_HTA) ++nA ;
	else if (cribType[czb->clusHash[i].hash] == CRIB_HTB) ++nB ;
      fprintf (outFile, " %d htA %d htB minShare %d %s\n", nA, nB, minShare[cz->i]-1,
	       cribText (b->clusHash[g[minShare[cz->i]-1]].hash)) ;
      if (cz->count < params.clusterThreshold) break ;
    }
  
  free (readMap) ; free(trueCluster) ; free(deadCluster) ;
  free (minShare) ; free (countShare) ; free (minShareCount) ;
  arrayDestroy (clusterMin) ;
}

/*************** end of file ***************/
