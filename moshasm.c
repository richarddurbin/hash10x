/*  File: moshasm.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 19 10:47 2018 (rd109)
 * Created: Tue Nov  6 18:30:49 2018 (rd109)
 *-------------------------------------------------------------------
 */

#include "moshset.h"
#include "seqio.h"

#ifdef OMP
#include <omp.h>
#endif

int numThreads = 1 ;		/* default to serial - reset if multi-threaded */
FILE *outFile ;			/* initialise to stdout at start of main() */
BOOL isVerbose = FALSE ;

/* drop the read names */

typedef struct {
  int len ;
  int nHit ;
  int nMiss ;
  int nCopy[4] ;		/* number of copy 0, 1, 2, M hits */
  int contained ;		/* the current read seems to be fully contained in this one */
} RSinfo ;

typedef struct {		/* data structure for long read set */
  Moshset *ms ;
  int   n ;
  Array info ;			/* of RSinfo */
  Array hits ;			/* of lists of U32 indexes into ms */
  Array pos ;			/* of lists of U32 positions of hits */
  U32   **inv ;			/* inverse array, for each mosh, a pointer to the list of reads */
} Readset ;

Readset *readsetCreate (Moshset *ms)
{
  Readset *rs = new0 (1, Readset) ;
  rs->ms = ms ;
  rs->info = arrayCreate (1<<16, RSinfo) ;
  rs->hits = arrayCreate (1<<16, U32*) ;
  rs->pos = arrayCreate (1<<16, U32*) ;
  return rs ;
}

void readsetDestroy (Readset *rs)
{
  arrayDestroy (rs->info) ;
  int i ;
  for (i = 0 ; i < rs->n ; ++i)
    if (arr(rs->hits,i,U32*)) { free(arr(rs->hits,i,U32*)) ; free(arr(rs->pos,i,U32*)) ; }
  free (rs->hits) ; free (rs->pos) ; free (rs) ;
}

void readsetWrite (Readset *rs, char *root)
{
  FILE *f ;
  if (!(f = fopenTag (root, "mosh", "w"))) die ("can't open file %s.mosh", root) ;
  moshsetWrite (rs->ms, f) ; fclose (f) ;
  if (!(f = fopenTag (root, "readset", "w"))) die ("can't open file %s.readset", root) ;
  if (fwrite("RSMSHv1",8,1,f) != 1) die ("failed to write readset header") ;
  if (fwrite (&rs->n,sizeof(U32),1,f) != 1) die ("failed to write n") ;
  arrayWrite (rs->info, f) ;
  int i, n ;
  for (i = 0 ; i < rs->n ; ++i)
    if ((n = arrp(rs->info,i,RSinfo)->nHit))
      { if (fwrite(arr(rs->hits,i,U32*),sizeof(U32),n,f) != n) die ("failed write hits %d", n) ;
	if (fwrite(arr(rs->pos,i,U32*),sizeof(U32),n,f) != n) die ("failed write pos %d", n) ;
      }
  fclose (f) ;
}

Readset *readsetRead (char *root)
{
  FILE *f ;
  if (!(f = fopenTag (root, "mosh", "r"))) die ("can't open file %s.mosh", root) ;
  Moshset *ms = moshsetRead (f) ; fclose (f) ;
  Readset *rs = readsetCreate (ms) ;
  if (!(f = fopenTag (root, "readset", "r"))) die ("can't open file %s.readset", root) ;
  char name[8] ;
  if (fread (name,8,1,f) != 1) die ("failed to read readset header") ;
  if (strcmp (name, "RSMSHv1")) die ("bad readset header") ;
  if (fread (&rs->n,sizeof(U32),1,f) != 1) die ("failed to read n") ;
  arrayDestroy (rs->info) ; rs->info = arrayRead (f) ;
  int i, n ;
  for (i = 0 ; i < rs->n ; ++i)
    if ((n = arrp(rs->info,i,RSinfo)->nHit))
      { arr(rs->hits,i,U32*) = new (n, U32) ;
	if (fread(arr(rs->hits,i,U32*),sizeof(U32),n,f) != n) die ("failed read hits %d", n) ;
	arr(rs->pos,i,U32*) = new (n, U32) ;
	if (fread(arr(rs->pos,i,U32*),sizeof(U32),n,f) != n) die ("failed write pos %d", n) ;
      }
  fclose (f) ;
  return rs ;
}

void readsetFileRead (Readset *rs, char *filename)
{
  char *seq ;			/* ignore the name for now */
  int len ;
  Array hitsA = arrayCreate (1024, U32) ; /* reuse these to build the lists of hits and posns */
  Array posA = arrayCreate (1024, U32) ;

  memset (rs->ms->depth, 0, (rs->ms->max+1)*sizeof(U16)) ; /* rebuild depth from this file */
  dna2indexConv['N'] = dna2indexConv['n'] = 0 ; /* to get 2-bit encoding */
  SeqIO *si = seqIOopen (filename, dna2indexConv, FALSE) ;
  while (seqIOread (si))
    { RSinfo *info = arrp(rs->info, rs->n, RSinfo) ;
      info->len = si->seqLen ;
      SeqhashRCiterator *mi = moshRCiterator (rs->ms->hasher, sqioSeq(si), si->seqLen) ;
      hitsA = arrayReCreate (hitsA, 1024, U32) ;
      posA = arrayReCreate (posA, 1024, U32) ;
      U64 hash ; int pos ;
      while (moshRCnext (mi, &hash, &pos))
	{ U32 index = moshsetIndexFind (rs->ms, hash, FALSE) ;
	  if (index)
	    { array(hitsA,info->nHit,U32) = index ;
	      array(posA,info->nHit,U32) = pos ;
	      ++info->nCopy[msCopy(rs->ms,index)] ;
	      ++info->nHit ;
	      U16 *di = &rs->ms->depth[index] ; ++*di ; if (!*di) *di = U16MAX ;
	    }
	  else ++info->nMiss ;
	}
      seqhashRCiteratorDestroy (mi) ;
      if (info->nHit)
	{ U32 *x = array(rs->hits,rs->n,U32*) = new(info->nHit,U32) ;
	  memcpy (x, arrp(hitsA,0,U32), info->nHit*sizeof(U32)) ;
	  x = array(rs->pos,rs->n,U32*) = new(info->nHit,U32) ;
	  memcpy (x, arrp(posA,0,U32), info->nHit*sizeof(U32)) ;
	}
      ++rs->n ;
    }
  seqIOclose (si) ;

  arrayDestroy (hitsA) ; arrayDestroy (posA) ;
}

void readsetStats (Readset *rs)
{
  if (!rs || !rs->n) { fprintf (stderr, "stats called on empty readset\n") ; return ; }

  seqhashReport (rs->ms->hasher, outFile) ;
  
  int i, j ;
  int nUnique0 = 0, nUnique1 = 0 ;
  U64 totLen = 0, totHit = 0, totMiss = 0, lenUnique0 = 0, lenUnique1 = 0 ;
  U64 totCopy[4] ; totCopy[0] = totCopy[1] = totCopy[2] = totCopy[3] = 0 ;
  
  for (i = 0 ; i < rs->n ; ++i)
    { RSinfo *info = arrp(rs->info, i, RSinfo) ;
      totLen += info->len ;
      totHit += info->nHit ;
      totMiss += info->nMiss ;
      for (j = 0 ; j < 4 ; ++j) totCopy[j] += info->nCopy[j] ;
      if (info->nCopy[1] == 0) { ++nUnique0 ; lenUnique0 += info->len ; }
      else if (info->nCopy[1] == 1) { ++nUnique1 ; lenUnique1 += info->len ; }
    }
  fprintf (outFile, "RS %d sequences, total length %llu (av %.1f)\n",
	   rs->n, totLen, totLen/(double)rs->n) ;
  fprintf (outFile, "RS %llu mosh hits, %.1f bp/hit, frac hit %.2f, av hits/read %.1f\n",
	   totHit, totLen/(double)totHit, totHit/(double)(totMiss+totHit),
	   totHit/(double)rs->n) ;
  fprintf (outFile, "RS hit distribution %.2f copy0, %.2f copy1, %.2f copy2, %.2f copyM\n",
	   totCopy[0]/(double)totHit, totCopy[1]/(double)totHit,
	   totCopy[2]/(double)totHit, totCopy[3]/(double)totHit) ;
  U32 nUniqueMulti = rs->n - nUnique0 - nUnique1 ;
  fprintf (outFile, "RS num reads and av_len with 0 copy1 hits %d %.1f with 1 copy1 hits %d %.1f"
	   " >1 copy1 hits %d %.1f av copy1 hits %.1f\n",
	   nUnique0, lenUnique0/(double)nUnique0, nUnique1, lenUnique1/(double)nUnique1,
	   nUniqueMulti, (totLen - lenUnique0 - lenUnique1)/(double)nUniqueMulti,
	   (totCopy[1]-nUnique1)/(double)nUniqueMulti) ;
  U32 nCopy[4], hitCopy[4], hit2Copy[4] ;
  U64 depthCopy[4] ;
  for (j = 0 ; j < 4 ; ++j) nCopy[j] = hitCopy[j] = hit2Copy[j] = depthCopy[j] = 0 ;
  for (i = 1 ; i <= rs->ms->max ; ++i)
    { j = msCopy(rs->ms,i) ;
      ++nCopy[j] ;
      if (rs->ms->depth[i] > 0) ++hitCopy[j] ;
      if (rs->ms->depth[i] > 1) { ++hit2Copy[j] ; depthCopy[j] += rs->ms->depth[i] ; }
    }
  fprintf (outFile, "RS mosh frac hit hit>1 av: copy0 %.3f %.3f %.1f copy1 %.3f %.3f %.1f copy2 %.3f %.3f %.1f copyM %.3f %.3f %.1f\n", 
   hitCopy[0]/(double)nCopy[0], hit2Copy[0]/(double)nCopy[0], depthCopy[0]/(double)hit2Copy[0], 
   hitCopy[1]/(double)nCopy[1], hit2Copy[1]/(double)nCopy[1], depthCopy[1]/(double)hit2Copy[1], 
   hitCopy[2]/(double)nCopy[2], hit2Copy[2]/(double)nCopy[2], depthCopy[2]/(double)hit2Copy[2], 
   hitCopy[3]/(double)nCopy[3], hit2Copy[3]/(double)nCopy[3], depthCopy[3]/(double)hit2Copy[3]) ;
}

/************************************************************/

typedef struct {
  int nHit ;
  int j0, jn ;			/* first and last hit offsets in query */
} Overlap ;

void findOverlaps (Readset *rs, U32 i)
{
  RSinfo *rsi = arrp(rs->info, i, RSinfo) ;
  int j, k ;
  U32 *h = arr(rs->hits,i,U32*), *p = arr(rs->hits,i,U32*) ;
  Overlap *olap = new0(rs->n,Overlap) ;
  for (j = 0 ; j < rsi->nHit ; ++j, ++h, ++p)
    if (msIsCopy1 (rs->ms, *h))
      { U32 *r2 = rs->inv[*h] ;
	for (k = 0 ; k < rs->ms->depth[*h] ; ++k, ++r2)
	  { Overlap *o = &olap[*r2] ;
	    if (!o->nHit++) o->j0 = j ;
	    o->jn = j ;
	  }
      }
}

/************************************************************/

void usage (void)
{ fprintf (stderr, "Usage: moshasm <commands>\n") ;
  fprintf (stderr, "Commands are executed in order - set parameters before using them!\n") ;
  fprintf (stderr, "  -v | --verbose : toggle verbose mode\n") ;
  fprintf (stderr, "  -t | --threads <number of threads for parallel ops> [%d]\n", numThreads) ;
  fprintf (stderr, "  -o | --output <output filename> : '-' for stdout\n") ;
  fprintf (stderr, "  -m | --moshset <mosh file>\n") ;
  fprintf (stderr, "  -f | --seqfile <file of reads: fasta/q, can be gzipped, or binary>\n") ;
  fprintf (stderr, "  -w | --write <file stem> : writes assembly files\n") ;
  fprintf (stderr, "  -r | --read <file stem> : read assembly files\n") ;
  fprintf (stderr, "  -S | --stats : give readset stats\n") ;
  fprintf (stderr, "  -a1 | --assembly1 : find maximal overlaps / containments &\n") ;
}

int main (int argc, char *argv[])
{
  --argc ; ++argv ;		/* eat program name */

  outFile = stdout ;
  
  timeUpdate (stdout) ;		/* initialise timer */
#ifdef OMP
  numThreads = omp_get_max_threads () ;
  omp_set_num_threads (numThreads) ;
#endif

  if (!argc) usage () ;

  int i ;			/* generically useful variables */
  FILE *f ;
  Moshset *ms = 0 ;
  Readset *rs = 0 ;

  while (argc) {
    if (**argv != '-')
      die ("option/command %s does not start with '-': run without arguments for usage", *argv) ;
    fprintf (stderr, "COMMAND %s", *argv) ;
    for (i = 1 ; i < argc && *argv[i] != '-' ; ++i) fprintf (stderr, " %s", argv[i]) ;
    fputc ('\n', stderr) ;
    
#define ARGMATCH(x,y,n)	((!strcmp (*argv, x) || (!strcmp (*argv,y))) && argc >= n && (argc -= n, argv += n))
    if (ARGMATCH("-t","--threads",2))
      {
#ifdef OMP
	numThreads = atoi(argv[-1]) ;
	if (numThreads > omp_get_max_threads ()) numThreads = omp_get_max_threads () ;
	omp_set_num_threads (numThreads) ;
#else
	fprintf (stderr, "  can't set thread number - not compiled with OMP\n") ;
#endif
      }
    else if (ARGMATCH("-v","--verbose",1)) isVerbose = !isVerbose ;
    else if (ARGMATCH("-o","--output",2))
      { if (!strcmp (argv[-1], "-"))
	  outFile = stdout ;
	else if (!(outFile = fopen (argv[-1], "w")))
	  { fprintf (stderr, "can't open output file %s - resetting to stdout\n", argv[-1]) ;
	    outFile = stdout ;
	  }
      }
    else if (ARGMATCH("-m","--moshset",2))
      { if (!(f = fopen (argv[-1], "r"))) die ("failed to open mosh file %s", argv[-1]) ;
	if (ms) moshsetDestroy (ms) ;
	ms = moshsetRead (f) ;
      }
    else if (ARGMATCH("-f","--seqfile",2))
      { if (ms)
	  { if (rs) readsetDestroy (rs) ;
	    rs = readsetCreate (ms) ;
	    readsetFileRead (rs, argv[-1]) ;
	    fclose (f) ;
	  }
	else fprintf (stderr, "** need to read a moshset before a fasta file\n") ;
      }
    else if (ARGMATCH("-r","--read",2))
      { if (rs) readsetDestroy (rs) ;
	rs = readsetRead (argv[-1]) ;
      }
    else if (ARGMATCH("-w","--write",2))
      readsetWrite (rs, argv[-1]) ;
    else if (ARGMATCH("-S","--stats",1)) readsetStats (rs) ;
    else die ("unkown command %s - run without arguments for usage", *argv) ;

    timeUpdate (outFile) ;
  }

  fprintf (outFile, "total resources used: ") ; timeTotal (outFile) ;
  if (outFile != stdout) { printf ("total resources used: ") ; timeTotal (stdout) ; }
}

/************* end of file ************/
