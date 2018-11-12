/*  File: moshasm.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 11 17:18 2018 (rd109)
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
  int contained ;		/* outermost copy 1 hits are contained in this */
} RSinfo ;

typedef struct {		/* data structure for long read set */
  Moshset *ms ;
  int n ;
  U32* depth ;			/* indexed by ms; size ms->size */
  Array info ;			/* of RSinfo */
  Array hits ;			/* of lists of U32 indexes into ms */
  Array pos ;			/* of lists of U32 positions of hits */
} Readset ;

Readset *readsetCreate (Moshset *ms)
{
  Readset *rs = new0 (1, Readset) ;
  rs->ms = ms ;
  rs->depth = new0 (ms->size, U32) ;
  rs->info = arrayCreate (1<<16, RSinfo) ;
  rs->hits = arrayCreate (1<<16, U32*) ;
  rs->pos = arrayCreate (1<<16, U32*) ;
  return rs ;
}

void readsetDestroy (Readset *rs)
{
  free (rs->depth) ;
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
  if (fwrite(rs->depth,sizeof(U32),rs->ms->size,f) != rs->ms->size) die ("failed write depth") ;
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
  if (fread(rs->depth,sizeof(U32),rs->ms->size,f) != rs->ms->size) die ("failed write depth") ;
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

void readsetFastaRead (Readset *rs, char *filename)
{
  char *seq ;			/* ignore the name for now */
  int len ;
  Array hitsA = arrayCreate (1024, U32) ; /* reuse these to build the lists of hits and posns */
  Array posA = arrayCreate (1024, U32) ;
  
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
	      ++rs->depth[index] ;
	    }
	  else ++info->nMiss ;
	}
      seqhashRCiteratorDestroy (mi) ;
      if (info->nHit)
	{ U32 *x = array(rs->hits,rs->n,U32*) = new(info->nHit,U32) ;
	  memcpy (x, arrp(hitsA,0,U32), info->nHit) ;
	  x = array(rs->pos,rs->n,U32*) = new(info->nHit,U32) ;
	  memcpy (x, arrp(posA,0,U32), info->nHit) ;
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
      for (j = 1 ; j < 4 ; ++j) totCopy[j] += info->nCopy[j] ;
      if (info->nCopy[1] == 0) { ++nUnique0 ; lenUnique0 += info->len ; }
      else if (info->nCopy[1] == 1) { ++nUnique1 ; lenUnique1 += info->len ; }
    }
  fprintf (outFile, "RS %d sequences, total length %lld (av %.1f)\n",
	   rs->n, totLen, totLen/(double)rs->n) ;
  fprintf (outFile, "RS %lld mosh hits (%.1f bp/hit), frac hit %.2f\n",
	   totHit, totLen/(double)totHit, totHit/(double)(totMiss+totHit)) ;
  fprintf (outFile, "RS hit distribution %.2f copy1, %.2f copy2, %.2f copyM\n",
	   totCopy[1]/(double)totHit, totCopy[2]/(double)totHit, totCopy[3]/(double)totHit) ;
  fprintf (outFile, "RS %d with 0 unique hits, av length %.1f\n",
	   nUnique0, lenUnique0/(double)nUnique0) ;
  fprintf (outFile, "RS %d with 1 unique hits, av length %.1f\n",
	   nUnique1, lenUnique1/(double)nUnique1) ;

  U32 nCopy[4], hitCopy[4], hit2Copy[4] ;
  U64 depthCopy[4] ;
  for (j = 0 ; j < 4 ; ++j) nCopy[j] = hitCopy[j] = hit2Copy[j] = depthCopy[j] = 0 ;
  for (i = 1 ; i <= rs->ms->max ; ++i)
    { j = msCopy(rs->ms,i) ;
      ++nCopy[j] ;
      if (rs->depth[i] > 0) ++hitCopy[j] ;
      if (rs->depth[i] > 1) { ++hit2Copy[j] ; depthCopy[j] += rs->depth[i] ; }
    }
  fprintf (outFile, "RS frac hit, hit>1 (av) per mosh: copy1 %.2f, %.2f (%.1f), copy2 %.2f, %.2f (%.1f), copyM %.2f, %.2f (%.1f)\n", 
   hitCopy[1]/(double)nCopy[1], hit2Copy[1]/(double)nCopy[1], depthCopy[1]/(double)hit2Copy[1], 
   hitCopy[2]/(double)nCopy[2], hit2Copy[2]/(double)nCopy[2], depthCopy[2]/(double)hit2Copy[2], 
   hitCopy[3]/(double)nCopy[3], hit2Copy[3]/(double)nCopy[3], depthCopy[3]/(double)hit2Copy[3]) ;
}

/************************************************************/

void usage (void)
{ fprintf (stderr, "Usage: moshasm <commands>\n") ;
  fprintf (stderr, "Commands are executed in order - set parameters before using them!\n") ;
  fprintf (stderr, "  -v | --verbose : toggle verbose mode\n") ;
  fprintf (stderr, "  -t | --threads <number of threads for parallel ops> [%d]\n", numThreads) ;
  fprintf (stderr, "  -o | --output <output filename> : '-' for stdout\n") ;
  fprintf (stderr, "  -m | --mosh <mosh file>\n") ;
  fprintf (stderr, "  -f | --fasta <fasta file of reads>\n") ;
  fprintf (stderr, "  -w | --write <file stem> : writes assembly files\n") ;
  fprintf (stderr, "  -r | --read <file stem> : read assembly files\n") ;
  fprintf (stderr, "  -S | --stats : give readset stats\n") ;
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
    else if (ARGMATCH("-m","--mosh",2))
      { if (!(f = fopen (argv[-1], "r"))) die ("failed to open mosh file %s", argv[-1]) ;
	ms = moshsetRead (f) ;
      }
    else if (ARGMATCH("-f","--seqfile",2))
      { if (ms)
	  { if (rs) readsetDestroy (rs) ;
	    rs = readsetCreate (ms) ;
	    readsetFastaRead (rs, argv[-1]) ;
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
