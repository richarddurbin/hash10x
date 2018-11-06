/*  File: moshmap.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Nov  6 22:10 2018 (rd109)
 * Created: Sat Oct 27 20:37:44 2018 (rd109)
 *-------------------------------------------------------------------
 */

#include "moshset.h"
#include "readseq.h"
#include <ctype.h>
#ifdef OMP
#include <omp.h>
#endif

int numThreads = 1 ;		/* default to serial - reset if multi-threaded */
FILE *outFile ;			/* initialise to stdout at start of main() */
BOOL isVerbose = FALSE ;

struct {
  int k ;
  int w ;
  int s ;
  int B ;
} params ;

/*******************************************************************/

/* object to hold a reference sequence */
typedef struct {
  Moshset *ms ;
  U32 size ;			/* size of index, offset, id */
  U32 max ;		      	/* number of indices in the reference */
  U32 *index ; 			/* consecutive Moshset index values */
  U32 *offset ;			/* base offset within relevant sequence */
  U32 *id ;			/* index into refDict */
  U32 *depth ;		  	/* number of times this index is seen; size ms->size */
  U32 *rev ;			/* reverse index from mosh to ref; size max */
  U32 *loc ;			/* offsets into rev for each mosh ; size ms->max */
  DICT *dict ;			/* set of reference names */
  Array len ;			/* of U32, indexed over refId */
} Reference ;
  
Reference *referenceCreate (Moshset *ms, U32 size)
{
  if (!ms || !ms->size) die ("moshset must be initialised before reference") ;
  if (!size) die ("refCreate must have size > 0") ;
  Reference *ref = new0 (1, Reference) ;
  ref->ms = ms ;
  ref->depth = ms->max ? new0 (ms->max, U32) : new0 (ms->size, U32) ;
  ref->size = size ;
  ref->index = new (size, U32) ;
  ref->offset = new (size, U32) ;
  ref->id = new (size, U32) ;
  ref->dict = dictCreate (1024) ;
  ref->len = arrayCreate (1024, U32) ;
  /* NB don't create loc and rev here */
  return ref ;
}

void referenceDestroy (Reference *ref)
{ free (ref->depth) ;
  if (ref->loc) free (ref->loc) ; if (ref->rev) free (ref->rev) ;
  free (ref->index) ; free (ref->offset) ; free (ref->id) ;
  dictDestroy (ref->dict) ; arrayDestroy (ref->len) ;
  free (ref) ;
}

void referencePack (Reference *ref)
{
  resize (ref->depth, ref->ms->max, U32) ;
  resize (ref->index, ref->max, U32) ;
  resize (ref->offset, ref->max, U32) ;
  resize (ref->id, ref->max, U32) ;
  ref->size = ref->max ;
  
  ref->rev = new (ref->size, U32) ;			/* now build rev and loc */
  ref->loc = new (ref->ms->max + 1, U32) ; ref->loc[0] = 0 ;
  int i ;
  for (i = 1 ; i <= ref->ms->max ; ++i)
    { ref->loc[i] = ref->loc[i-1] + ref->depth[i] ; ref->depth[i] = 0 ; } /* reuse depth below */
  U32 *ri = ref->index ;
  for (i = 0 ; i < ref->max ; ++i, ++ri) ref->rev[ref->loc[*ri] + ref->depth[*ri]++] = i ;
}

void referenceFastaRead (Reference *ref, FILE *f, BOOL isAdd)
{
  char *seq, *seqName ;
  int len ;
  U64 totLen = 0 ;
  
  dna2indexConv['N'] = dna2indexConv['n'] = 0 ; /* to get 2-bit encoding */
  while (readSequence (f, dna2indexConv, &seq, &seqName, 0, &len))
    { int id ;
      if (!dictAdd (ref->dict, seqName, &id)) die ("duplicate ref sequence name %s", seqName) ;
      array (ref->len, id, int) = len ;
      totLen += len ;
      SeqhashRCiterator *mi = moshRCiterator (ref->ms->hasher, seq, len) ;
      U64 hash ; int pos ;
      while (moshRCnext (mi, &hash, &pos))
	{ U32 index = moshsetIndexFind (ref->ms, hash, isAdd) ;
	  if (index)
	    { if (ref->max+1 >= ref->size) die ("reference size overflow") ;
	      ref->index[ref->max] = index ;
	      ++ref->depth[index] ;
	      ref->offset[ref->max] = pos ;
	      ref->id[ref->max] = id ;
	      ++ref->max ;
	    }
	}
      seqhashRCiteratorDestroy (mi) ;
      free (seq) ; free (seqName) ;
    }
  fclose (f) ;

  fprintf (outFile, "  %d hashes from %d reference sequences, total length %lld\n",
	   ref->max, dictMax(ref->dict), totLen) ;
  int i ; U32 *d = &ref->depth[1] ; U32 n1 = 0, n2 = 0, nM = 0 ;
  for (i = 1 ; i <= ref->ms->max ; ++i, ++d)
    if (*d == 1) { msSetCopy1 (ref->ms, i) ; ++n1 ; }
    else if (*d == 2) { msSetCopy2 (ref->ms, i) ; ++n2 ; }
    else { msSetCopyM (ref->ms, i) ; ++nM ; }
  fprintf (outFile, "  %d copy 1, %d copy 2, %d multiple\n", n1, n2, nM) ;

  if (isAdd) moshsetPack (ref->ms) ;
  referencePack (ref) ;
}

void referenceWrite (Reference *ref, char *root)
{
  FILE *f ;
  if (!(f = fopenTag (root, "mosh", "w"))) die ("failed to open %s.mosh to write", root) ;
  moshsetWrite (ref->ms, f) ;
  fclose (f) ;
  if (!(f = fopenTag (root, "ref", "w"))) die ("failed to open %s.ref to write", root) ;
  if (fwrite ("RFMSHv1",8,1,f) != 1) die ("failed to write reference header") ;
  U32 size = ref->max ; if (fwrite (&size,sizeof(U32),1,f) != 1) die ("failed to write size") ;
  if (fwrite (&ref->max,sizeof(U32),1,f) != 1) die ("failed to write max") ;
  if (fwrite (ref->index,sizeof(U32),size,f) != size) die ("failed write ref index") ;
  if (fwrite (ref->offset,sizeof(U32),size,f) != size) die ("failed write ref offset") ;
  if (fwrite (ref->id,sizeof(U32),size,f) != size) die ("failed write ref id") ;
  if (fwrite (ref->depth,sizeof(U32),ref->ms->max+1,f) != ref->ms->max+1) die ("fail depth") ;
  if (fwrite (ref->rev,sizeof(U32),size,f) != size) die ("fail rev") ;
  if (fwrite (ref->loc,sizeof(U32),ref->ms->max+1,f) != ref->ms->max+1) die ("fail loc") ;
  if (!arrayWrite (ref->len, f)) die ("failed write ref len") ;
  if (!dictWrite (ref->dict, f)) die ("failed write ref dict") ;
  /* write the dict last, because it can terminate in arbitrary byte alignment */
  fclose (f) ;
}  

Reference *referenceRead (char *root)
{
  FILE *f ;
  if (!(f = fopenTag (root, "mosh", "r"))) die ("failed to open %s.mosh to read", root) ;
  Moshset *ms = moshsetRead (f) ;
  fclose (f) ;
  if (!(f = fopenTag (root, "ref", "r"))) die ("failed to open %s.ref to read", root) ;
  char name[8] ;
  if (fread (name,8,1,f) != 1) die ("failed to read reference header") ;
  if (strcmp (name, "RFMSHv1")) die ("bad reference header") ;
  U32 size ; if (fread (&size,sizeof(U32),1,f) != 1) die ("failed to read size") ;
  Reference *ref = referenceCreate (ms, size) ;
  if (fread (&ref->max,sizeof(U32),1,f) != 1) die ("failed to read max") ;
  if (fread (ref->index,sizeof(U32),size,f) != size) die ("failed read ref index") ;
  if (fread (ref->offset,sizeof(U32),size,f) != size) die ("failed read ref offset") ;
  if (fread (ref->id,sizeof(U32),size,f) != size) die ("failed read ref id") ;
  if (fread (ref->depth,sizeof(U32),ref->ms->max+1,f) != ref->ms->max+1) die ("fail depth") ;
  ref->rev = new (size, U32) ; ref->loc = new (ref->ms->max+1, U32) ;
  if (fread (ref->rev,sizeof(U32),size,f) != size) die ("fail rev") ;
  if (fread (ref->loc,sizeof(U32),ref->ms->max+1,f) != ref->ms->max+1) die ("fail loc") ;
  arrayDestroy (ref->len) ; if (!(ref->len = arrayRead (f))) die ("failed read ref len") ;
  dictDestroy (ref->dict) ; if (!(ref->dict = dictRead (f))) die ("failed read ref dict") ;
  fclose (f) ;
  return ref ;
}

/************************************************************/

typedef struct { U32 index ; U32 pos ; } Seed ;

void queryProcess (Reference *ref, FILE *f)
{
  char *seq, *seqName ;
  int len ;
  Array seeds = 0 ;

  dna2indexConv['N'] = dna2indexConv['n'] = 0 ; /* to get 2-bit encoding */
  while (readSequence (f, dna2indexConv, &seq, &seqName, 0, &len))
    { SeqhashRCiterator *mi = moshRCiterator (ref->ms->hasher, seq, len) ;
      U64 hash ; int pos ;
      seeds = arrayReCreate (seeds, 1024, Seed) ;
      int missed = 0, unique = 0 ;
      while (moshRCnext (mi, &hash, &pos))
	{ U32 index = moshsetIndexFind (ref->ms, hash, FALSE) ;
	  Seed *s = arrayp(seeds,arrayMax(seeds),Seed) ;
	  s->index = index ; s->pos = pos ;
	  if (index && ref->depth[index] == 1) ++unique ;
	  else if (!index) ++missed ;
	}
      seqhashRCiteratorDestroy (mi) ;
      fprintf (outFile, "Q\t%s\t%d\t%d miss, %d unique, %d multi, %.2f hit\n",
	       seqName, len, missed, unique, arrayMax(seeds)-unique-missed,
	       (arrayMax(seeds)-missed)/(double)(arrayMax(seeds))) ;

      int i ;
      U32 loc0 = 0, locN, i0, iN ;
      int nHit = 0 ;
      for (i = 0 ; i < arrayMax(seeds) ; ++i)
	{ Seed *s = arrp (seeds, i, Seed) ;
	  if (s->index && ref->depth[s->index] == 1) /* only consider unique hits for now */
	    { U32 loc = ref->rev[ref->loc[s->index]] ;
	      if (isVerbose)
		printf ("  %6d\t%s\t%8d\t%8d\t%8d\n", s->pos, dictName(ref->dict,ref->id[loc]),
			ref->offset[loc], ref->offset[loc] - s->pos, ref->offset[loc] + s->pos) ;
	      BOOL endBlock = FALSE ;
	      if (!loc0 || ref->id[loc] != ref->id[loc0]) endBlock = TRUE ;
	      else if (loc0 < locN)
		{ if (loc < locN) endBlock = TRUE ;
		  int d = locN - loc0 - iN + i0 ; if (d > 50 || d < -50) endBlock = TRUE ;
		}
	      else if (loc0 > locN)
		{ if (loc > locN) endBlock = TRUE ;
		  int d = loc0 - locN - iN + i0 ; if (d > 50 || d < -50) endBlock = TRUE ;
		}
	      if (endBlock)
		{ if (nHit > 2)
		    fprintf (outFile, "M\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%.2f\t%.2f\n",
			     seqName, arrp(seeds,i0,Seed)->pos, arrp(seeds,iN,Seed)->pos, len,
			     dictName(ref->dict,ref->id[loc0]),
			     ref->offset[loc0], ref->offset[locN],
			     nHit, nHit / (double)((locN>loc0) ? (locN-loc0) : (loc0-locN)),
			     nHit / (double)unique) ;
		  nHit = 0 ; loc0 = loc ; i0 = i ;
		}
	      ++nHit ; locN = loc ; iN = i ;
	    }
	}
      if (nHit > 2)
	fprintf (outFile, "M\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%.2f\t%.2f\n",
		 seqName, arrp(seeds,i0,Seed)->pos, arrp(seeds,iN,Seed)->pos, len,
		 dictName(ref->dict,ref->id[loc0]),
		 ref->offset[loc0], ref->offset[locN],
		 nHit, nHit / (double)((locN>loc0) ? (locN-loc0) : (loc0-locN)),
		 nHit / (double) unique) ;

      free (seq) ; free (seqName) ;
    }

  arrayDestroy (seeds) ;
}

/************************************************************/

void usage (void)
{ fprintf (stderr, "Usage: moshmap <commands>\n") ;
  fprintf (stderr, "Commands are executed in order - set parameters before using them!\n") ;
  fprintf (stderr, "  -K | --kmer <kmer size> [%d]\n", params.k) ;
  fprintf (stderr, "  -W | --window <window> [%d]\n", params.w) ;
  fprintf (stderr, "  -S | --seed <random number seed> [%d]\n", params.s) ;
  fprintf (stderr, "  -B | --tableBits <hash index table bitcount> [%d]\n", params.B) ;
  fprintf (stderr, "  -v | --verbose : toggle verbose mode\n") ;
  fprintf (stderr, "  -t | --threads <number of threads for parallel ops> [%d]\n", numThreads) ;
  fprintf (stderr, "  -o | --output <output filename> : '-' for stdout\n") ;
  fprintf (stderr, "  -f | --referenceFasta <reference fasta file>\n") ;
  fprintf (stderr, "  -w | --referenceWrite <file stem> : writes reference hash files\n") ;
  fprintf (stderr, "  -r | --referenceRead <file stem> : read reference hash files\n") ;
  fprintf (stderr, "  -q | --query <query fasta file>\n") ;
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

  /* command line parameters */
  params.k = 19 ;
  params.w = 31 ;
  params.s = 17 ;
  params.B = 28 ;

  if (!argc) usage () ;

  int i ;			/* generically useful variables */
  FILE *f ;

  Reference *ref = 0 ;

  while (argc) {
    if (**argv != '-')
      die ("option/command %s does not start with '-': run without arguments for usage", *argv) ;
    fprintf (stderr, "COMMAND %s", *argv) ;
    for (i = 1 ; i < argc && *argv[i] != '-' ; ++i) fprintf (stderr, " %s", argv[i]) ;
    fputc ('\n', stderr) ;
    
#define ARGMATCH(x,y,n)	((!strcmp (*argv, x) || (!strcmp (*argv,y))) && argc >= n && (argc -= n, argv += n))
    if (ARGMATCH("-K","--kmer",2)) params.k = atoi(argv[-1]) ;
    else if (ARGMATCH("-W","--window",2)) params.w = atoi(argv[-1]) ;
    else if (ARGMATCH("-S","--seed",2)) params.s = atoi(argv[-1]) ;
    else if (ARGMATCH("-B","--tableBits",2)) params.B = atoi(argv[-1]) ;
    else if (ARGMATCH("-t","--threads",2))
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
    else if (ARGMATCH("-f","--referenceFasta",2))
      { if (!(f = fopen (argv[-1], "r"))) die ("failed to open fasta file %s", argv[-1]) ;
	srandom (params.s) ;
	if (params.k <= 0 || params.w <= 0) die ("k %d, w %d must be > 0", params.k, params.w) ;
	Seqhash *hasher = seqhashCreate (params.k, params.w) ;
	fprintf (outFile, "  moshmap initialised with k = %d, w = %d, random seed = %d\n",
		 params.k, params.w, params.s) ;
	Moshset *ms = moshsetCreate (hasher, params.B, 0) ;
	ref = referenceCreate (ms, 1 << 26) ;
	referenceFastaRead (ref, f, TRUE) ;
      }
    else if (ARGMATCH("-q","--query",2))
      { if (!ref) die ("need to read a reference before processing query sequences") ;
	if (!(f = fopen (argv[-1], "r"))) die ("failed to open query file %s", argv[-1]) ;
	queryProcess (ref, f) ;
      }
    else if (ARGMATCH("-r","--referenceRead",2))
      { if (ref) referenceDestroy (ref) ;
	ref = referenceRead (argv[-1]) ;
      }
    else if (ARGMATCH("-w","--referenceWrite",2))
      referenceWrite (ref, argv[-1]) ;
    timeUpdate (outFile) ;
  }

  fprintf (outFile, "total resources used: ") ; timeTotal (outFile) ;
  if (outFile != stdout) { printf ("total resources used: ") ; timeTotal (stdout) ; }
}

/************* end of file ************/
