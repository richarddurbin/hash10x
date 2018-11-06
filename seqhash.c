/*  File: seqhash.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description: seqhash package - uses random bit patterns for bases and rotate/XOR
 	compile with -DTEST to test, and -DDEBUG to debug, -DNDEBUG turns off asserts
	see test main() at end for standard usage pattern
 * Exported functions: see seqhash.h
 * HISTORY:
 * Last edited: Nov  6 11:39 2018 (rd109)
 * Created: Sat Feb 24 19:20:18 2018 (rd)
 *-------------------------------------------------------------------
 */

#include "seqhash.h"

static inline U64 rotateLeft (U64 x) { return (x << 2) | (x >> 62) ; }
static inline U64 rotateRight (U64 x) { return (x << 62) | (x >> 2) ; }

Seqhash *seqhashCreate (int k, int w)
{
  assert (sizeof (U64) == 8) ;
  Seqhash *sh = new0 (1, Seqhash) ;
  sh->k = k ; if (k < 1 || k >= 32) die ("seqhash k %d must be between 1 and 32\n", k) ;
  sh->w = w ; if (w < 1) die ("seqhash w %d must be positive\n", w) ;
  sh->mask = ((U64)1 << (2*k)) - 1 ;
  int i ;
  
  sh->factor1 = (random() << 32) | random() | 0x01 ;
  sh->shift1 = 64 - 2*k ;
  sh->factor2 = (random() << 32) | random() | 0x01 ;
  sh->shift2 = 2*k ;
  for (i = 0 ; i < 4 ; ++i) { sh->patternRC[i] = (3-i) ; sh->patternRC[i] <<= 2*(k-1) ; }
  return sh ;
}

#include <stdio.h>

void seqhashWrite (Seqhash *sh, FILE *f)
{ if (fwrite ("SQHSHv1",8,1,f) != 1) die ("failed to write seqhash header") ;
  if (fwrite (sh,sizeof(Seqhash),1,f) != 1) die ("failed to write seqhash") ;
}

Seqhash *seqhashRead (FILE *f)
{ Seqhash *sh = new (1, Seqhash) ;
  char name[8] ;
  if (fread (name,8,1,f) != 1) die ("failed to read seqhash header") ;
  if (strcmp (name, "SQHSHv1")) die ("seqhash read mismatch") ;
  if (fread (sh,sizeof(Seqhash),1,f) != 1) die ("failed to read seqhash") ;
  return sh ;
}

/************** basic hash functions *************/

static inline U64 hashFunc (U64 x, Seqhash *sh)
{ return ((x * sh->factor1) >> sh->shift1) ; } // + ((x * sh->factor2) >> sh->shift2) ; }

static inline U64 hashRC (SeqhashRCiterator *mi)
{ U64 hashF = hashFunc (mi->h, mi->sh) ;
  U64 hashR = hashFunc (mi->hRC, mi->sh) ;
#ifdef DEBUG
  printf ("hashRC: h %lx hRC %lx hashF %lx hashR %lx\n", mi->h, mi->hRC, hashF, hashR) ;
#endif
  if (hashF < hashR) return hashF ; else return hashR ;
}

static inline U64 advanceHashRC (SeqhashRCiterator *mi)
{ Seqhash *sh = mi->sh ;
  if (mi->s < mi->sEnd)
    { mi->h = ((mi->h << 2) & sh->mask) | *(mi->s) ;
      mi->hRC = (mi->hRC >> 2) | sh->patternRC[*(mi->s)] ;
      return hashRC (mi) ;
    }
  else
    return U64MAX ;
}

/************ iterator to run across a sequence ***********/

SeqhashRCiterator *minimizerRCiterator (Seqhash *sh, char *s, int len)
{
  assert (s && len >= 0) ;
  SeqhashRCiterator *mi = (SeqhashRCiterator*) mycalloc (sizeof(SeqhashRCiterator), 1) ;
  mi->sh = sh ;
  mi->s = s ; mi->sEnd = s + len ;
  mi->hashBuf = (U64*) myalloc (sh->w * sizeof(U64)) ;
  if (len < sh->k) { mi->isDone = TRUE ; return mi ; } /* edge case */

  int i ;			/* preinitialise the hashes for the first kmer */
  for (i = 0 ; i < sh->k ; ++i, ++mi->s)
    { mi->h = (mi->h << 2) | *mi->s ;
      mi->hRC = (mi->hRC >> 2) | sh->patternRC[*(mi->s)] ;
    }
  
  /* store first w hashes in hashBuf and set ->iMin */
  U64 min = hashRC (mi) ;
  mi->iMin = 0 ;
  for (i = 1 ; i < sh->w ; ++i, ++mi->s)
    { mi->hashBuf[i] = advanceHashRC (mi) ;
      if (mi->hashBuf[i] < min) { min = mi->hashBuf[i] ; mi->iMin = i ; }
    }

  return mi ;
}

BOOL minimizerRCnext (SeqhashRCiterator *mi, U64 *u, int *pos) /* returns (u, pos) */
{
  if (mi->isDone) return FALSE ; /* we are done */

#ifdef DEBUG
  printf ("base %d, iStart %d, iMin %d\n", mi->base, mi->iStart, mi->iMin) ;
  int j ; for (j = 0 ; j < mi->sh->w ; ++j) printf ("  %x", mi->hashBuf[j]) ;
  printf ("\n") ;
#endif

  assert (u && pos) ;
  *pos = mi->base + mi->iMin ; if (mi->iMin < mi->iStart) *pos += mi->sh->w ;
  *u = mi->hashBuf[mi->iMin] ;
  if (mi->s >= mi->sEnd) { mi->isDone = TRUE ; return TRUE ; }

  int i ;	    		/* next update hashBuf splitting into two cases */
  U64 min = *u ;    /* save this here for end case - see below */
  if (mi->iMin >= mi->iStart)
    for (i = mi->iStart ; i <= mi->iMin ; ++i, ++mi->s) mi->hashBuf[i] = advanceHashRC (mi) ;
  else
    { for (i = mi->iStart ; i < mi->sh->w ; ++i, ++mi->s) mi->hashBuf[i] = advanceHashRC (mi) ;
      mi->base += mi->sh->w ;
      for (i = 0 ; i <= mi->iMin ; ++i, ++mi->s) mi->hashBuf[i] = advanceHashRC (mi) ;
    }
  mi->iStart = mi->iMin + 1 ;
  if (mi->iStart == mi->sh->w) { mi->iStart = 0 ; mi->base += mi->sh->w ; }

  /* finally find new min to set up for next call */
  if (mi->hashBuf[mi->iMin] != U64MAX) /* there was a full new window */
    min = U64MAX ;
  else				/* otherwise, keep the last min */
    mi->iMin = -1 ;
  for (i = 0 ; i < mi->sh->w ; ++i)
    if (mi->hashBuf[i] < min) { min = mi->hashBuf[i] ; mi->iMin = i ; }
  if (mi->iMin == -1)		/* our old min was not beaten - we are done */
    mi->isDone = TRUE ;
  
  return TRUE ;
}

SeqhashRCiterator *moshRCiterator (Seqhash *sh, char *s, int len)
{
  assert (s && len >= 0) ;
  SeqhashRCiterator *mi = (SeqhashRCiterator*) mycalloc (sizeof(SeqhashRCiterator), 1) ;
  mi->sh = sh ;
  mi->s = s ; mi->sEnd = s + len ;
  mi->hashBuf = (U64*) myalloc (sh->w * sizeof(U64)) ;
  if (len < sh->k) { mi->isDone = TRUE ; return mi ; } /* edge case */

  int i ;			/* preinitialise the hashes for the first kmer */
  for (i = 0 ; i < sh->k ; ++i, ++mi->s)
    { mi->h = (mi->h << 2) | *mi->s ;
      mi->hRC = (mi->hRC >> 2) | sh->patternRC[*(mi->s)] ;
    }
  
  U64 hash = hashRC(mi) ;
  while ((hash % sh->w) && mi->s < mi->sEnd)
    { hash = advanceHashRC (mi) ; ++mi->iMin ; ++mi->s ; }
  if (!(hash % sh->w)) *mi->hashBuf = hash ;
  else mi->isDone = TRUE ;

  return mi ;
}

BOOL moshRCnext (SeqhashRCiterator *mi, U64 *u, int *pos) /* returns (u, pos) */
{
  if (mi->isDone) return FALSE ; /* we are done */

  *u = *mi->hashBuf ; *pos = mi->iMin ;

  if (mi->s >= mi->sEnd) { mi->isDone = TRUE ; return TRUE ; }
    
  U64 hash = advanceHashRC (mi) ; ++mi->iMin ; ++mi->s ;
  int w = mi->sh->w ;
  while ((hash % w) && mi->s < mi->sEnd)
    { hash = advanceHashRC (mi) ; ++mi->iMin ; ++mi->s ; }
  if (!(hash % w)) *mi->hashBuf = hash ;
  else mi->isDone = TRUE ;

  return TRUE ;
}

/************** short test program, illustrating standard usage *************/

#ifdef TEST

#include "readseq.h"

int main (int argc, char *argv[])
{
  char *seq, *id ;
  int len ;
  U64 u ; int pos ;
  
  Seqhash *sh = seqhashCreate (16,32) ;
  while (readSequence (stdin, dna2indexConv, &seq, &id, 0, &len))
    { printf ("\nread sequence %s length %d\n", id, len) ;
      SeqhashRCiterator *si = moshRCiterator (sh, seq, len) ;
      while (moshRCnext (si, &u, &pos)) printf ("\t%08lx\t%d\n", u, pos) ;
      seqhashRCiteratorDestroy (si) ;
    }
}

#endif

/**************** end of file ****************/
