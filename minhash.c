/*  File: minhash.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2018
 *-------------------------------------------------------------------
 * Description: minhash package - uses random bit patterns for bases and rotate/XOR
 	compile with -DTEST to test, and -DDEBUG to debug, -DNDEBUG turns off asserts
 * Exported functions: see minhash.h
 * HISTORY:
 * Last edited: Jul 24 09:55 2018 (rd)
 * Created: Sat Feb 24 19:20:18 2018 (rd)
 *-------------------------------------------------------------------
 */

#include "minhash.h"

static inline U64 rotateLeft (U64 x) { return (x << 2) | (x >> 62) ; }
static inline U64 rotateRight (U64 x) { return (x << 62) | (x >> 2) ; }

MinHash *minHashCreate (int k, int w)
{
  assert (sizeof (U64) == 8) ;
  MinHash *m = (MinHash*) mycalloc (1, sizeof(MinHash)) ;
  m->k = k ; if (k < 1 || k >= 32) die ("minHash k %d must be between 1 and 32\n", k) ;
  m->w = w ; if (w < 1) die ("minHash w %d must be positive\n", w) ;
  m->mask = ((U64)1 << (2*k)) - 1 ;
  int i ;
  
  m->factor1 = (random() << 32) | random() | 0x01 ;
  m->shift1 = 64 - 2*k ;
  m->factor2 = (random() << 32) | random() | 0x01 ;
  m->shift2 = 2*k ;
  for (i = 0 ; i < 4 ; ++i) { m->patternRC[i] = (3-i) ; m->patternRC[i] <<= 2*(k-1) ; }
  return m ;
}

static inline U64 hashFunc (U64 x, MinHash *m)
{ return ((x * m->factor1) >> m->shift1) ; } // + ((x * m->factor2) >> m->shift2) ; }

static inline U64 hashRC (MinHashRCiterator *mi)
{ U64 hashF = hashFunc (mi->h, mi->m) ;
  U64 hashR = hashFunc (mi->hRC, mi->m) ;
#ifdef DEBUG
  printf ("hashRC: h %lx hRC %lx hashF %lx hashR %lx\n", mi->h, mi->hRC, hashF, hashR) ;
#endif
  if (hashF < hashR) return hashF ;
  else return hashR ;
}

static inline U64 advanceHashRC (MinHashRCiterator *mi)
{ MinHash *m = mi->m ;
  if (mi->s < mi->sEnd)
    { mi->h = ((mi->h << 2) & m->mask) | *(mi->s) ;
      mi->hRC = (mi->hRC >> 2) | m->patternRC[*(mi->s)] ;
      return hashRC (mi) ;
    }
  else
    return U64MAX ;
}

MinHashRCiterator *minHashRCiterator (MinHash *m, char *s, int len)
{
  assert (s && len >= 0) ;
  MinHashRCiterator *mi = (MinHashRCiterator*) mycalloc (sizeof(MinHashRCiterator), 1) ;
  mi->m = m ;
  mi->s = s ; mi->sEnd = s + len ;
  mi->hashBuf = (U64*) myalloc (m->w * sizeof(U64)) ;
  if (len < m->k) { mi->isDone = TRUE ; return mi ; } /* edge case */

  int i ;			/* preinitialise the hashes for the first kmer */
  for (i = 0 ; i < m->k ; ++i, ++mi->s)
    { mi->h = (mi->h << 2) | *mi->s ;
      mi->hRC = (mi->hRC >> 2) | m->patternRC[*(mi->s)] ;
    }
  
  /* store first w hashes in hashBuf and set ->iMin */
  U64 min = hashRC (mi) ;
  mi->iMin = 0 ;
  for (i = 1 ; i < m->w ; ++i, ++mi->s)
    { mi->hashBuf[i] = advanceHashRC (mi) ;
      if (mi->hashBuf[i] < min) { min = mi->hashBuf[i] ; mi->iMin = i ; }
    }

  return mi ;
}

BOOL minHashRCnext (MinHashRCiterator *mi, U64 *u, int *pos) /* returns (u, pos) */
{
  if (mi->isDone) return FALSE ; /* we are done */

#ifdef DEBUG
  printf ("base %d, iStart %d, iMin %d\n", mi->base, mi->iStart, mi->iMin) ;
  int j ; for (j = 0 ; j < mi->m->w ; ++j) printf ("  %x", mi->hashBuf[j]) ;
  printf ("\n") ;
#endif

  assert (u && pos) ;
  *pos = mi->base + mi->iMin ; if (mi->iMin < mi->iStart) *pos += mi->m->w ;
  *u = mi->hashBuf[mi->iMin] ;
  if (mi->s >= mi->sEnd) { mi->isDone = TRUE ; return TRUE ; }

  int i ;	    		/* next update hashBuf splitting into two cases */
  U64 min = *u ;    /* save this here for end case - see below */
  if (mi->iMin >= mi->iStart)
    for (i = mi->iStart ; i <= mi->iMin ; ++i, ++mi->s) mi->hashBuf[i] = advanceHashRC (mi) ;
  else
    { for (i = mi->iStart ; i < mi->m->w ; ++i, ++mi->s) mi->hashBuf[i] = advanceHashRC (mi) ;
      mi->base += mi->m->w ;
      for (i = 0 ; i <= mi->iMin ; ++i, ++mi->s) mi->hashBuf[i] = advanceHashRC (mi) ;
    }
  mi->iStart = mi->iMin + 1 ;
  if (mi->iStart == mi->m->w) { mi->iStart = 0 ; mi->base += mi->m->w ; }

  /* finally find new min to set up for next call */
  if (mi->hashBuf[mi->iMin] != U64MAX) /* there was a full new window */
    min = U64MAX ;
  else				/* otherwise, keep the last min */
    mi->iMin = -1 ;
  for (i = 0 ; i < mi->m->w ; ++i)
    if (mi->hashBuf[i] < min) { min = mi->hashBuf[i] ; mi->iMin = i ; }
  if (mi->iMin == -1)		/* our old min was not beaten - we are done */
    mi->isDone = TRUE ;
  
  return TRUE ;
}

MinHashRCiterator *ranHashRCiterator (MinHash *m, char *s, int len)
{
  assert (s && len >= 0) ;
  MinHashRCiterator *mi = (MinHashRCiterator*) mycalloc (sizeof(MinHashRCiterator), 1) ;
  mi->m = m ;
  mi->s = s ; mi->sEnd = s + len ;
  mi->hashBuf = (U64*) myalloc (m->w * sizeof(U64)) ;
  if (len < m->k) { mi->isDone = TRUE ; return mi ; } /* edge case */

  int i ;			/* preinitialise the hashes for the first kmer */
  for (i = 0 ; i < m->k ; ++i, ++mi->s)
    { mi->h = (mi->h << 2) | *mi->s ;
      mi->hRC = (mi->hRC >> 2) | m->patternRC[*(mi->s)] ;
    }
  
  U64 hash = hashRC(mi) ;
  while ((hash % m->w) && mi->s < mi->sEnd)
    { hash = advanceHashRC (mi) ; ++mi->iMin ; ++mi->s ; }
  if (!(hash % m->w)) *mi->hashBuf = hash ;
  else mi->isDone = TRUE ;

  return mi ;
}

BOOL ranHashRCnext (MinHashRCiterator *mi, U64 *u, int *pos) /* returns (u, pos) */
{
  if (mi->isDone) return FALSE ; /* we are done */

  *u = *mi->hashBuf ; *pos = mi->iMin ;

  if (mi->s >= mi->sEnd) { mi->isDone = TRUE ; return TRUE ; }
    
  U64 hash = advanceHashRC (mi) ; ++mi->iMin ; ++mi->s ;
  int w = mi->m->w ;
  while ((hash % w) && mi->s < mi->sEnd)
    { hash = advanceHashRC (mi) ; ++mi->iMin ; ++mi->s ; }
  if (!(hash % w)) *mi->hashBuf = hash ;
  else mi->isDone = TRUE ;

  return TRUE ;
}

#ifdef TEST

#include "readseq.h"

int main (int argc, char *argv[])
{
  char *seq, *id ;
  int len ;
  U64 u ; int pos ;
  
  MinHash *m = minHashCreate (16,32) ;
  while (readSequence (stdin, dna2indexConv, &seq, &id, 0, &len))
    { printf ("\nread sequence %s length %d\n", id, len) ;
      MinHashRCiterator *mi = ranHashRCiterator (m, seq, len) ;
      while (ranHashRCnext (mi, &u, &pos)) printf ("\t%08lx\t%d\n", u, pos) ;
      minHashRCiteratorDestroy (mi) ;
    }
}

#endif

/**************** end of file ****************/
