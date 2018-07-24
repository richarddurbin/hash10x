/*  File: minhash.h
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2018
 *-------------------------------------------------------------------
 * Description: header file for minhash package
 * Exported functions: see below
 * HISTORY:
 * Last edited: Apr 15 12:30 2018 (rd)
 * Created: Mon Mar  5 08:43:45 2018 (rd)
 *-------------------------------------------------------------------
 */

#include "utils.h"

typedef struct {
  int k ;			/* kmer */
  int w ;			/* window */
  U64 mask ;			/* 2*k bits */
  int shift1, shift2 ;
  U64 factor1, factor2 ;
  U64 patternRC[4] ;		/* one per base */
} MinHash ;

typedef struct {
  MinHash *m ;
  char *s, *sEnd ;     		/* sequence currently being hashed, end marker */
  U64 h, hRC ;			/* current hash values */
  U64 *hashBuf ;		/* buffer of length w holding hashes for current window */
  int base ;			/* start of buf in sequence */
  int iStart, iMin ;		/* position in buf of start of current window, next min */
  BOOL isDone ;
} MinHashRCiterator ;

MinHash *minHashCreate (int k, int w) ;
static void minHashDestroy (MinHash *m) { free (m) ; }

/* iterator to extract minHashes from a sequence */
/* NB sequence must continue to exist through the life of the iterator */
MinHashRCiterator *minHashRCiterator (MinHash *m, char *s, int len) ;
BOOL minHashRCnext (MinHashRCiterator *mi, U64 *u, int *pos) ; /* returns (u, pos) */
MinHashRCiterator *ranHashRCiterator (MinHash *m, char *s, int len) ;
BOOL ranHashRCnext (MinHashRCiterator *mi, U64 *u, int *pos) ; /* returns (u, pos) */
static void minHashRCiteratorDestroy (MinHashRCiterator *mi) { free (mi->hashBuf) ; free (mi) ; }

/******* end of file ********/
