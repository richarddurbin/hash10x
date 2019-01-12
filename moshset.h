/*  File: moshset.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Jan 12 12:55 2019 (rd109)
 * Created: Tue Nov  6 17:31:35 2018 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "seqhash.h"

/* object to hold sets of moshes */
typedef struct {
  Seqhash *hasher ;
  int tableBits ;		/* max 34 so size < 2^32 so index is 32bit */
  U32 size ;			/* size of value array, and other arrays over this set */
  U64 tableSize ; 		/* = 1 << tableBits */
  U64 tableMask ;		/* = tableSize - 1 */
  U32 *index ;			/* this is the primary table - size tableSize */
  U64 *value ;			/* the hashed values */
  U16 *depth ;			/* depth at each index */
  U8  *info ;			/* bits for various things */
  U32 max ;			/* number of entries in the set - must be less than size */
} Moshset ;

Moshset *moshsetCreate (Seqhash *sh, int bits, U32 size) ;
void moshsetDestroy (Moshset *ms) ; 
void moshsetWrite (Moshset *ms, FILE *f) ;
Moshset *moshsetRead (FILE *f) ;

/* this is the key low level function, both to insert new hashes and find existing ones */
U32 moshsetIndexFind (Moshset *ms, U64 hash, int isAdd) ;

/* the following act on the whole set */
void moshsetSummary (Moshset *ms, FILE *f) ;
BOOL moshsetPack (Moshset *ms)	; /* reduce size to max+1 and compress value; TRUE if changes */
void moshsetDepthPrune (Moshset *ms, int min, int max) ;
BOOL moshsetMerge (Moshset *ms1, Moshset *ms2) ;

/* info fields */
static inline void msSetCopy0 (Moshset *ms, U32 i) { ms->info[i] &= 0xfc ; }
static inline void msSetCopy1 (Moshset *ms, U32 i) { ms->info[i] = (ms->info[i] & 0xfc) | 1 ; }
static inline void msSetCopy2 (Moshset *ms, U32 i) { ms->info[i] = (ms->info[i] & 0xfc) | 2 ; }
static inline void msSetCopyM (Moshset *ms, U32 i) { ms->info[i] |= 3 ; }
static inline BOOL msIsCopy0  (Moshset *ms, U32 i) { return ((ms->info[i] & 3) == 0) ; }
static inline BOOL msIsCopy1  (Moshset *ms, U32 i) { return ((ms->info[i] & 3) == 1) ; }
static inline BOOL msIsCopy2  (Moshset *ms, U32 i) { return ((ms->info[i] & 3) == 2) ; }
static inline BOOL msIsCopyM  (Moshset *ms, U32 i) { return ((ms->info[i] & 3) == 3) ; }
static inline int  msCopy (Moshset *ms, U32 i) { return (ms->info[i] & 3) ; }

/*************************/
