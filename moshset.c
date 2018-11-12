/*  File: moshset.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description: package to handle sets of "mosh" sequence hashes
 * Exported functions:
 * HISTORY:
 * Last edited: Nov  6 21:58 2018 (rd109)
 * Created: Tue Nov  6 17:31:14 2018 (rd109)
 *-------------------------------------------------------------------
 */

#include "moshset.h"
  
Moshset *moshsetCreate (Seqhash *sh, int bits, U32 size)
{
  if (bits < 20 || bits > 34) die ("table bits %d must be between 20 and 34", bits) ;
  Moshset *ms = new0 (1, Moshset) ;
  ms->hasher = sh ;
  ms->tableBits = bits ;
  ms->tableSize = 1 << ms->tableBits ;
  ms->tableMask = ms->tableSize - 1 ;
  ms->index = new0 (ms->tableSize, U32) ;
  if (size >= (ms->tableSize >> 2)) die ("Moshset size %lld is too big for %d bits", size, bits) ;
  else if (size) ms->size = size ;
  else ms->size = (ms->tableSize >> 2) - 1 ;
  ms->value = new (ms->size, U64) ;
  ms->info = new0 (ms->size, U8) ;
  return ms ;
}

void moshsetDestroy (Moshset *ms)
{ free (ms->index) ; free (ms->value) ; free (ms->info) ; free (ms) ; }

void moshsetPack (Moshset *ms)	/* compress per-item arrays */
{ ms->size = ms->max+1 ;
  resize (ms->value, ms->size, U64) ;
  resize (ms->info, ms->size, U8) ;
}

U32 moshsetIndexFind (Moshset *ms, U64 hash, int isAdd)
{
  U64 diff = 0 ;
  U64 offset = hash & ms->tableMask ;
  U32 index = ms->index[offset] ;
  while (index && (ms->value[index] != hash))
    { if (!diff) diff = ((hash >> ms->tableBits) & ms->tableMask) | 1 ; /* odd so comprime */
      offset = (offset + diff) & ms->tableMask ;
      index = ms->index[offset] ;
    }
  if (!index && isAdd)
    { index = ms->index[offset] = ++ms->max ;
      if (ms->max >= ms->size) die ("hashTableSize is too small") ;
      ms->value[index] = hash ;
    }
  return index ;
}

void moshsetDepthPrune (Moshset *ms, U32 *depth, int min, int max)
{
  U32 i ;
  U32 N = ms->max ; ms->max = 0 ;
  memset (ms->index, 0, ms->tableSize*sizeof(U32)) ;
  for (i = 1 ; i <= N ; ++i)	/* NB index runs from 1..max */
    if (depth[i] >= min && (!max || depth[i] <= max))
      { moshsetIndexFind (ms, ms->value[i], TRUE) ;
	depth[ms->max] = depth[i] ;
	ms->info[ms->max] = ms->info[i] ;
      }
  fprintf (stderr, "  pruned Moshset from %d to %d with min %d max %d\n", N, ms->max, min, max) ;
}

void moshsetWrite (Moshset *ms, FILE *f)
{ if (fwrite ("MSHSTv1",8,1,f) != 1) die ("failed to write moshset header") ;
  if (fwrite (&ms->tableBits,sizeof(int),1,f) != 1) die ("failed to write bits") ;
  U32 size = ms->max+1 ; if (fwrite (&size,sizeof(U32),1,f) != 1) die ("failed to write size") ;
  seqhashWrite (ms->hasher, f) ;
  if (fwrite (ms->index,sizeof(U32),ms->tableSize,f) != ms->tableSize) die ("fail write index") ;
  if (fwrite (ms->value,sizeof(U64),ms->max+1,f) != ms->max+1) die ("failed to write value") ;
  if (fwrite (ms->info,sizeof(U8),ms->max+1,f) != ms->max+1) die ("failed to write info") ;
}

Moshset *moshsetRead (FILE *f)
{ char name[8] ;
  if (fread (name,8,1,f) != 1) die ("failed to read moshset header") ;
  if (strcmp (name, "MSHSTv1")) die ("bad reference header") ;
  int bits ; if (fread (&bits,sizeof(int),1,f) != 1) die ("failed to read bits") ;
  U32 size ; if (fread (&size,sizeof(U32),1,f) != 1) die ("failed to read size") ;
  Seqhash *sh = seqhashRead (f) ;
  Moshset *ms = moshsetCreate (sh, bits, size) ;
  if (fread (ms->index,sizeof(U32),ms->tableSize,f) != ms->tableSize) die ("failed read index") ;
  if (fread (ms->value,sizeof(U64),size,f) != size) die ("failed to read value") ;
  if (fread (ms->info,sizeof(U8),size,f) != size) die ("failed to read info") ;
  ms->max = size - 1 ;
  return ms ;
}

/***************************************************/
