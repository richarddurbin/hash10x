/*  File: composition.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 12 21:09 2018 (rd109)
 * Created: Sun Nov 11 17:21:40 2018 (rd109)
 *-------------------------------------------------------------------
 */

#include "seqio.h"
#include <ctype.h>

static char *typeName[] = { "FASTA", "FASTQ", "BINARY" } ;

int main (int argc, char *argv[])
{
  --argc ; ++argv ;

  timeUpdate (stdout) ;
  
  SeqIO *si = seqIOopen (*argv, 0, TRUE) ;
  if (!si) die ("failed to open sequence file %s\n", *argv) ;

  U64 totLen = 0, n = 0 ;
  U64 *totBase = new0 (256, U64) ;
  U64 *totQual = new0 (256, U64) ;
  while (seqIOread (si))
    { char *s = sqioSeq(si), *e = s + si->seqLen ;
      while (s < e) ++totBase[*s++] ;
      totLen += si->seqLen ;
      if (si->isQual)
	{ char *q = sqioQual(si), *e = q + si->seqLen ;
	  while (q < e) ++totQual[*q++] ;
	}
    }
  printf ("%s file, %lld sequences >= 0, %lld total, %.2f average\n",
	  typeName[si->type], si->nSeq, totLen, totLen / (double) si->nSeq) ;
  int i ;
  U64 totUnprint = 0 ;
  printf ("bases\n") ;
  for (i = 0 ; i < 256 ; ++i)
    if (totBase[i])
      { if (isprint(i)) printf ("  %c %lld %4.1f %%\n", i, totBase[i], totBase[i]*100.0/totLen) ;
	else totUnprint += totBase[i] ;
      }
  if (totUnprint) printf (" unprintable %lld %4.1f %%\n", totUnprint, totUnprint*100.0/totLen) ;

  if (si->isQual)
    { printf ("qualities\n") ;
      for (i = 0 ; i < 256 ; ++i)
	if (totQual[i]) printf (" %3d %lld %4.1f %%\n", i, totQual[i], totQual[i]*100.0/totLen) ;
    }
  
  seqIOclose (si) ;
  timeTotal (stdout) ;
}

/****************/