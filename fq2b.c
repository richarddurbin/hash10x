/*  File: fq2b.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2018
 *-------------------------------------------------------------------
 * Description: converts fastq to binary packed sequence. Can correct 10x barcodes.
 * Exported functions:
 * HISTORY:
 * Last edited: Mar 24 22:40 2018 (rd)
 * Created: Sun Mar 11 09:13:43 2018 (rd)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "readseq.h"
#include <zlib.h>

BOOL gzReadFastq (gzFile f, char *idBuf, int *idLen, char *sBuf, char *qBuf, int *sLen) ;

/********** some routines to pack sequence and qualities **************/

static BOOL isPackInitialised = 0 ;
static U32 sPack[256], qPack[256] ;
void packInit (void)
{
  if (isPackInitialised) return ;
  int i ;
  static int sIndex[] = {'a','c','g','t','A','C','G','T'} ;
  for (i = 0 ; i < 8 ; ++i) sPack[sIndex[i]] = i%4 ;
  for (i = '$' + 20 ; i < 256 ; ++i) qPack[i] = 1 ;
  isPackInitialised = TRUE ;
}

void seqPack (char *s, U32 *u, int len) /* compress s into (len+)/4 U32  */
{
  assert (isPackInitialised) ;
  int i ;
  while (len > 16)
    { *u = 0 ; for (i = 0 ; i < 16 ; ++i) *u = (*u << 2) | sPack[*s++] ;
      len -= 16 ; ++u ;
    }
  *u = 0 ; for (i = 0 ; i < len ; ++i) *u = (*u << 2) | sPack[*s++] ;
}

char *uSeq (U32 u)
{
  static char s[17], map[4] ;
  if (!*map) { map[0] = 'A' ; map[1] = 'C' ; map[2] = 'G' ; map[3] = 'T' ; }
  int i ; for (i = 0 ; i < 16 ; ++ i) s[15-i] = map[(u >> (2*i)) & 3] ;
  return s ;
}

void qualPack (char *q, U32 *u, int len)
{
  assert (isPackInitialised) ;
  int i ;
  while (len > 32)
    { *u = 0 ; for (i = 0 ; i < 32 ; ++i) *u = (*u << 1) | qPack[*q++] ;
      len -= 32 ; ++u ;
    }
  *u = 0 ; for (i = 0 ; i < len ; ++i) *u = (*u << 1) | qPack[*q++] ;
}

/**************** code to correct 10x barcodes on packed data ****************/

U8 *table ;
U32 mask[64], bits[64] ;

static inline U32 switchBase (U32 u, int i) /* map (i-1)/4th base to (i-1)%4 */
{ --i ; u &= mask[i] ; u |= bits[i] ; return u ; }

void read10xWhitelist (char *fileName)
{ FILE *f = fopen (fileName, "r") ;
  if (!f) die ("failed to open 10x whitelist file %s\n", fileName) ;
  char s[32] ;
  U32 u ; 
  int i,j,n = 0 ;

  if (!(table = calloc ((size_t)1<<32, 1))) die ("can't allocate barcode table") ;
  for (i = 0 ; i < 16 ; ++i)
    for (j = 0 ; j < 4 ; ++j)
      { mask[i*4 + j] = ~(3 << (2*i)) ; bits[i*4 + j] = j << (2*i) ; }
  
  while (!feof (f) && fscanf (f, "%s\n", s))
    { ++n ;
      if (strlen(s) != 16) die ("bad barcode line % in %s: %s", n, fileName, s) ;
      seqPack (s, &u, 16) ;
      for (i = 0 ; i < 16 ; ++i)
	{ U8 ui = 1 + i*4 + ((u >> (2*i)) & 0x3) ; /* code to reset i'th base to u[i] */
	  for (j = 0 ; j < 4 ; ++j) table[switchBase(u,1+i*4+j)] = ui ;
	}
    }
  fclose (f) ;
  fprintf (stderr, "read %d barcodes from file %s\n", n, fileName) ;
}

static int nBad = 0, nFixed = 0, nFixBase[16] ;

BOOL find10xBarcode (U32* u)
{
  if (!table[*u]) { ++nBad ; return FALSE ; }
  U32 v = switchBase (*u, table[*u]) ;
  if (v != *u) { ++nFixed ; ++nFixBase[15 - (table[*u]-1)/4] ; *u = v ; }
  return TRUE ;
}

/***************** now the main function *********************/

int main (int argc, char *argv[])
{
  BOOL isCheckId = TRUE ;
  FILE *fout = stdout ;

  packInit () ;
  
  --argc ; ++argv ;		/* eat program name */
  while (argc > 2 && *argv[0] == '-')
    if (!strcmp (*argv, "-10x")) { read10xWhitelist (argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-checkId")) { isCheckId = TRUE ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-o"))
      { if (!(fout = fopen (argv[1], "wb"))) die ("failed to open output file %s", argv[1]) ;
	argc -= 2 ; argv += 2 ;
      }
    else die ("Unknown arg %s for fq2b - run without args for usage", *argv) ;
    
  if (argc < 1 || argc > 2)
    die ("Usage: fq2fqb [opts] <fastq.gz> [<fastq.gz>]\n"
	 "  Converts fastq to binary with 2 bits per base, converting N to A (!).\n"
	 "  Throws out read names and qualities. [TODO handle qualities]\n"
	 "  If two fastq files are given they are interleaved.\n"
	 "Opts: -10x <whitelist file>\n"
	 "      -checkId  checks whether id lines match in first and second files\n"
	 "      -o <outfile> [standard output]\n"
	 "  10x option matches first16bp barcode of read 1 to whitelist.\n"
	 "  Only outputs an entry if there is a match after correcting for 1 mismatch\n") ;
  
  gzFile f1 = gzopen (argv[0], "r") ; if (!f1) die ("failed to open %s", argv[0]) ;
  gzFile f2 = 0 ;
  if (argc == 2)
    { f2 = gzopen (argv[1], "r") ; if (!f1) die ("failed to open %s", argv[0]) ; }

  char id1Buf[1024], id2Buf[1024], s1Buf[1024], s2Buf[1024], q1Buf[1024], q2Buf[1024] ;
  U32 u1s[256], u2s[256], u1q[128], u2q[128] ;
  int idLen = 0, s1Len = 0, s2Len = 0 ;
  int n = 0 ;
  while (gzReadFastq (f1, id1Buf, &idLen, s1Buf, q1Buf, &s1Len))
    { seqPack (s1Buf, u1s, s1Len) ; qualPack (q1Buf, u1q, s1Len) ;
      if (f2)
	{ if (!gzReadFastq (f2, id2Buf, &idLen, s2Buf, q2Buf, &s2Len))
	    die ("second fastq file terminated early at %d", n) ;
	  assert (!strcmp (id1Buf, id2Buf)) ;
	}
      if (table && !find10xBarcode (u1s)) continue ;
      if (f2) { seqPack (s2Buf, u2s, s2Len) ; qualPack (q2Buf, u2q, s2Len) ; }
      fwrite (u1s, 4, (s1Len+15)/16, fout) ; fwrite (u1q, 4, (s1Len+31)/32, fout) ;
      if (f2) { fwrite (u2s,4,(s2Len+15)/16,fout) ; fwrite (u2q,4,(s2Len+31)/32,fout) ;}
      ++n ;
    }
  if (f2)
    fprintf (stderr, "written %d read pairs %d + %d bp packed in %d word records\n",
	     n, s1Len, s2Len, (s1Len+15)/16 + (s1Len+31)/32 + (s2Len+15)/16 + (s2Len+31)/32) ;
  else
     fprintf (stderr, "written %d reads %d bp packed in %d word records\n",
	     n, s1Len, (s1Len+15)/16 + (s1Len+31)/32) ;
  if (table)
    { fprintf (stderr, "%d (%.1f%%) not matching barcodes were dropped\n",
	       nBad, 100.0 * nBad / (double)(nBad + n)) ;
      fprintf (stderr, "%d (%.1f%%) of those that matched were error corrected\n",
	       nFixed, 100.0 * nFixed/(double)n) ;
      fprintf (stderr, "by base position:") ;
      int i ; for (i = 0 ; i < 16 ; ++i) fprintf (stderr, " %d", nFixBase[i]) ;
      fprintf (stderr, "\n") ;
    }
}

BOOL gzReadFastq (gzFile f, char *idBuf, int *idLen, char *sBuf, char *qBuf, int *sLen)
{
  static char *buf = 0 ;
  static int nEntry = 0 ;
  char *s ;

  ++nEntry ;
  *idLen = 0 ;	/* id lines are variable length, so need to read char by char */
  for (s = idBuf ; !gzeof(f) && (*s++ = gzgetc(f)) != '\n' ; ++*idLen) ;
  if (gzeof (f)) return FALSE ;
  if (*idBuf != '@') die ("fastq id line for entry %d does not start with @", nEntry) ;
  if (idBuf[*idLen] != '\n') die ("fastq entry %d id line does not end in \\n", nEntry) ;
  idBuf[*idLen] = 0 ;

  if (*sLen)
    { if (gzread (f, sBuf, *sLen+1) != *sLen+1) die ("bad seq gzread entry %d", nEntry) ; }
  else
    for (s = sBuf ; !gzeof(f) && (*s++ = gzgetc(f)) != '\n' ; ++*sLen) ;
  if (sBuf[*sLen] != '\n') die ("fastq entry %d seq line does not end in \\n", nEntry) ;
  sBuf[*sLen] = 0 ;

  if (gzgetc(f) != '+' || gzgetc(f) != '\n') die ("bad + fastq line entry %d", nEntry) ;

  if (gzread (f, qBuf, *sLen+1) != *sLen+1) die ("bad qual gzread entry %d", nEntry) ;
  if (qBuf[*sLen] != '\n') die ("fastq entry %d qual line does not end in \\n", nEntry) ;
  qBuf[*sLen] = 0 ;
  
  return TRUE ;
}

/******** end of file ********/
