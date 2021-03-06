/*  File: dict.h
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2003-2008
 *-------------------------------------------------------------------
 * Description: header file for DICT package, string hash tables
                developed from the corresponding functions in acedb
		Jean Thierry-Mieg and Richard Durbin 1989-
 * Exported functions:
 * HISTORY:
 * Last edited: Nov  6 12:43 2018 (rd109)
 * Created: Sat Dec 20 09:34:14 2008 (rd)
 *-------------------------------------------------------------------
 */

#ifndef DICT_DEFINED
#define DICT_DEFINED

typedef struct {
  char* *names ;
  int *table ;
  int max ;			/* current number of entries */
  int dim ;
  int size ;			/* 2^dim = size of tables */
} DICT ;

DICT *dictCreate (int size) ;
void dictDestroy (DICT *dict) ;
BOOL dictWrite (DICT *dict, FILE *f) ; /* return success or failure */
DICT *dictRead (FILE *f) ;	       /* return 0 on failure */
int dictAdd (DICT *dict, char* string, int *index) ;
int dictFind (DICT *dict, char *string, int *index) ;
char* dictName (DICT *dict, int i) ;

#define dictMax(dict)  ((dict)->max)

#endif

/*********** end of file ***********/
