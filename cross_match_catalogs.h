/*
 *  \file    cross_match_catalogs.h
 *  \brief   Header file for cross_match_catalogs.c
 *
 */

#include "astrolibc.h"
#include "allvars.h"

typedef struct Options
{
  int    iverbose;
  int    numFinders;
  char   simFormat   [NAME_LENGTH];
  char   simPath     [NAME_LENGTH];
  char   simPrefix   [NAME_LENGTH];
  int    simNumFiles;
} Options;


 int options_crossMatchCatalogs (int argc, char ** argv);
