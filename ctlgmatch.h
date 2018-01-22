/*
 *
 *  \file    ctlgmatch.h
 *  \brief   Header file for cross_match_catalogs.c
 *
 *
 */


#ifndef CTLGMATCH_H
#define CTLGMATCH_H


#include "astrolibc.h"
#include "allvars.h"


typedef struct Options
{
  int         verbose;
  Archive     param;
  int         numCatalogs;
  Catalog  *  catalog;
  Archive     arx;
  char        outprefix   [NAME_LENGTH];
} Options;


int options_ctlgMatch (int argc, char ** argv, Options * opt);

void usage_ctlgMatch  (int opt, char ** argv);


#endif    /*  CTLGMATCH_H  */
