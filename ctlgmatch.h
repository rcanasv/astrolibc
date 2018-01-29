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
#include "stf.h"
#include "halomaker.h"
#include "ramses.h"


typedef struct Options
{
  int            verbose;
  Archive        param;
  Archive        mtree;
  Archive        output;
  int            numCatalogs;
  Catalog     *  catalog;
  Simulation  *  simulation;
} Options;


int   ctlgMatch_options (int argc, char ** argv, Options * opt);
void  ctlgMatch_usage   (int opt, char ** argv);
void  ctlgMatch_params  (Options * opt);


#endif    /*  CTLGMATCH_H  */
