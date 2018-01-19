/*
 *
 *  \file    convert.h
 *  \brief   Header file for convert.c
 *
 *
 */


#ifndef CONVERT_H
#define CONVERT_H


#include "astrolibc.h"
#include "allvars.h"


#define  SIMULATION  0
#define  FINDER      1


typedef struct Options
{
  int         verbose;
  Archive     param;
  Catalog     icatalog;
  Catalog     ocatalog;
  char        type   [NAME_LENGTH];
} Options;



int options_convert (int argc, char ** argv, Options * opt);

void usage_convert  (int opt, char ** argv);

int convert_finder (Options * opt);

#endif    /*  CTLGMATCH_H  */
