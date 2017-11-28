/*
 * \file allvars.h
 * \brief Header file with global variables declarations.
 *
 *
 */

#ifndef ALLVARS_H
#define ALLVARS_H


#include "base.h"


//typedef struct Options
//{
//  int    iverbose;
//  int    inputname [NAME_LENGTH];
//  int    outname   [NAME_LENGTH];
//  int    num;
//} Options;


gheader header1;


pdata_s ** P, * p;


int * ID;
int NumPart;
int outType;


char   buffer       [NAME_LENGTH];
char   longbuffer   [LONG_LENGTH];
char   stfprefix    [NAME_LENGTH];
char   galprefix    [NAME_LENGTH];
char   outprefix    [NAME_LENGTH];
char   directory    [NAME_LENGTH];
char   format       [NAME_LENGTH];
char   output_fname [NAME_LENGTH];


int NUMFILES;


#endif    /*  ALLVARS_H  */
