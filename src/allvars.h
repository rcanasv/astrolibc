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


//gheader header1;


//pdata_s ** P, * p;


int * ID;
int NumPart;
int outType;


//char   buffer       [NAME_LENGTH];
//char   longbuffer   [LONG_LENGTH];

int NUMFILES;


// Misc Macros
#define   ALL     0

// Spherical Overdensity Macros
#define   R200C     1
#define   R200B     2
#define   R500C     3
#define   RBN98     4

#endif    /*  ALLVARS_H  */
