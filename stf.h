/*
 *  \file stf.h
 *  \brief  This file contains all VELOCIraptor-stf related stuff
 *
 */

#ifndef STF_H
#define STF_H


#include "base.h"
#include "particle.h"
#include "structure.h"


typedef struct stfOutput
{
  char           prefix [NAME_LENGTH];
  int            nstruct;
  int            nprocs;
  int            iprops;
  int            iparts;
  objProps     * strctProps;
  pdata_s     ** strctParts;
} stfOutput;


#endif    /*  FINDERS_H  */
