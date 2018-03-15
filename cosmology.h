/*
 *  \file cosmology.h
 *  \brief  header file for particle properties structure
 *
 */

#ifndef COSMOLOGY_H
#define COSMOLOGY_H


#include "base.h"


typedef struct Cosmology
{
  double   HubbleParam;
  double   OmegaM;
  double   OmegaL;
  double   OmegaB;
  double   OmegaK;
} Cosmology;


#endif    /*  COSMOLOGY_H  */
