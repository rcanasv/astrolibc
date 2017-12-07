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
  float   Aexp;
  float   Redshift;
  float   OmegaM;
  float   OmegaL;
  float   OmegaB;
  float   OmegaK;
  float   AgeUniv;
} Cosmology;


#endif    /*  COSMOLOGY_H  */
