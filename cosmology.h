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
  double   aexp;
  double   z;
  double   H0;
  double   OmegaM;
  double   OmegaL;
  double   OmegaB;
  double   OmegaK;
  double   AgeUniv;
  double   Lbox;
} Cosmology;


#endif    /*  COSMOLOGY_H  */
