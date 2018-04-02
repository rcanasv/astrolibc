/*
 *  \file gadget.h
 *  \brief  This file contains all GADGET related stuff
 *
 */

#ifndef GADGET_H
#define GADGET_H


#include "base.h"
#include "typedef.h"


//
//  Gadget-2  header
//
typedef struct gheader
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam;
  char     fill[256 - 6*4 - 6*8 - 2*8 - 2*4 - 6*4 - 2*4 - 4*8];      /* fills to 256 Bytes */
} gheader;


void gadget_load_particles (Simulation * gdt, int filenum, Particle ** part);
void gadget_write_snapshot (Particle * P, int NPartTot, gheader * header, Archive * output);


#endif    /*  GADGET_H  */
