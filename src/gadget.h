/*
 *  \file gadget.h
 *  \brief  This file contains all GADGET related stuff
 *
 */

#ifndef GADGET_H
#define GADGET_H


#include "base.h"
#include "typedef.h"
#include "format.h"


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

//
//  Info block
//
typedef struct ginfo
{
  char     name[5];
  char     type[9];
  int      ndim;
  int      flag[6];
  int      read;
} ginfo;


void gadget_init (Simulation * gdt);
void gadget_load_particles (Simulation * gdt, int filenum, Particle ** part);
void gadget_write_snapshot (Particle * P, int NPartTot, gheader * header, Archive * output);

void gadget_get_header    (Simulation * gdt, FILE * f, gheader * header);
long gadget_get_property  (Simulation * gdt, FILE * f, Particle * P, char * prop);
int  gadget_get_info      (Simulation * gdt, FILE * f, ginfo ** info);
int gadget_get_noinfo     (Simulation * gdt, FILE * f, ginfo ** info, int nblocks);

int  gadget_get_npart_ThisFile (Simulation * gdt, int filenum);

#endif    /*  GADGET_H  */
