/*
 *  \file structure.h
 *  \brief  header file containing stuff related to structures in simulations
 *  e.g.  Haloes, Subhaloes, Galaxies, Streams, etc ...
 *
 */

#ifndef STRUCTURE_H
#define STRUCTURE_H


#include "base.h"


typedef struct Structure
{
  int      ID;
  int      DirectHostID;
  int      HostID;
  int      NumSubs;
  int      Type;
  int      NumPart;
  double   TotMass;
  double   Pos[3];
  double   Vel[3];
  double   Efrac;
  double   Rsize;
  double   RHalfMass;
  double   Vmax;
  double   Rvmax;
  double   Vdisp;
  double   Lambda;
  double   L[3];
  double   Ix[3];
  double   Ekin;
  double   Epot;
  double   Etot;
  double   Mvir;
  double   Rvir;
  double   Tvir;
  double   Csvir;
  double   RhoNFW;
  double   ReNFW;
  int      Timestep;
  int      Level;
  int      IDfirstSub;
  int      IDnextSub;
  int    * SubIDs;
  int      NumFiles;
  int    * FilesOfGroup;
  int      NumProg;
  int    * ProgIDs;
  double * ProgMrrts;
  int      dummy;
  int      dummyi;
  double   dummyd;
  int    * PIDs;
} Structure;


#endif    /*  STRUCTURE_H  */
