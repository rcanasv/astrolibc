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
  int        ID;
  int        DirectHostID;
  int        HostID;
  int        NumSubs;
  int        Type;
  int        NumPart;
  double     TotMass;
  double     Pos[3];
  double     Vel[3];
  double     Efrac;
  double     Rsize;
  double     RHalfMass;
  double     Vmax;
  double     Rvmax;
  double     Vdisp;
  double     Lambda;
  double     L[3];
  double     Ix[3];
  double     Ekin;
  double     Epot;
  double     Etot;
  double     Mvir;
  double     Rvir;
  double     Tvir;
  double     Csvir;
  double     RhoNFW;
  double     ReNFW;
  double     SFR20;
  double     SFR50;
  double     SFR100;
  double     R20;
  double     R50;
  double     R100;

  // Additioanl data
  int        Timestep;
  int        Level;
  int        loaded;

  // structure particle's Info
  int        iPart;
  Particle * Part;

  // Particle IDs
  int        iIDs;
  int      * PIDs;

  // Substructure Info
  int        IDfirstSub;
  int        IDnextSub;
  int        iSubs;
  int      * SubIDs;

  // Files over which structure is distributed
  int        NumFiles;
  int        iFiles;
  int      * FilesOfGroup;

  // Matching structure/parent structure
  int        NumMatch;
  int        iMatch;
  int      * MatchIDs;
  float    * MatchMrrts;

  // Dummy variables
  int        dummy;
  int        dummyi;
  double     dummyd;
} Structure;


#endif    /*  STRUCTURE_H  */
