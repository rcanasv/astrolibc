/*
 *  \file halomaker.h
 *  \brief This file contains all HaloMaker stuff
 *
 */

#ifndef HALOMAKER_H
#define HALOMAKER_H


#include "base.h"
#include "particle.h"
#include "structure.h"


typedef struct hmProps
{
  int      nparts;
  int   *  partIDs;
  int      id;
  int      timestep;
  int      level;
  int      parentid;
  int      idfirstchild;
  int      nsubs;
  int      idnextgal;
  float    mass;
  float    pos[3];
  float    vel[3];
  float    ang[3];
  float    rsize;
  float    eigval[3];
  float    Ekin;
  float    Epot;
  float    Etot;
  float    lambda;
  float    vdisp;
  float    bvdisp;
  float    bmass;
  float    rvir;
  float    mvir;
  float    tvir;
  float    cs;
  float    nfw_rho;
  float    nfw_rad;
  int      nbins;
  float *  rbin;
  float *  sbin;
  float    mhalf;
  float    rhalf;
} hmProps;


typedef struct hmOutput
{
  char         tbricks     [NAME_LENGTH];
  char         galpath     [NAME_LENGTH];
  char         snappath    [NAME_LENGTH];
  char         snapprefix  [NAME_LENGTH];
  int          nparts;
  float        massres;
  float        aexp;
  float        OmegaM;
  float        AgeUniv;
  int          nsubs;
  int          nstructs;
  int          ntotal;
  hmProps    * strctProps;
  pdata_s   ** strctParts;
} hmOutput;


#endif    /*  FINDERS_H  */