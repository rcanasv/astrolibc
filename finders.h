/*
 *  \file finders.h
 *  \brief  header file for code finders structures
 *
 */

#ifndef FINDERS_H
#define FINDERS_H


typedef struct objProps
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
  int    * SubIDs;
  int      NumFiles;
  int    * FilesOfGroup;
  int      NumProg;
  int    * ProgIDs;
  double * ProgMrrts;
  int      dummy;
} objProps;


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
