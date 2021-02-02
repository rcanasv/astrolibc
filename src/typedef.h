/*
 *  \file typedef.h
 *  \brief This file contains headers for astronomical simulation codes
 *         e.g. Gadget-2, RAMSES, etc.
 *
 */


#ifndef TYPEDEF_H
#define TYPEDEF_H


typedef struct Archive
{
  FILE * file;
  char   name     [NAME_LENGTH];
  char   format   [NAME_LENGTH];
  char   path     [NAME_LENGTH];
  char   prefix   [NAME_LENGTH];
  int    nfiles;
} Archive;


typedef struct Cosmology
{
  double   HubbleParam;
  double   OmegaM;
  double   OmegaL;
  double   OmegaB;
  double   OmegaK;
} Cosmology;


typedef struct Particle
{
  double      Pos[3];
  double      Vel[3];
  double      Mass;
  int         Id;
  int         Level;
  float       Age;
  float       Metal;
  float       Chem;
  int         Type;
  int         dummyi;
  int         HostID;
  int         DirectHostID;
  float       Radius;
} Particle;



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
  double     Rx;
  double     R20;
  double     R90;
  int        mbpID;
  double     mbpPos[3];
  double     mbpVel[3];
  double     sigmaPosEval[3];
  double     sigmaVelEval[3];
  double     j[4];
  double     sigma;
  double     R200;
  double     M200;
  double     M30;
  double     M100;
  double     M2R50;
  int        inR200;
 
  double     Mstar200;
  double     Mnogal200;
  int        Nstar200;
  int        Nnogal200;

  // Isolation
  char       Isolated;
  char       Looselyint;
  char       Highlyint;
  char       Central;

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

  // Flags
  char flg_PartRadius;
  char flg_SortedByRadius;
  char flg_CorrectedPeriodicity;
  char flg_ShiftedCM;

  // oTask
  int oTask;
} Structure;


typedef struct Catalog
{
  Archive        archive;
  Cosmology      cosmology;
  int            format;
  int            nstruct;
  int            nprocs;
  int            iprops;
  int            iparts;
  Structure    * strctProps;
  Particle    ** strctParts;

  //
  double         a;
  double         z;
  double         h;
  double         AgeUniv;
  double         Lbox;
  double         Time;
} Catalog;



typedef struct Simulation
{
  Archive    archive;
  Cosmology  cosmology;
  int        format;

  double     a;
  double     z;
  double     h;
  double     AgeUniv;
  double     Lbox;
  double     Time;
  double     LookBackTime;

  // For RAMSES files
  double     unit_l;
  double     unit_d;
  double     unit_v;
  double     unit_t;
  double     unit_m;
  int        ndim;
  int        ncpu;
  int        npart;
  int        seed[4];
  int        nstarTot;
  double     mstarTot;
  double     mstarLst;
  int        nsink;
  int      * npartinfile;

  // For GADGET
  double     Ez;
  int        Cooling;
  int        Feedback;
  int        IcInfo;
  int        Metals;
  int        SFR;
  int        Age;
  int        Hz;
  int        HubbleParam;
  int        NfilesPerSnapshot;
  double     MassTable     [6];
  int        NpartThisFile [6];
  int        NpartTot      [6];
  int        NpartTotHW    [6];
  char       RunLabel      [NAME_LENGTH];
} Simulation;

// Specifically for RAMSES
//
// May also be aplicable for nonRamses AMR codes
// no need to change this for the moment
//
typedef struct
{
  int     myIndex;
  int     nextIndex;
  int     prevIndex;
  double  Pos       [3];
  int     fatherIndex;
  int     nborIndex [6];
  int     sonIndex  [8];
  int     cpuMap    [8];
  int     refMap    [8];
  double  octRho    [8];
  double  octPos    [8][3];
  double  dx;
  int     okOct     [8];
} Cell;

typedef struct
{
  Cell * cell;
  int    num;
} Level;

typedef struct
{
  int     alloc_ngrid;
  int     alloc_level;
  int     nlevelmax;
  int     ncpu;
  int  ** ngrid;
  Level * level;
} Grid;




#endif    /*  TYPEDEF_H  */
