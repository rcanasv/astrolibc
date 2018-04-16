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
  float       Pos[3];
  float       Vel[3];
  float       Mass;
  int         Id;
  int         Level;
  float       Age;
  float       Metal;
  float       Chem;
  int         Type;
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
  double     R20;
  double     R90;

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


#endif    /*  TYPEDEF_H  */
