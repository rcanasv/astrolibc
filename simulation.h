/*
 *  \file simulation.h
 *  \brief This file contains headers for astronomical simulation codes
 *         e.g. Gadget-2, RAMSES, etc.
 *
 */


#ifndef SIMULATION_H
#define SIMULATION_H


#include "base.h"
#include "archive.h"
#include "particle.h"
#include "format.h"
#include "cosmology.h"


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
  double     NpartThisFile [6];
  double     NpartTot      [6];
  char       RunLabel      [NAME_LENGTH];
} Simulation;


void Simulation_init           (Simulation * sim);
void Simulation_load_particles (Simulation * sim, int filenum, Particle ** part);


#endif    /*  SIMULATION_H  */
