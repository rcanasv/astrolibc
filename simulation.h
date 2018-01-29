/*
 *  \file simcodes.h
 *  \brief This file contains headers for astronomical simulation codes
 *         e.g. Gadget-2, RAMSES, etc.
 *
 */


#ifndef SIMULATION_H
#define SIMULATION_H


#include "base.h"
#include "archive.h"
#include "particle.h"


#define   GADGET   1
#define   RAMSES   2


typedef struct Simulation
{
  Archive    archive;
  int        format;
  double     z;
  double     H0;
  double     OmegaM;
  double     OmegaL;
  double     OmegaB;
  double     OmegaK;
  double     AgeUniv;
  double     Lbox;
  double     Time;
  double     Aexp;
  double     unit_l;
  double     unit_d;
  double     unit_v;
  double     unit_t;
  double     unit_m;
  int        ncpu;
  int        npart;
  int        seed[4];
  int        nstarTot;
  double     mstarTot;
  double     mstarLst;
  int        nsink;
} Simulation;


void Simulation_load_particles (Simulation * sim, int filenum, Particle ** part);


#endif    /*  SIMULATION_H  */
