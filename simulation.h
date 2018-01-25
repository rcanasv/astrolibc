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
#include "cosmology.h"
#include "particle.h"


typedef struct Simulation
{
  Archive    archive;
  int        format;
  Cosmology  cosmology;
  //
  // Basically the header for Gadget files
  //
} Simulation;


void Simulation_load_particles (Simulation * sim, int filenum, Particle ** part);


#endif    /*  SIMULATION_H  */
