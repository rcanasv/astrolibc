/*
 *  \file simulation.c
 *  \brief This file contains functions related to simulations actions
 *
 */

#include "particle.h"
#include "structure.h"
#include "catalog.h"
#include "simulation.h"
#include "format.h"

#include "ramses.h"


void Simulation_load_particles (Simulation * sim, int filenum, Particle ** part)
{
  if (sim->format == RAMSES) ramses_load_particles (sim, filenum, part);
  //if (sim->format == GADGET) gadget_load_particles (sim, filenum, part);
}
