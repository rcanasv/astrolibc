/*
 *  \file simulation.c
 *  \brief This file contains functions related to simulations actions
 *
 */

#include "particle.h"
#include "structure.h"
#include "catalog.h"
#include "simulation.h"



void Simulation_load_particles (Simulation * sim, int filenum, Particle ** part)
{
  if (sim->format == RAMSES) ramses_load_particles (sim, i, part);
  if (sim->format == GADGET) gadget_load_particles (sim, i, part);
}
