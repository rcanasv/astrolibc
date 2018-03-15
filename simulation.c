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


void Simulation_init (Simulation * sim)
{
  if ((!strcmp(sim->archive.format, "Gadget")) ||  \
      (!strcmp(sim->archive.format, "GADGET")) ||  \
      (!strcmp(sim->archive.format, "gadget")))
    sim->format = GADGET;
  else
  if ((!strcmp(sim->archive.format, "ramses")) ||  \
      (!strcmp(sim->archive.format, "RAMSES")) ||  \
      (!strcmp(sim->archive.format, "Ramses")))
  {
    sim->format = RAMSES;
    ramses_init (sim);
  }
  else
  if ((!strcmp(sim->archive.format, "galfile")) ||  \
      (!strcmp(sim->archive.format, "Galfile")) ||  \
      (!strcmp(sim->archive.format, "GALFILE")))
    sim->format = GALFILE;
  else
  {
    printf ("Format %s not supported\n", sim->archive.format);
    printf ("Exiting...\n");
    exit (0);
  }
}

void Simulation_load_particles (Simulation * sim, int filenum, Particle ** part)
{
  if (sim->format == RAMSES) ramses_load_particles (sim, filenum, part);
  //if (sim->format == GADGET) gadget_load_particles (sim, filenum, part);
}
