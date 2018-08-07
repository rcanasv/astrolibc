/*
 *  \file simulation.c
 *  \brief This file contains functions related to simulations actions
 *
 */


#include "simulation.h"


void Simulation_init (Simulation * sim)
{
  if ((!strcmp(sim->archive.format, "Gadget")) ||  \
      (!strcmp(sim->archive.format, "GADGET")) ||  \
      (!strcmp(sim->archive.format, "gadget")))
  {
    sim->format = GADGET;
  }
  else
  if ((!strcmp(sim->archive.format, "ramses")) ||  \
      (!strcmp(sim->archive.format, "RAMSES")) ||  \
      (!strcmp(sim->archive.format, "Ramses")))
  {
    sim->format = RAMSES;
    ramses_init (sim);
  }
  else
  if ((!strcmp(sim->archive.format, "ramses_star")) ||  \
      (!strcmp(sim->archive.format, "RAMSES_STAR")) ||  \
      (!strcmp(sim->archive.format, "Ramses_Star")))
  {
    sim->format = RAMSES_STAR;
    ramses_init (sim);
  }
  else
  if ((!strcmp(sim->archive.format, "galfile")) ||  \
      (!strcmp(sim->archive.format, "Galfile")) ||  \
      (!strcmp(sim->archive.format, "GALFILE")))
  {
    sim->format = GALFILE;
    ramses_init (sim);
  }
  else
  if ((!strcmp(sim->archive.format, "eagle")) ||  \
      (!strcmp(sim->archive.format, "Eagle")) ||  \
      (!strcmp(sim->archive.format, "EAGLE")))
  {
    sim->format = EAGLE;
    hdf5_sim_init (sim);
  }
  else
  if ((!strcmp(sim->archive.format, "illustris")) ||  \
      (!strcmp(sim->archive.format, "Illustris")) ||  \
      (!strcmp(sim->archive.format, "ILLUSTRIS")))
  {
    sim->format = ILLUSTRIS;
  }
  else
  {
    printf ("Format %s not supported\n", sim->archive.format);
    printf ("Exiting...\n");
    exit (0);
  }
}

void Simulation_load_particles (Simulation * sim, int filenum, Particle ** part)
{
  if (sim->format == RAMSES)      ramses_load_particles   (sim, filenum, part);
  if (sim->format == RAMSES_STAR) ramses_load_particles   (sim, filenum, part);
  if (sim->format == EAGLE)       hdf5_sim_load_particles (sim, filenum, part);
  //if (sim->format == GADGET) gadget_load_particles (sim, filenum, part);
}
