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
  if ((!strcmp(ctlg->archive.format, "Gadget")) ||  \
      (!strcmp(ctlg->archive.format, "GADGET")) ||  \
      (!strcmp(ctlg->archive.format, "gadget")))
    ctlg->format = GADGET;
  else
    if ((!strcmp(ctlg->archive.format, "ramses")) ||  \
        (!strcmp(ctlg->archive.format, "RAMSES")) ||  \
        (!strcmp(ctlg->archive.format, "Ramses")))
      ctlg->format = RAMSES;
    else
    {
      printf ("Format %s not supported\n", ctlg->archive.format);
      printf ("Exiting...\n");
      exit (0);
    }
}

void Simulation_load_particles (Simulation * sim, int filenum, Particle ** part)
{
  if (sim->format == RAMSES) ramses_load_particles (sim, filenum, part);
  //if (sim->format == GADGET) gadget_load_particles (sim, filenum, part);
}
