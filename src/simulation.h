/*
 *  \file simulation.h
 *  \brief This file contains headers for astronomical simulation codes
 *         e.g. Gadget-2, RAMSES, etc.
 *
 */


#ifndef SIMULATION_H
#define SIMULATION_H


#include "base.h"
#include "typedef.h"
#include "format.h"
#include "ramses.h"
#include "gadget.h"
#include "hdf5sim.h"


void Simulation_init               (Simulation * sim);
void Simulation_load_particles     (Simulation * sim, int filenum, Particle ** part);
void Simuation_tag_neighbour_files (Simulation * sim, int * file_list);
int  Simulation_get_npart_ThisFile (Simulation * sim, int filenum);

#endif    /*  SIMULATION_H  */
