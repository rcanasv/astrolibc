/*
 *  \file ramses.h
 *  \brief This file all RAMSES related stuff
 *
 */


#ifndef RAMSES_H
#define RAMSES_H

#include "base.h"
#include "particle.h"
#include "simulation.h"


void ramses_load_particles (Simulation * rmss, int filenum, Particle ** part);


#endif    /*  RAMSES_H  */
