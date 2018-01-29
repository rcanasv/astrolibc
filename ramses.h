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


#define   RMSSSKIP   dummy = 0; fread (&dummy, sizeof(int), 1, f);


void ramses_load_particles (Simulation * ramses, int filenum, Particle ** part);


void read_gal_file (char * filename);
void read_ramses_snapshot(char * dir, char * snapshot, int numfile);

void free_ramses_arrays   (void);


double friedman(double Omega0, double OmegaL, double OmegaK, double alpha, double axp_min, double ** axp_out, double ** hexp_out, double ** tau_out, double ** t_out, int ntable);
double dadtau (double axp_tau, double Omega0, double OmegaL, double OmegaK);
double dadt (double axp_t, double Omega0, double OmegaL, double OmegaK);

#endif    /*  RAMSES_H  */
