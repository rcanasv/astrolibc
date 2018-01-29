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


#define RMSSSKIP   dummy=0; fread (&dummy, sizeof(int), 1, f);


void ramses_load_particles (Simulation * ramses, int filenum, Particle ** part);




double ** ramses_pos;
double ** ramses_vel;
double ** ramses_met;
double  * ramses_age;
double  * ramses_mass;
int     * ramses_id;
int     * ramses_lvl;




int     ndim;
int     gal_number;
int     gal_level;
double  gal_mass;
double  gal_pos[3];
double  gal_vel[3];
double  gal_ang[3];
int     nlist;


void read_gal_file (char * filename);
void read_ramses_snapshot(char * dir, char * snapshot, int numfile);

void free_ramses_arrays   (void);


double friedman(double Omega0, double OmegaL, double OmegaK, double alpha, double axp_min, double ** axp_out, double ** hexp_out, double ** tau_out, double ** t_out, int ntable);
double dadtau (double axp_tau, double Omega0, double OmegaL, double OmegaK);
double dadt (double axp_t, double Omega0, double OmegaL, double OmegaK);

#endif    /*  RAMSES_H  */
