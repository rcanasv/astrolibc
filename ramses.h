/*
 *  \file ramses.h
 *  \brief This file all RAMSES related stuff
 *
 */


#ifndef RAMSES_H
#define RAMSES_H


#include "base.h"
#include "format.h"
#include "typedef.h"


#define   RMSSSKIP   dummy = 0; fread (&dummy, sizeof(int), 1, f);


void  ramses_load_particles                (Simulation * ramses, int filenum, Particle ** part);
void  ramses_structure_calculate_star_age  (Simulation * ramses, Structure * strct);
void  ramses_catalog_calculate_star_age    (Simulation * ramses, Catalog * ctlg);
void  ramses_init                          (Simulation * ramses);


double friedman (double Omega0,  double OmegaL, double OmegaK, double alpha, double axp_min, double ** axp_out, double ** hexp_out, double ** tau_out, double ** t_out, int ntable);
double dadtau   (double axp_tau, double Omega0, double OmegaL, double OmegaK);
double dadt     (double axp_t,   double Omega0, double OmegaL, double OmegaK);

#endif    /*  RAMSES_H  */
