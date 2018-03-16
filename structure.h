/*
 *  \file structure.h
 *  \brief  header file containing stuff related to structures in simulations
 *  e.g.  Haloes, Subhaloes, Galaxies, Streams, etc ...
 *
 */

#ifndef STRUCTURE_H
#define STRUCTURE_H


#include "base.h"
#include "typedef.h"


void  Structure_correct_periodicity       (Structure * strct, Simulation * sim);
void  Structure_calculate_centre_of_mass  (Structure * strct);
void  Structure_shift_to_centre_of_mass   (Structure * strct);
void  Structure_get_particle_radius       (Structure * strct);


#endif    /*  STRUCTURE_H  */
