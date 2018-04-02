/*
 *  \file structure.h
 *  \brief  header file containing stuff related to structures in simulations
 *  e.g.  Haloes, Subhaloes, Galaxies, Streams, etc ...
 *
 */

#ifndef STRUCTURE_H
#define STRUCTURE_H


#include "base.h"
#include "particle.h"
#include "typedef.h"
#include "format.h"
#include "stf.h"


void  Structure_correct_periodicity       (Structure * strct, Simulation * sim);
void  Structure_calculate_centre_of_mass  (Structure * strct);
void  Structure_shift_to_centre_of_mass   (Structure * strct);
void  Structure_get_particle_radius       (Structure * strct);
void  Structure_get_particle_properties   (Catalog * ctlg, Simulation * sim, int * strct_to_get);


#endif    /*  STRUCTURE_H  */
