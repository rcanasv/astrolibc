/*
 *  \file particle.h
 *  \brief  header file for particle properties structure
 *
 */

#ifndef PARTICLE_H
#define PARTICLE_H


#include "base.h"
#include "typedef.h"

int   Particle_rad_compare (const void * a, const void * b);

void  Particle_copy        (Particle * src, Particle * dst);
void  Particle_get_radius  (Particle * P);
void  Particle_correct_periodicity (Particle * P, double lbox_2);


#endif    /*  PARTICLE_H  */
