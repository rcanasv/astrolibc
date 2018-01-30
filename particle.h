/*
 *  \file particle.h
 *  \brief  header file for particle properties structure
 *
 */

#ifndef PARTICLE_H
#define PARTICLE_H


#include "base.h"


typedef struct Particle
{
  float   Pos[3];
  float   Vel[3];
  float   Mass;
  int     Id;
  int     Level;
  float   Age;
  float   Metal;
  float   Chem;
  int     Type;
  float   Radius;
} Particle;


int   Particle_rad_compare (const void * a, const void * b);
void  Particle_copy        (Particle * src, Particle * dst);
void  Particle_get_radius  (Particle * P);


#endif    /*  PARTICLE_H  */
