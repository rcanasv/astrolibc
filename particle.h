/*
 *  \file particle.h
 *  \brief  header file for particle properties structure
 *
 */

#ifndef PARTICLE_H
#define PARTICLE_H


#include "base.h"


/*
typedef struct pdata_d
{
  double   Pos[3];
  double   Vel[3];
  double   Mass;
  int      Id;
  double   Age;
  double   Metal;
  double   Chem;
  int      Type;
} pdata_d;
*/

typedef struct Particle
{
  float   Pos[3];
  float   Vel[3];
  float   Mass;
  int     Id;
  float   Age;
  float   Metal;
  float   Chem;
  int     Type;
  float   Radius;
} Particle;


int Particle_rad_compare (const void * a, const void * b);
int Particle_copy        (Particle * src, Particle * dst);

#endif    /*  PARTICLE_H  */
