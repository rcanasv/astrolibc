/*
 *  \file particle.c
 *  \brief This file contains functions related to struct Particle.
 *
 */

#include "particle.h"



int Particle_rad_compare (const void * a, const void * b)
{
  Particle * Part1 = (Particle *) a;
  Particle * Part2 = (Particle *) b;

  if (Part1->Radius >  Part2->Radius)  return  1;
  if (Part1->Radius == Part2->Radius)  return  0;
  if (Part1->Radius <  Part2->Radius)  return -1;
}


int Particle_id_compare (const void * a, const void * b)
{
  Particle * Part1 = (Particle *) a;
  Particle * Part2 = (Particle *) b;

  if (Part1->Id >  Part2->Id)  return  1;
  if (Part1->Id == Part2->Id)  return  0;
  if (Part1->Id <  Part2->Id)  return -1;
}


void Particle_copy (Particle * src, Particle * dst)
{
  memcpy (dst, src, sizeof(Particle));
  /*
  dst->Pos[0] = src->Pos[0];
  dst->Pos[1] = src->Pos[1];
  dst->Pos[2] = src->Pos[2];
  dst->Vel[0] = src->Vel[0];
  dst->Vel[1] = src->Vel[1];
  dst->Vel[2] = src->Vel[2];
  dst->Mass   = src->Mass;
  dst->Id     = src->Id;
  dst->Age    = src->Age;
  dst->Metal  = src->Metal;
  dst->Chem   = src->Chem;
  dst->Type   = src->Type;
  dst->Radius = src->Radius;
  dst->HostID = src->HostID;
  dst->DirectHostID = src->DirectHostID;
  */
}


void Particle_get_radius (Particle * P)
{
  P->Radius = sqrt (P->Pos[0]*P->Pos[0] + \
                    P->Pos[1]*P->Pos[1] + \
                    P->Pos[2]*P->Pos[2]);
}

void Particle_correct_periodicity (Particle * P, double lbox_2)
{
  while (P->Pos[0] >  lbox_2) P->Pos[0] -= 2*lbox_2;
  while (P->Pos[1] >  lbox_2) P->Pos[1] -= 2*lbox_2;
  while (P->Pos[2] >  lbox_2) P->Pos[2] -= 2*lbox_2;

  while (P->Pos[0] < -lbox_2) P->Pos[0] += 2*lbox_2;
  while (P->Pos[1] < -lbox_2) P->Pos[1] += 2*lbox_2;
  while (P->Pos[2] < -lbox_2) P->Pos[2] += 2*lbox_2;
}
