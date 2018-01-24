/*
 *  \file particle.c
 *  \brief This file contains functions related to struct Particle.
 *
 */

#include "particle.h"



int rad_compare (const void * a, const void * b)
{
  Particle * Part1 = (Particle *) a;
  Particle * Part2 = (Particle *) b;

  if (Part1->Radius > Part2->Radius)
    return 1;
  if (Part1->Radius == Part2->Radius)
    return 0;
  if (Part1->Radius < Part2->Radius)
    return -1;
}
