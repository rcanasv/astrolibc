/*
 *  \file   structure.c
 *  \brief  This file contains structure
 *
 *
 */

#include "base.h"
#include "structure.h"
#include "particle.h"


void Structure_correct_periodicity (Structure * strct, Simulation * sim)
{
  int     i, j, k;
  double  Lbox;
  double  hbox;

  Lbox = sim->Lbox;
  hbox = Lbox / 2.0;

  for (i = 0; i < strct->NumPart; i++)
  {
    strct->Part[i].Pos[0] -= strct->Pos[0];
    strct->Part[i].Pos[1] -= strct->Pos[1];
    strct->Part[i].Pos[2] -= strct->Pos[2];

    while (strct->Part[i].Pos[0] >  hbox)   strct->Part[i].Pos[0] -= Lbox;
    while (strct->Part[i].Pos[1] >  hbox)   strct->Part[i].Pos[1] -= Lbox;
    while (strct->Part[i].Pos[2] >  hbox)   strct->Part[i].Pos[2] -= Lbox;

    while (strct->Part[i].Pos[0] < -hbox)   strct->Part[i].Pos[0] += Lbox;
    while (strct->Part[i].Pos[1] < -hbox)   strct->Part[i].Pos[1] += Lbox;
    while (strct->Part[i].Pos[2] < -hbox)   strct->Part[i].Pos[2] += Lbox;
  }
}



void Structure_calculate_centre_of_mass (Structure * strct)
{
  int i;

  strct->Pos[0]  = 0.0;
  strct->Pos[1]  = 0.0;
  strct->Pos[2]  = 0.0;
  strct->TotMass = 0.0;

  for (i = 0; 1 < strct->NumPart; i++)
  {
    strct->Pos[0]  += strct->Part[i].Pos[0] * strct->Part[i].Mass;
    strct->Pos[1]  += strct->Part[i].Pos[1] * strct->Part[i].Mass;
    strct->Pos[2]  += strct->Part[i].Pos[2] * strct->Part[i].Mass;
    strct->TotMass += strct->Part[i].Mass;
  }

  strct->Pos[0] /= strct->TotMass;
  strct->Pos[1] /= strct->TotMass;
  strct->Pos[2] /= strct->TotMass;
}



void Structure_shift_to_centre_of_mass (Structure * strct)
{
  int     i;
  double  cmpos [3];
  double  cmvel [3];
  double  totmass;

  cmpos[0]  = 0.0;
  cmpos[1]  = 0.0;
  cmpos[2]  = 0.0;

  cmvel[0]  = 0.0;
  cmvel[1]  = 0.0;
  cmvel[2]  = 0.0;

  totmass   = 0.0;

  for (i = 0; i < strct->NumPart; i++)
  {
    cmpos[0]  += strct->Part[i].Pos[0] * strct->Part[i].Mass;
    cmpos[1]  += strct->Part[i].Pos[1] * strct->Part[i].Mass;
    cmpos[2]  += strct->Part[i].Pos[2] * strct->Part[i].Mass;

    cmvel[0]  += strct->Part[i].Pos[0] * strct->Part[i].Mass;
    cmvel[1]  += strct->Part[i].Pos[1] * strct->Part[i].Mass;
    cmvel[2]  += strct->Part[i].Pos[2] * strct->Part[i].Mass;

    totmass   += strct->Part[i].Mass;
  }

  cmpos[0] /= totmass;
  cmpos[1] /= totmass;
  cmpos[2] /= totmass;

  cmvel[0] /= totmass;
  cmvel[1] /= totmass;
  cmvel[2] /= totmass;

  for (i = 0; i < strct->NumPart; i++)
  {
    strct->Part[i].Pos[0] -= cmpos[0];
    strct->Part[i].Pos[1] -= cmpos[1];
    strct->Part[i].Pos[2] -= cmpos[2];

    strct->Part[i].Vel[0] -= cmvel[0];
    strct->Part[i].Vel[1] -= cmvel[1];
    strct->Part[i].Vel[2] -= cmvel[2];
  }
}


void Structure_get_particle_radius (Structure * strct)
{
  int i;
  for (i = 0; i < strct->NumPart; i++)
    Particle_get_radius (&strct->Part[i]);
}
