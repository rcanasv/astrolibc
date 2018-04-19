/*
 *  \file   structure.c
 *  \brief  This file contains structure
 *
 *
 */


#include "structure.h"


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


void Structure_calculate_surface_density (Structure * strct, double * rotation, double ledge, double redge, int nbins, double ** bins, double ** Sigma)
{
  // This function assumes bins of same size
  //
  //
  // rotation - is assumed to be a double array with size 3
  //            containing the rotation angles for x y and z
  //            if rotation == NULL  no rotation is done
  //
  int      i, j, k;
  double   deltar;
  float  * tmppos;
  int      bob;

  double * surface_density = *(Sigma);
  double * radius = *(bins);


  // Rotate
  if (rotation != NULL)
  {
    ;
  }

  // Radii array
  if (radius == NULL)
  {
    radius  = (double *) malloc ((nbins + 1) * sizeof(double));
    // Particles are assumed to be sorted by radius
    if (ledge == 0 && redge == 0)
    {
      ledge = 0;
      redge = strct->Part[strct->NumPart-1].Radius;
    }

    if (ledge == 0)
      deltar = redge / (double) nbins;
    else
      deltar = (redge - ledge) / (double) nbins;

    for (i = 0; i <= nbins; i++)
      radius[i] = i * deltar;
  }
  else
  {
    deltar = radius[1] - radius[0];
  }

  // Copy z-coordinate
  tmppos = (float *) malloc (strct->NumPart * sizeof(float));
  for (k = 0; k < strct->NumPart; k++)
  {
    tmppos[k] = strct->Part[k].Pos[2];
    strct->Part[k].Pos[2] = 0.0;
  }

  Structure_get_particle_radius (strct);
  qsort (strct->Part, strct->NumPart, sizeof(Particle), Particle_rad_compare);

  // Surface density
  if  (surface_density == NULL)
    surface_density = (double *) malloc (nbins * sizeof(double));
  for (i = 0; i < nbins; i++)
    surface_density[i] = 0;

  for (i = 0; i < strct->NumPart; i++)
  {
    bob = (int) (strct->Part[i].Radius / deltar);
    if ((bob < nbins) && (bob >= 0))
      surface_density[bob] += strct->Part[i].Mass;
  }

  for (i = 0; i < nbins; i++)
   surface_density[i] /= (acos(-1) * (radius[i+1]*radius[i+1] - radius[i]*radius[i]));
//   surface_density[i] /= (acos(-1) * (radius[i+1]*radius[i+1]*radius[i+1] - radius[i]*radius[i]*radius[i]));

  for (i = 0; i < strct->NumPart; i++)
    strct->Part[i].Pos[2] = tmppos[i];

  // Free memory
  if (*(Sigma) == NULL || *(bins) == NULL)
  {
    *(Sigma) = surface_density;
    *(bins)  = radius;
  }
  free (tmppos);
}



void Structure_get_particle_radius (Structure * strct)
{
  int i;
  for (i = 0; i < strct->NumPart; i++)
    Particle_get_radius (&strct->Part[i]);
}



void Structure_get_particle_properties (Catalog * ctlg, Simulation * sim, int * strct_to_get)
{
  if (ctlg->format == STF)
    stf_structure_get_particle_properties (ctlg, sim, strct_to_get);
}
