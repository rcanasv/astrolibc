/*
 *  \file   structure.c
 *  \brief  This file contains structure
 *
 *
 */


#include "structure.h"


void Structure_correct_periodicity (Structure * strct, Simulation * sim)
{
  if (!strct->flg_CorrectedPeriodicity)
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
    strct->flg_CorrectedPeriodicity = 1;
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
  if (!strct->flg_ShiftedCM)
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

      cmvel[0]  += strct->Part[i].Vel[0] * strct->Part[i].Mass;
      cmvel[1]  += strct->Part[i].Vel[1] * strct->Part[i].Mass;
      cmvel[2]  += strct->Part[i].Vel[2] * strct->Part[i].Mass;

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
    strct->flg_ShiftedCM = 1;
  }
}


void Structure_sort_by_radius (Structure * strct)
{
  if (!strct->flg_SortedByRadius)
  {
    qsort (strct->Part, strct->NumPart, sizeof(Particle), Particle_rad_compare);
    strct->flg_SortedByRadius = 1;
  }
}


void Structure_get_particle_radius (Structure * strct)
{
  int i;
  if (!strct->flg_PartRadius)
  {
    for (i = 0; i < strct->NumPart; i++)
      Particle_get_radius (&strct->Part[i]);
    strct->flg_PartRadius = 1;
  }
}


void Structure_get_particle_properties (Catalog * ctlg, Simulation * sim, int * strct_to_get)
{
  if (ctlg->format == STF)        stf_structure_get_particle_properties       (ctlg, sim, strct_to_get);
  if (ctlg->format == HALOMAKER)  halomaker_structure_get_particle_properties (ctlg, sim, strct_to_get);
}


void Structure_calculate_fmass_radius    (Catalog * ctlg, Simulation * sim, int * strct_to_get, double fraction)
{
  int     i, j, k;
  double  mtot;
  double  fmass;
  Structure * strct;

  for (i = 1; i <= ctlg->nstruct; i++)
  {
    if (strct_to_get[i])
    {
      mtot  = 0.0;
      fmass = 0.0;

      strct = &ctlg->strctProps[i];
      if (!strct->iPart)
      {
        printf ("Particles have not been loaded...Exiting\n");
        exit (0);
      }

      Structure_correct_periodicity       (strct, sim);
      //Structure_shift_to_centre_of_mass   (strct);
      Structure_get_particle_radius       (strct);
      Structure_sort_by_radius            (strct);

      // Recalculate Total Mass
      for (k = 0; k < strct->NumPart; k++)
        mtot += strct->Part[k].Mass;
      fmass = fraction * mtot;
      mtot = 0;

      for (k = 0; k < strct->NumPart; k++)
      {
        mtot += strct->Part[k].Mass;
        if (mtot >= fmass)
        {
          strct->Rx = strct->Part[k].Radius;
          break;
        }
      }
    }
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
   //surface_density[i] /= (acos(-1) * (radius[i+1]*radius[i+1]*radius[i+1] - radius[i]*radius[i]*radius[i]));

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


void Structure_calculate_spherical_density (Structure * strct, double ledge, double redge, int nbins, double deltar, double ** bins, double ** Rho)
{
  //
  // This function assumes bins of same size
  //
  int      i, j, k;
  float  * tmppos;
  int      bob;

  double * density = *(Rho);
  double * radius = *(bins);

  // Radii array
  if (radius == NULL)
  {
    radius  = (double *) malloc ((nbins + 1) * sizeof(double));

    if (deltar == 0)
    {
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
      ledge = 0;
      redge = nbins * deltar;
      for (i = 0; i <= nbins; i++)
        radius[i] = i * deltar;
    }
  }
  else
  {
    deltar = radius[1] - radius[0];
  }

  // Surface density
  if  (density == NULL)
    density = (double *) malloc (nbins * sizeof(double));
  for (i = 0; i < nbins; i++)
    density[i] = 0;

  for (i = 0; i < strct->NumPart; i++)
  {
    //printf("%e  %e  %e  %e  %e\n", strct->Part[i].Radius, strct->Part[i].Pos[0],\
    strct->Part[i].Pos[1], strct->Part[i].Pos[2], strct->Part[i].Mass);
    bob = (int) (strct->Part[i].Radius / deltar);
    if ((bob < nbins) && (bob >= 0))
      //density[bob] += 1;
      density[bob] += strct->Part[i].Mass;
  }

  for (i = 0; i < nbins; i++)
    density[i] /= (4*acos(-1)/3.0*(radius[i+1]*radius[i+1]*radius[i+1] - radius[i]*radius[i]*radius[i]));

  // Free memory
  if (*(Rho) == NULL || *(bins) == NULL)
  {
    *(Rho)  = density;
    *(bins) = radius;
  }
}


int Structure_mass_compare (const void * a, const void * b)
{
  Structure * strct1 = (Structure *) a;
  Structure * strct2 = (Structure *) b;

  if (strct1->TotMass > strct2->TotMass)
    return  1;
  if (strct1->TotMass == strct2->TotMass)
    return  0;
  if (strct1->TotMass < strct2->TotMass)
    return -1;
}


int Structure_dummyd_compare (const void * a, const void * b)
{
  Structure * strct1 = (Structure *) a;
  Structure * strct2 = (Structure *) b;

  if (strct1->dummyd > strct2->dummyd)
    return  1;
  if (strct1->dummyd == strct2->dummyd)
    return  0;
  if (strct1->dummyd < strct2->dummyd)
    return -1;
}
