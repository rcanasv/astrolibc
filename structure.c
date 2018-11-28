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

      strct->Part[i].Vel[0] -= strct->Vel[0];
      strct->Part[i].Vel[1] -= strct->Vel[1];
      strct->Part[i].Vel[2] -= strct->Vel[2];

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
  if (ctlg->format == STF_HDF5)   stf_structure_get_particle_properties       (ctlg, sim, strct_to_get);
  if (ctlg->format == HALOMAKER)  halomaker_structure_get_particle_properties (ctlg, sim, strct_to_get);
}



void Structure_calculate_fmass_radius (Catalog * ctlg, Simulation * sim, int * strct_to_get, double fraction)
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


void Structure_calculate_j_r (Structure * strct, double radius)
{
  int     i;
  double  mass = 0;

  Particle * P;

  Structure_get_particle_radius (strct);

  strct->j[0] = 0.0;
  strct->j[1] = 0.0;
  strct->j[2] = 0.0;

  for (i = 0; i < strct->NumPart; i++)
  {
    P = &strct->Part[i];
    if (P->Radius <= radius)
    {
      strct->j[0] += P->Mass * (P->Pos[1]*P->Vel[2] - P->Pos[2]*P->Vel[1]);
      strct->j[1] += P->Mass * (P->Pos[2]*P->Vel[0] - P->Pos[0]*P->Vel[2]);
      strct->j[2] += P->Mass * (P->Pos[0]*P->Vel[1] - P->Pos[1]*P->Vel[0]);
      mass += P->Mass;
    }
  }

  strct->j[0] /= mass;
  strct->j[1] /= mass;
  strct->j[2] /= mass;

  strct->j[3] = sqrt(strct->j[0]*strct->j[0] + \
                     strct->j[1]*strct->j[1] + \
                     strct->j[2]*strct->j[2]);
}


void Structure_calculate_sigma_v_r (Structure * strct, double radius)
{
  int     i;
  double  mass   = 0;
  double  sigmax = 0;
  double  sigmay = 0;
  double  sigmaz = 0;

  Particle * P;

  Structure_get_particle_radius (strct);
  for (i = 0; i < strct->NumPart; i++)
  {
    P = &strct->Part[i];
    if (P->Radius <= radius)
    {
      sigmax += P->Mass * P->Vel[0];
      sigmay += P->Mass * P->Vel[1];
      sigmaz += P->Mass * P->Vel[2];
      mass   += P->Mass;
    }
  }

  sigmax /= mass;
  sigmay /= mass;
  sigmaz /= mass;

  strct->sigma = sqrt(sigmax*sigmax + sigmay*sigmay + sigmaz*sigmaz);
}


void Structure_calculate_sfr (Structure * strct)
{
  int k;
  
  if (strct->Type > 7)
  {
    strct->SFR20  = 0;
    strct->SFR50  = 0;
    strct->SFR100 = 0;
    for (k = 0; k < strct->NumPart; k++)
    {
      if ((strct->Part[k].Age > 0.0))// && (strct->Part[k].Radius < opt.Rlim))
      {
        if (strct->Part[k].Age <  20e6)   strct->SFR20  += strct->Part[k].Mass;
        if (strct->Part[k].Age <  50e6)   strct->SFR50  += strct->Part[k].Mass;
        if (strct->Part[k].Age < 100e6)   strct->SFR100 += strct->Part[k].Mass;
      }
    }
    //
    // SFR in  Msun / yr
    //
    strct->SFR20  /=  20e6;
    strct->SFR50  /=  50e6;
    strct->SFR100 /= 100e6;
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



void Structure_calculate_disp_tensor_pos (Catalog * ctlg, Simulation * sim, int * strct_to_get)
{
  int     i, j, k;
  double  mtot;
  double  fmass;
  int     ndim = 3;
  double  tensor[ndim * 2];
  Structure * strct;

  gsl_matrix * disp = gsl_matrix_alloc (ndim, ndim);
  gsl_vector * eval = gsl_vector_alloc (ndim);
  gsl_matrix * evec = gsl_matrix_alloc (ndim, ndim);
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (ndim);

  for (k = 1; k <= ctlg->nstruct; k++)
  {
    if (strct_to_get[k])
    {
      strct = &ctlg->strctProps[k];
      Structure_correct_periodicity       (strct, sim);

      for (i = 0; i < ndim * 2; i++)
          tensor[i] = 0.0;

      mtot = 0;

      for (i = 0; i < strct->NumPart; i++)
      {
        tensor[0] += strct->Part[i].Mass * strct->Part[i].Pos[0] * strct->Part[i].Pos[0];
        tensor[1] += strct->Part[i].Mass * strct->Part[i].Pos[1] * strct->Part[i].Pos[1];
        tensor[2] += strct->Part[i].Mass * strct->Part[i].Pos[2] * strct->Part[i].Pos[2];

        tensor[3] += strct->Part[i].Mass * strct->Part[i].Pos[0] * strct->Part[i].Pos[1];
        tensor[4] += strct->Part[i].Mass * strct->Part[i].Pos[0] * strct->Part[i].Pos[2];
        tensor[5] += strct->Part[i].Mass * strct->Part[i].Pos[1] * strct->Part[i].Pos[2];

        mtot += strct->Part[i].Mass;
      }

      gsl_matrix_set (disp, 0, 0, tensor[0]/mtot);
      gsl_matrix_set (disp, 1, 1, tensor[1]/mtot);
      gsl_matrix_set (disp, 2, 2, tensor[2]/mtot);
      gsl_matrix_set (disp, 0, 1, tensor[3]/mtot);
      gsl_matrix_set (disp, 1, 0, tensor[3]/mtot);
      gsl_matrix_set (disp, 0, 2, tensor[4]/mtot);
      gsl_matrix_set (disp, 2, 0, tensor[4]/mtot);
      gsl_matrix_set (disp, 1, 2, tensor[5]/mtot);
      gsl_matrix_set (disp, 2, 1, tensor[5]/mtot);

      /*
      for (i = 0; i < ndim; i++)
      {
        for (j = 0; j < ndim; j++)
          printf ("%2.5f  ", gsl_matrix_get(disp, i, j));
        printf ("\n");
      }
      */

      gsl_eigen_symmv (disp, eval, evec, w);
      gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);

      /*
      for (i = 0; i < ndim; i++)
      {
        printf ("%2.5f   ", gsl_vector_get(eval, i));
        for (j = 0; j < ndim; j++)
          printf ("%2.5f  ", gsl_matrix_get(evec, j, i));
        printf ("\n");
      }
      */

      /*
      for (i = 0; i < ndim; i++)
      {
        double eval_i  = gsl_vector_get (eval, i);
        gsl_vector_view evec_i = gsl_matrix_column (evec, i);

        printf ("eigenvalue = %g\n", eval_i);
        printf ("eigenvector = \n");
        gsl_vector_fprintf (stdout, &evec_i.vector, "%g");
      }
      */

      strct->sigmaPosEval[0] = gsl_vector_get(eval, 0);
      strct->sigmaPosEval[1] = gsl_vector_get(eval, 1);
      strct->sigmaPosEval[2] = gsl_vector_get(eval, 2);

//------------------------>
/*
      for (i = 0; i < ndim * 2; i++)
          tensor[i] = 0.0;

      mtot = 0;

      double ex, ey, ez;
      double x, y, z;
      double angle, angle2, angle3;

      ex = gsl_matrix_get (evec, 0, 0);
      ey = gsl_matrix_get (evec, 1, 0);
      ez = gsl_matrix_get (evec, 2, 0);

      angle = -atan(ey/ex);
      x = ex * cos(angle) - ey * sin(angle);
      y = ey * cos(angle) + ex * sin(angle);
      z = ez;

      printf ("%2.5f  %2.5f  %2.5f\n", ex, ey, ez);
      printf ("%2.5f  %2.5f  %2.5f\n", x, y, z);

      angle2 = atan(z/x);
      ex = x * cos(angle2) + z * sin(angle2);
      ez = z * cos(angle2) - x * sin(angle2);
      ey = y;

      printf ("%2.5f  %2.5f  %2.5f\n\n", ex, ey, ez);

      //----
      ex = gsl_matrix_get (evec, 0, 1);
      ey = gsl_matrix_get (evec, 1, 1);
      ez = gsl_matrix_get (evec, 2, 1);

      x = ex * cos(angle) - ey * sin(angle);
      y = ey * cos(angle) + ex * sin(angle);
      z = ez;

      printf ("%2.5f  %2.5f  %2.5f\n", ex, ey, ez);
      printf ("%2.5f  %2.5f  %2.5f\n", x, y, z);

      ex = x * cos(angle2) + z * sin(angle2);
      ez = z * cos(angle2) - x * sin(angle2);
      ey = y;

      angle3 = -atan(ez/ey);
      y = ey * cos(angle3) - ez * sin(angle3);
      z = ey * sin(angle3) + ez * cos(angle3);
      x = ex;

      printf ("%2.5f  %2.5f  %2.5f\n\n", x, y, z);


      Structure_rotate_position_z (strct, angle);
      Structure_rotate_position_y (strct, angle2);
      Structure_rotate_position_x (strct, angle3);

      for (i = 0; i < strct->NumPart; i++)
      {
        tensor[0] += strct->Part[i].Mass * strct->Part[i].Pos[0] * strct->Part[i].Pos[0];
        tensor[1] += strct->Part[i].Mass * strct->Part[i].Pos[1] * strct->Part[i].Pos[1];
        tensor[2] += strct->Part[i].Mass * strct->Part[i].Pos[2] * strct->Part[i].Pos[2];

        tensor[3] += strct->Part[i].Mass * strct->Part[i].Pos[0] * strct->Part[i].Pos[1];
        tensor[4] += strct->Part[i].Mass * strct->Part[i].Pos[0] * strct->Part[i].Pos[2];
        tensor[5] += strct->Part[i].Mass * strct->Part[i].Pos[1] * strct->Part[i].Pos[2];

        mtot += strct->Part[i].Mass;
      }

      gsl_matrix_set (disp, 0, 0, tensor[0]/mtot);
      gsl_matrix_set (disp, 1, 1, tensor[1]/mtot);
      gsl_matrix_set (disp, 2, 2, tensor[2]/mtot);
      gsl_matrix_set (disp, 0, 1, tensor[3]/mtot);
      gsl_matrix_set (disp, 1, 0, tensor[3]/mtot);
      gsl_matrix_set (disp, 0, 2, tensor[4]/mtot);
      gsl_matrix_set (disp, 2, 0, tensor[4]/mtot);
      gsl_matrix_set (disp, 1, 2, tensor[5]/mtot);
      gsl_matrix_set (disp, 2, 1, tensor[5]/mtot);

      for (i = 0; i < ndim; i++)
      {
        for (j = 0; j < ndim; j++)
          printf ("%2.5f  ", gsl_matrix_get(disp, i, j));
        printf ("\n");
      }

      gsl_eigen_symmv (disp, eval, evec, w);
      gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);


      for (i = 0; i < ndim; i++)
      {
        printf ("%2.5f   ", gsl_vector_get(eval, i));
        for (j = 0; j < ndim; j++)
          printf ("%2.5f  ", gsl_matrix_get(evec, j, i));
        printf ("\n");
     }
*/
//--------------------------->

    }
  }
  gsl_eigen_symmv_free (w);
  gsl_vector_free (eval);
  gsl_matrix_free (evec);
  gsl_matrix_free (disp);
}


void Structure_calculate_disp_tensor_vel (Catalog * ctlg, Simulation * sim, int * strct_to_get)
{
  int     i, j, k;
  double  mtot;
  double  fmass;
  int     ndim = 3;
  double  tensor[ndim * 2];
  Structure * strct;

  gsl_matrix * disp = gsl_matrix_alloc (ndim, ndim);
  gsl_vector * eval = gsl_vector_alloc (ndim);
  gsl_matrix * evec = gsl_matrix_alloc (ndim, ndim);
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (ndim);

  for (k = 1; k <= ctlg->nstruct; k++)
  {
    if (strct_to_get[k])
    {
      strct = &ctlg->strctProps[k];
      Structure_correct_periodicity       (strct, sim);

      for (i = 0; i < ndim * 2; i++)
          tensor[i] = 0.0;

      mtot = 0.0;

      for (i = 0; i < strct->NumPart; i++)
      {
        tensor[0] += strct->Part[i].Mass * strct->Part[i].Vel[0] * strct->Part[i].Vel[0];
        tensor[1] += strct->Part[i].Mass * strct->Part[i].Vel[1] * strct->Part[i].Vel[1];
        tensor[2] += strct->Part[i].Mass * strct->Part[i].Vel[2] * strct->Part[i].Vel[2];

        tensor[3] += strct->Part[i].Mass * strct->Part[i].Vel[0] * strct->Part[i].Vel[1];
        tensor[4] += strct->Part[i].Mass * strct->Part[i].Vel[0] * strct->Part[i].Vel[2];
        tensor[5] += strct->Part[i].Mass * strct->Part[i].Vel[1] * strct->Part[i].Vel[2];

        mtot += strct->Part[i].Mass;
      }

      gsl_matrix_set (disp, 0, 0, tensor[0]/mtot);
      gsl_matrix_set (disp, 1, 1, tensor[1]/mtot);
      gsl_matrix_set (disp, 2, 2, tensor[2]/mtot);
      gsl_matrix_set (disp, 0, 1, tensor[3]/mtot);
      gsl_matrix_set (disp, 1, 0, tensor[3]/mtot);
      gsl_matrix_set (disp, 0, 2, tensor[4]/mtot);
      gsl_matrix_set (disp, 2, 0, tensor[4]/mtot);
      gsl_matrix_set (disp, 1, 2, tensor[5]/mtot);
      gsl_matrix_set (disp, 2, 1, tensor[5]/mtot);

      /*
      for (i = 0; i < ndim; i++)
      {
        for (j = 0; j < ndim; j++)
          printf ("%e  ", gsl_matrix_get(disp, i, j));
        printf ("\n");
      }
      */

      gsl_eigen_symmv (disp, eval, evec, w);
      gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);

      /*
      for (i = 0; i < ndim; i++)
      {
        for (j = 0; j < ndim; j++)
          printf ("%e  ", gsl_matrix_get(disp, i, j));
        printf ("\n");
      }

      for (i = 0; i < ndim; i++)
      {
        printf ("%e   ", gsl_vector_get(eval, i));
        for (j = 0; j < ndim; j++)
          printf ("%e  ", gsl_matrix_get(evec, i, j));
        printf ("\n");
      }
      */

      strct->sigmaVelEval[0] = gsl_vector_get(eval, 0);
      strct->sigmaVelEval[1] = gsl_vector_get(eval, 1);
      strct->sigmaVelEval[2] = gsl_vector_get(eval, 2);
    }
  }
  gsl_eigen_symmv_free (w);
  gsl_vector_free (eval);
  gsl_matrix_free (evec);
  gsl_matrix_free (disp);
}


void Structure_rotate_position_x (Structure * strct, double angle)
{
  int    i;
  float  r0;
  float  r1;
  float  r2;

  for (i = 0; i < strct->NumPart; i++)
  {
    r0 = strct->Part[i].Pos[0];
    r1 = strct->Part[i].Pos[1];
    r2 = strct->Part[i].Pos[2];

    strct->Part[i].Pos[1] =  r1 * cos(angle) - r2 * sin(angle);
    strct->Part[i].Pos[2] =  r1 * sin(angle) + r2 * cos(angle);
  }
}


void Structure_rotate_position_y (Structure * strct, double angle)
{
  int    i;
  float  r0;
  float  r1;
  float  r2;

  for (i = 0; i < strct->NumPart; i++)
  {
    r0 = strct->Part[i].Pos[0];
    r1 = strct->Part[i].Pos[1];
    r2 = strct->Part[i].Pos[2];

    strct->Part[i].Pos[0] = r0 * cos(angle) + r2 * sin(angle);
    strct->Part[i].Pos[2] = r2 * cos(angle) - r0 * sin(angle);
  }
}


void Structure_rotate_position_z (Structure * strct, double angle)
{
  int    i;
  float  r0;
  float  r1;
  float  r2;

  for (i = 0; i < strct->NumPart; i++)
  {
    r0 = strct->Part[i].Pos[0];
    r1 = strct->Part[i].Pos[1];
    r2 = strct->Part[i].Pos[2];

    strct->Part[i].Pos[0] = r0 * cos(angle) - r1 * sin(angle);
    strct->Part[i].Pos[1] = r1 * cos(angle) + r0 * sin(angle);
  }
}


void Structure_rotate_velocity_x (Structure * strct, double angle)
{
  int    i;
  float  r0;
  float  r1;
  float  r2;

  for (i = 0; i < strct->NumPart; i++)
  {
    r0 = strct->Part[i].Vel[0];
    r1 = strct->Part[i].Vel[1];
    r2 = strct->Part[i].Vel[2];

    strct->Part[i].Vel[1] =  r1 * cos(angle) - r2 * sin(angle);
    strct->Part[i].Vel[2] =  r1 * sin(angle) + r2 * cos(angle);
  }
}


void Structure_rotate_velocity_y (Structure * strct, double angle)
{
  int    i;
  float  r0;
  float  r1;
  float  r2;

  for (i = 0; i < strct->NumPart; i++)
  {
    r0 = strct->Part[i].Vel[0];
    r1 = strct->Part[i].Vel[1];
    r2 = strct->Part[i].Vel[2];

    strct->Part[i].Vel[0] = r0 * cos(angle) + r2 * sin(angle);
    strct->Part[i].Vel[2] = r2 * cos(angle) - r0 * sin(angle);
  }
}


void Structure_rotate_velocity_z (Structure * strct, double angle)
{
  int    i;
  float  r0;
  float  r1;
  float  r2;

  for (i = 0; i < strct->NumPart; i++)
  {
    r0 = strct->Part[i].Vel[0];
    r1 = strct->Part[i].Vel[1];
    r2 = strct->Part[i].Vel[2];

    strct->Part[i].Vel[0] = r0 * cos(angle) - r1 * sin(angle);
    strct->Part[i].Vel[1] = r1 * cos(angle) + r0 * sin(angle);
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
