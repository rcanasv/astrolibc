/*
 *  \file   ramses.c
 *  \brief  This file contains RAMSES input-output routines
 *
 *
 */

#include "base.h"
#include "ramses.h"


void ramses_init (Simulation * ramses)
{
  int     i;
  FILE  * f;
  char    fname  [NAME_LENGTH];
  char    dummys [NAME_LENGTH];
  char    buffer [NAME_LENGTH];


  sprintf (fname, "%s/info_%s.txt", ramses->archive.path, ramses->archive.prefix);
  if ((f = fopen (fname, "r")) == NULL)
  {
    printf ("Couldn't open file %s\n", fname);
    exit (0);
  }
  for (i = 0; i < 7; i++)
    fgets (buffer, 100, f);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->Lbox);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->Time);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->a);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->cosmology.HubbleParam);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->cosmology.OmegaM);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->cosmology.OmegaL);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %s ", dummys, dummys, dummys);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->cosmology.OmegaB);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->unit_l);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->unit_d);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->unit_t);
  fclose (f);

  ramses->unit_m = ramses->unit_d * ramses->unit_l * ramses->unit_l * ramses->unit_l;  // in grams
  ramses->unit_m = ramses->unit_m / 1.989e+33;                                         // in solar masses

  ramses->unit_v = ramses->unit_l / ramses->unit_t;                                    // in cm / s
  ramses->unit_v = ramses->unit_v / 100000.0;                                          // in km / s

  ramses->unit_l = ramses->unit_l / 3.08e+21;                                          // in kpc
  ramses->h = ramses->cosmology.HubbleParam / 100.0;

  // Box is now in kpc
  printf ("Lbox    %lf\n", ramses->Lbox);
  printf ("unit_l  %lf\n", ramses->unit_l);
  ramses->Lbox *= ramses->unit_l;
  printf ("Lbox    %lf\n", ramses->Lbox);
}


//
//  Read RAMSES Particle FILE
//
void ramses_load_particles (Simulation * ramses, int filenum, Particle ** part)
{

  int     i, j;
  char    fname  [NAME_LENGTH];
  char    buffer [NAME_LENGTH];

  int     dummy;
  int     dummyi;
  float   dummyf;
  double  dummyd;
  char    dummys [NAME_LENGTH];

  FILE * f;

  Particle * P;

  //
  //  Read Info file to get Simulation info
  //
  sprintf (fname, "%s/info_%s.txt", ramses->archive.path, ramses->archive.prefix);
  if ((f = fopen (fname, "r")) == NULL)
  {
    printf ("Couldn't open file %s\n", fname);
    exit (0);
  }
  for (i = 0; i < 7; i++)
    fgets (buffer, 100, f);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->Lbox);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->Time);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->a);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->cosmology.HubbleParam);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->cosmology.OmegaM);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->cosmology.OmegaL);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %s ", dummys, dummys, dummys);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->cosmology.OmegaB);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->unit_l);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->unit_d);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->unit_t);
  fclose (f);

  ramses->unit_m = ramses->unit_d * ramses->unit_l * ramses->unit_l * ramses->unit_l;  // in grams
  ramses->unit_m = ramses->unit_m / 1.989e+33;                                         // in solar masses

  ramses->unit_v = ramses->unit_l / ramses->unit_t;                                    // in cm / s
  ramses->unit_v = ramses->unit_v / 100000.0;                                          // in km / s

  ramses->unit_l = ramses->unit_l / 3.08e+21;                                          // in kpc

  // Box is now in kpc
  ramses->Lbox *= ramses->unit_l;

  // H0 -> h
  ramses->h = ramses->cosmology.HubbleParam / 100.0;

  //
  //  Read Particle file to get Simulation info
  //
  sprintf (fname, "%s/part_%s.out%05d", ramses->archive.path, ramses->archive.prefix, filenum+1);
  f = fopen(fname,"r");
  if ((f = fopen (fname, "r")) == NULL)
  {
    printf ("Couldn't open file %s\n", fname);
    exit (0);
  }

  //!--- Header
  RMSSSKIP  fread(&ramses->ncpu,     sizeof(int),    1, f);  RMSSSKIP
  RMSSSKIP  fread(&ramses->ndim,     sizeof(int),    1, f);  RMSSSKIP
  RMSSSKIP  fread(&ramses->npart,    sizeof(int),    1, f);  RMSSSKIP
  RMSSSKIP  fread(&ramses->seed[0],  sizeof(int),    4, f);  RMSSSKIP
  RMSSSKIP  fread(&ramses->nstarTot, sizeof(int),    1, f);  RMSSSKIP
  RMSSSKIP  fread(&ramses->mstarTot, sizeof(double), 1, f);  RMSSSKIP
  RMSSSKIP  fread(&ramses->mstarLst, sizeof(double), 1, f);  RMSSSKIP
  RMSSSKIP  fread(&ramses->nsink,    sizeof(int),    1, f);  RMSSSKIP

/*
   printf ("NumProcs        %d\n", ramses->ncpu);
   printf ("Num Dims        %d\n", ramses->ndim);
   printf ("Npart           %d\n", ramses->npart);
   for (i = 0; i < 4; i++)
     printf ("LocalSeed[%d]    %d\n", i, ramses->seed[i]);
   printf ("NstarTot        %d\n", ramses->nstarTot);
   printf ("Mstar_tot       %g\n", ramses->mstarTot);
   printf ("Mstar_lost      %g\n", ramses->mstarLst);
   printf ("NumSink         %d\n", ramses->nsink);
*/

  if ((*(part) = (Particle *) malloc (ramses->npart * sizeof(Particle))) == NULL)
  {
    printf ("Couldn't allocate memory for Particle array\n");
    exit(0);
  }

  P = *(part);

  //--- Pos
  for (i = 0; i < ramses->ndim; i++)
  {
    RMSSSKIP
    for (j = 0; j < ramses->npart; j++)
    {
      fread(&dummyd, sizeof(double), 1, f);
      P[j].Pos[i] = dummyd;
    }
    RMSSSKIP
  }

  //--- Vel
  for (i = 0; i < ramses->ndim; i++)
  {
    RMSSSKIP
    for (j = 0; j < ramses->npart; j++)
    {
      fread(&dummyd, sizeof(double), 1, f);
      P[j].Vel[i] = dummyd;
    }
    RMSSSKIP
  }

  //--- Mass
  RMSSSKIP
  for (j = 0; j < ramses->npart; j++)
  {
    fread(&dummyd, sizeof(double), 1, f);
    P[j].Mass = dummyd;
  }
  RMSSSKIP

  //--- Id
  RMSSSKIP
  for (j = 0; j < ramses->npart; j++)
  {
    fread(&dummyi, sizeof(int), 1, f);
    P[j].Id = dummyi;
  }
  RMSSSKIP

  //--- Level
  RMSSSKIP
  for (j = 0; j < ramses->npart; j++)
  {
    fread(&dummyi, sizeof(int), 1, f);
    P[j].Level = dummyi;
  }
  RMSSSKIP

  //--- Birth Epoch
  RMSSSKIP
  for (j = 0; j < ramses->npart; j++)
  {
    fread(&dummyd, sizeof(double), 1, f);
    P[j].Age = dummyd;
  }
  RMSSSKIP

  //--- Metallicity if ((STAR || SINK) && (METAL))
  //
  //  Skip for the moment.
  //
  //for (i = 0; i < 11; i++)
  //{
  //    SKIP  fread(&ramses_met[i][0], sizeof(double), npart, f);  SKIP
  //}
  fclose (f);

  //
  // Convert to human readable units
  //
  for (i = 0; i < ramses->npart; i++)
  {
    P[i].Pos[0] *= ramses->unit_l;
    P[i].Pos[1] *= ramses->unit_l;
    P[i].Pos[2] *= ramses->unit_l;

    P[i].Vel[0] *= ramses->unit_v;
    P[i].Vel[1] *= ramses->unit_v;
    P[i].Vel[2] *= ramses->unit_v;

    P[i].Mass   *= ramses->unit_m;
  }
}


double friedman(double Omega0, double OmegaL, double OmegaK, double alpha, double axp_min, double ** axp_out, double ** hexp_out, double ** tau_out, double ** t_out, int ntable)
{
  /*
  ! ######################################################!
  ! This subroutine assumes that axp = 1 at z = 0 (today) !
  ! and that t and tau = 0 at z = 0 (today).              !
  ! axp is the expansion factor, hexp the Hubble constant !
  ! defined as hexp=1/axp*daxp/dtau, tau the conformal    !
  ! time, and t the look-back time, both in unit of 1/H0. !
  ! alpha is the required accuracy and axp_min is the     !
  ! starting expansion factor of the look-up table.       !
  ! ntable is the required size of the look-up table.     !
  ! ######################################################!
  */
  double axp_tau = 1.0;
  double axp_t   = 1.0;
  double tau     = 0.0;
  double t       = 0.0;
  double age_tot;

  double dtau;
  double dt;
  double axp_tau_pre;
  double axp_t_pre;

  int nstep = 0;
  int nskip;
  int nout;

  *axp_out  = (double *) malloc (ntable * sizeof(double));
  *hexp_out = (double *) malloc (ntable * sizeof(double));
  *tau_out  = (double *) malloc (ntable * sizeof(double));
  *t_out    = (double *) malloc (ntable * sizeof(double));

  while ((axp_tau >= axp_min) || (axp_t >= axp_min))
  {
    nstep++;
    dtau        = alpha * axp_tau / dadtau(axp_tau, Omega0, OmegaL, OmegaK);
    axp_tau_pre = axp_tau - dadtau(axp_tau, Omega0, OmegaL, OmegaK) * dtau / 2.0;
    axp_tau     = axp_tau - dadtau(axp_tau_pre, Omega0, OmegaL, OmegaK) * dtau;
    tau         = tau - dtau;

    dt          = alpha * axp_t / dadt(axp_t, Omega0, OmegaL, OmegaK);
    axp_t_pre   = axp_t - dadt(axp_t, Omega0, OmegaL, OmegaK) * dt / 2.0;
    axp_t       = axp_t - dadt(axp_t_pre, Omega0, OmegaL, OmegaK) * dt;
    t           = t - dt;
  }

  age_tot =-t;
//   printf ("Age of the universe (in unit of 1/H0)=%e \n", -t);

  nskip = nstep / ntable;

  axp_t   = 1.0;
  axp_tau = 1.0;
  tau     = 0.0;
  t       = 0.0;

  nstep = 0;
  nout  = 0;

  t_out   [0][nout] = t;
  tau_out [0][nout] = tau;
  axp_out [0][nout] = axp_tau;
  hexp_out[0][nout] = dadtau (axp_tau, Omega0, OmegaL, OmegaK) / axp_tau;


  while ((axp_tau >= axp_min) || (axp_t >= axp_min))
  {
    nstep++;
    dtau        = alpha * axp_tau / dadtau (axp_tau, Omega0, OmegaL, OmegaK);
    axp_tau_pre = axp_tau - dadtau(axp_tau, Omega0, OmegaL, OmegaK) * dtau/2.0;
    axp_tau     = axp_tau - dadtau(axp_tau_pre, Omega0, OmegaL, OmegaK) * dtau;
    tau         = tau - dtau;

    dt          = alpha * axp_t / dadt(axp_t, Omega0, OmegaL, OmegaK);
    axp_t_pre   = axp_t - dadt(axp_t, Omega0, OmegaL, OmegaK) * dt / 2.0;
    axp_t       = axp_t - dadt(axp_t_pre, Omega0, OmegaL, OmegaK) * dt;
    t           = t -dt;

    if ((nstep%nskip) == 0)
    {
      nout = nout + 1;
      t_out   [0][nout] = t;
      tau_out [0][nout] = tau;
      axp_out [0][nout] = axp_tau;
      hexp_out[0][nout] = dadtau(axp_tau, Omega0, OmegaL, OmegaK) / axp_tau;
    }
  }

  t_out   [0][ntable-1] = t;
  tau_out [0][ntable-1] = tau;
  axp_out [0][ntable-1] = axp_tau;
  hexp_out[0][ntable-1] = dadtau(axp_tau, Omega0, OmegaL, OmegaK) / axp_tau;

  return age_tot;
}



double dadtau (double axp_tau, double Omega0, double OmegaL, double OmegaK)
{
  return sqrt(axp_tau*axp_tau*axp_tau * (Omega0 + OmegaL*axp_tau*axp_tau*axp_tau + OmegaK*axp_tau));
}



double dadt (double axp_t, double Omega0, double OmegaL, double OmegaK)
{
  return sqrt((1.0/axp_t) * (Omega0 + OmegaL*axp_t*axp_t*axp_t + OmegaK*axp_t));
}



void  ramses_structure_calculate_star_age (Simulation * ramses, Structure * strct)
{
  ;
}



void  ramses_catalog_calculate_star_age (Simulation * ramses, Catalog * ctlg)
{

  int      i, j, k;

  double * axp_frw;
  double * hexp_frw;
  double * tau_frw;
  double * t_frw;
  int      n_frw = 1000;
  double   t;

  double   time_tot;
  double   time_simu;

  FILE   * f;
  char     fname  [NAME_LENGTH];
  char     dummys [NAME_LENGTH];
  char     buffer [NAME_LENGTH];
  double   dummyd;

  Structure * strct;

  //
  //  Read Info file to get Simulation info
  //
  sprintf (fname, "%s/info_%s.txt", ramses->archive.path, ramses->archive.prefix);
  if ((f = fopen (fname, "r")) == NULL)
  {
    printf ("Couldn't open file %s\n", fname);
    exit (0);
  }
  for (i = 0; i < 7; i++)
    fgets (buffer, 100, f);
  fgets (buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &dummyd);
  fgets (buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &dummyd);
  fgets (buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &dummyd);
  fgets (buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &dummyd);
  fgets (buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &dummyd);
  fgets (buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &dummyd);
  fgets (buffer, 100, f);   sscanf (buffer, "%s %s %s ", dummys, dummys, dummys);
  fgets (buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &dummyd);
  fgets (buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &dummyd);
  fgets (buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &dummyd);
  fgets (buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &dummyd);
  fclose (f);

  //
  //  Fill frw arrays
  //
  time_tot = friedman(ramses->cosmology.OmegaM, ramses->cosmology.OmegaL, 0.0, 1e-6, 1e-3, &axp_frw, &hexp_frw, &tau_frw, &t_frw, n_frw);

  // Find neighbouring conformal time
  i = 1;
  while (tau_frw[i] > ramses->Time  && i < n_frw)
    i = i+1;

  // Interpolate time
  time_simu = t_frw[i]   * (ramses->Time - tau_frw[i-1]) / (tau_frw[i]   - tau_frw[i-1]) + \
              t_frw[i-1] * (ramses->Time - tau_frw[i])   / (tau_frw[i-1] - tau_frw[i]);



  printf ("Time simu    %lf\n", (time_tot + time_simu) / (ramses->cosmology.HubbleParam*1e5/3.08e24) / (365*24*3600*1e9));
  printf ("Hubble time  %lf\n", time_tot / (ramses->cosmology.HubbleParam*1e5/3.08e24) / (365*24*3600*1e9));
  printf ("i               %d\n", i);
  printf ("time            %e\n", ramses->Time);
  printf ("time_tot        %e\n", time_tot);
  printf ("time_simu       %e\n", time_simu);
  printf ("t_tot + t_simu  %e\n", time_tot + time_simu);


  for (i = 1; i <= ctlg->nstruct; i++)
  {
    if (ctlg->strctProps[i].iPart)
    {
      strct = &ctlg->strctProps[i];
      for (j = 0; j < strct->NumPart; j++)
      {
        if (strct->Part[j].Age != 0)
        {
          k = 1;
          while (tau_frw[k] > strct->Part[j].Age  && k < n_frw)
            k++;

          t = t_frw[k]   * (strct->Part[j].Age - tau_frw[k-1]) / (tau_frw[k]   - tau_frw[k-1]) + \
              t_frw[k-1] * (strct->Part[j].Age - tau_frw[k])   / (tau_frw[k-1] - tau_frw[k]);

          // Age in years
          strct->Part[j].Age = (time_simu - t) / (ramses->cosmology.HubbleParam*1e5/3.08e24) / (365*24*3600.0);
        }
      }
    }
    else
    {
      printf ("Somethin weird happened. Particle array for Structure has not been initialized\n");
      exit (0);
    }
  }


  free (axp_frw);
  free (hexp_frw);
  free (tau_frw);
  free (t_frw);
}
