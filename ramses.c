/*
 *  \file   ramsesio.c
 *  \brief  This file contains RAMSES input-output routines
 *
 *
 */

#include "base.h"
#include "ramses.h"


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
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->cosmology.Lbox);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->Time);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->cosmology.aexp);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->cosmology.H0);
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


  //
  //  Read Particle file to get Simulation info
  //
  sprintf (fname, "%s/part_%s.out%05d", ramses->archive.path, ramses->archive.prefix, filenum+1);
  f = fopen(fname,"r");

  //!--- Header
  RMSSSKIP  fread(&ramses->ncpu,     sizeof(int),    1, f);  RMSSSKIP
  RMSSSKIP  fread(&ramses->ndim,     sizeof(int),    1, f);  RMSSSKIP
  RMSSSKIP  fread(&ramses->npart,    sizeof(int),    1, f);  RMSSSKIP
  RMSSSKIP  fread(&ramses->seed[0],  sizeof(int),    4, f);  RMSSSKIP
  RMSSSKIP  fread(&ramses->nstarTot, sizeof(int),    1, f);  RMSSSKIP
  RMSSSKIP  fread(&ramses->mstarTot, sizeof(double), 1, f);  RMSSSKIP
  RMSSSKIP  fread(&ramses->mstarLst, sizeof(double), 1, f);  RMSSSKIP
  RMSSSKIP  fread(&ramses->nsink,    sizeof(int),    1, f);  RMSSSKIP

   printf ("NumProcs        %d\n", ramses->ncpu);
   printf ("Num Dims        %d\n", ramses->ndim);
   printf ("Npart           %d\n", ramses->npart);
   for (i = 0; i < 4; i++)
     printf ("LocalSeed[%d]    %d\n", i, ramses->seed[i]);
   printf ("NstarTot        %d\n", ramses->nstarTot);
   printf ("Mstar_tot       %g\n", ramses->mstarTot);
   printf ("Mstar_lost      %g\n", ramses->mstarLst);
   printf ("NumSink         %d\n", ramses->nsink);

  if ((P = (Particle *) malloc (ramses->npart * sizeof(Particle))) == NULL)
  {
    printf ("Couldn't allocate memory for Particle array\n");
    exit(0);
  }

  //--- Pos
  for (i = 0; i < ramses->ndim; i++)
  {
    RMSSSKIP
    for (j = 0; j < ramses->npart; j++)
    {
      fread(&dummyd, sizeof(double), 1, f);
      P[j].Pos[i];
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
      P[j].Vel[i];
    }
    RMSSSKIP
  }

  //--- Mass
  RMSSSKIP
  for (j = 0; j < ramses->npart; j++)
  {
    fread(&dummyd, sizeof(double), 1, f);
    P[j].Mass;
  }
  RMSSSKIP

  //--- Id
  RMSSSKIP
  for (j = 0; j < ramses->npart; j++)
  {
    fread(&dummyi, sizeof(int), 1, f);
    P[j].Id;
  }
  RMSSSKIP

  //--- Level
  RMSSSKIP
  for (j = 0; j < ramses->npart; j++)
  {
    fread(&dummyi, sizeof(int), 1, f);
    P[j].Level;
  }
  RMSSSKIP

  //--- Birth Epoch
  RMSSSKIP
  for (j = 0; j < ramses->npart; j++)
  {
    fread(&dummyd, sizeof(double), 1, f);
    P[j].Age;
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

  *(part) = P;
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

  int    nstep   = 0;
  int nskip;
  int nout;

  *axp_out  = (double *) malloc (ntable * sizeof(double));
  *hexp_out = (double *) malloc (ntable * sizeof(double));
  *tau_out  = (double *) malloc (ntable * sizeof(double));
  *t_out    = (double *) malloc (ntable * sizeof(double));

  while ((axp_tau >= axp_min) || (axp_t >= axp_min))
  {
    nstep++;
    dtau = alpha * axp_tau / dadtau(axp_tau, Omega0, OmegaL, OmegaK);
    axp_tau_pre = axp_tau - dadtau(axp_tau, Omega0, OmegaL, OmegaK) * dtau / 2.0;
    axp_tau = axp_tau - dadtau(axp_tau_pre, Omega0, OmegaL, OmegaK) * dtau;
    tau = tau - dtau;

    dt = alpha * axp_t / dadt(axp_t, Omega0, OmegaL, OmegaK);
    axp_t_pre = axp_t - dadt(axp_t, Omega0, OmegaL, OmegaK) * dt / 2.0;
    axp_t = axp_t - dadt(axp_t_pre, Omega0, OmegaL, OmegaK) * dt;
    t = t - dt;
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
    dtau = alpha * axp_tau / dadtau (axp_tau, Omega0, OmegaL, OmegaK);
    axp_tau_pre = axp_tau - dadtau(axp_tau, Omega0, OmegaL, OmegaK) * dtau/2.0;
    axp_tau = axp_tau - dadtau(axp_tau_pre, Omega0, OmegaL, OmegaK) * dtau;
    tau = tau - dtau;

    dt = alpha * axp_t / dadt(axp_t, Omega0, OmegaL, OmegaK);
    axp_t_pre = axp_t - dadt(axp_t, Omega0, OmegaL, OmegaK) * dt / 2.0;
    axp_t = axp_t - dadt(axp_t_pre, Omega0, OmegaL, OmegaK) * dt;
    t = t -dt;

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


//   double * axp_frw;
//   double * hexp_frw;
//   double * tau_frw;
//   double * t_frw;
//   int      n_frw = 1000;
//   double   time_tot = friedman(Omega_m, Omega_l, 0.0, 1e-6, 1e-3, &axp_frw, &hexp_frw, &tau_frw, &t_frw, n_frw);
//
//   // Find neighbouring conformal time
//   i = 1;
//   while (tau_frw[i] > time  && i < n_frw)
//     i = i+1;
//
//   // Interpolate time
//   double time_simu = t_frw[i]   * (time - tau_frw[i-1]) / (tau_frw[i]   - tau_frw[i-1]) + \
//                      t_frw[i-1] * (time - tau_frw[i])   / (tau_frw[i-1] - tau_frw[i]);
//
//   printf ("Time simu    %lf\n", (time_tot + time_simu) / (H0*1e5/3.08e24) / (365*24*3600*1e9));
//   printf ("Hubble time  %lf\n", time_tot / (H0*1e5/3.08e24) / (365*24*3600*1e9));
//
//   printf ("i               %d\n", i);
//   printf ("time            %e\n", time);
//   printf ("time_tot        %e\n", time_tot);
//   printf ("time_simu       %e\n", time_simu);
//   printf ("t_tot + t_simu  %e\n", time_tot + time_simu);
//     if (ramses_age[i] != 0)
//     {
//       j = 1;
//       while (tau_frw[j] > ramses_age[i]  && j < n_frw)
//         j++;
//
//       printf ("%e  ", ramses_age[i]);
//
//       t = t_frw[j]   * (ramses_age[i] - tau_frw[j-1]) / (tau_frw[j]   - tau_frw[j-1]) + \
//           t_frw[j-1] * (ramses_age[i] - tau_frw[j])   / (tau_frw[j-1] - tau_frw[j]);
//       ramses_age[i] = (time_simu - t) / (H0*1e5/3.08e24) / (365*24*3600.0);
//       printf ("%e  %e\n", t, ramses_age[i]);
//     }
//   free (axp_frw);
//   free (hexp_frw);
//   free (tau_frw);
//   free (t_frw);
