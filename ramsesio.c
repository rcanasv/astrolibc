/*
 *  \file   ramsesio.c
 *  \brief  This file contains RAMSES input-output routines
 *
 *
 */

#include "base.h"
#include "ramses.h"
#include "ramsesio.h"




//
//  Read RAMSES Particle FILE
//
void read_ramses_snapshot(char * dir, char * snapshot, int numfile)
{
  char name[100];
  char buffer[100];
  char dummys[100];
  int i,j;

  double boxlen;
  double time;
  double aexp;
  double H0;
  double Omega_m;
  double Omega_l;
  double Omega_b;
  double unit_l;
  double unit_d;
  double unit_v;
  double unit_t;
  double unit_m;

  sprintf (name, "%s/info_%s.txt", dir, snapshot);
  FILE * ff = fopen(name,"r");
  for (i = 0; i < 7; i++)
    fgets(buffer, 100, ff);

  fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &boxlen);
  fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &time);
  fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &aexp);
  fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &H0);
  fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &Omega_m);
  fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &Omega_l);
  fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %s ", dummys, dummys, dummys);
  fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &Omega_b);
  fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &unit_l);
  fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &unit_d);
  fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &unit_t);
  fclose (ff);

  unit_m = unit_d * unit_l * unit_l * unit_l;  // in grams
  unit_m = unit_m / 1.989e+33;                 // in solar masses

  unit_v = unit_l / unit_t;                    // in cm / s
  unit_v = unit_v / 100000.0;                  // in km / s

  unit_l = unit_l / 3.08e+21;
//   printf("DM mass  %g\n", 7.755314938e-10 * unit_m);

  sprintf (name, "%s/part_%s.out%05d", dir, snapshot, numfile+1);
//   printf ("Reading particle file   %s\n", name);

  FILE * f = fopen(name,"r");

//   fread (ptr, num, nmemb, stream)
  int dummy;
  int info[100];

#define SSKIP   dummy=0; fread (&dummy, sizeof(int), 1, f);
#define PSSKIP  printf ("%d  ", dummy);
#define PSSKIP2 printf ("%d\n", dummy);

  int lala;

  int size = 0;

  int     ncpu;
  int     npart;
  int     seed[4];
  int     nstarTot;
  double  mstarTot;
  double  mstarLst;
  int     nsink;

  //!--- Header
  SSKIP  fread(&ncpu,     sizeof(int),    1, f);  SSKIP
  SSKIP  fread(&ndim,     sizeof(int),    1, f);  SSKIP
  SSKIP  fread(&npart,    sizeof(int),    1, f);  SSKIP
  SSKIP  fread(&seed[0],  sizeof(int),    4, f);  SSKIP
  SSKIP  fread(&nstarTot, sizeof(int),    1, f);  SSKIP
  SSKIP  fread(&mstarTot, sizeof(double), 1, f);  SSKIP
  SSKIP  fread(&mstarLst, sizeof(double), 1, f);  SSKIP
  SSKIP  fread(&nsink,    sizeof(int),    1, f);  SSKIP

//   printf ("NumProcs        %d\n", ncpu);
//   printf ("Num Dims        %d\n", ndim);
//   printf ("Npart           %d\n", npart);
//   for (i = 0; i < 4; i++)
//     printf ("LocalSeed[%d]    %d\n", i, seed[i]);
//   printf ("NstarTot        %d\n", nstarTot);
//   printf ("Mstar_tot       %g\n", mstarTot);
//   printf ("Mstar_lost      %g\n", mstarLst);
//   printf ("NumSink         %d\n", nsink);

  //!--- Allocate for DM and Stars
  ramses_pos = (double **) malloc (ndim * sizeof(double *));
  ramses_vel = (double **) malloc (ndim * sizeof(double *));
  ramses_met = (double **) malloc (11   * sizeof(double *));

  for (i = 0; i < ndim; i++) ramses_pos[i] = (double *) malloc (npart * sizeof(double));
  for (i = 0; i < ndim; i++) ramses_vel[i] = (double *) malloc (npart * sizeof(double));
  for (i = 0; i < 11;   i++) ramses_met[i] = (double *) malloc (npart * sizeof(double));

  ramses_age  = (double *) malloc (npart * sizeof(double));
  ramses_mass = (double *) malloc (npart * sizeof(double));
  ramses_id   = (int    *) malloc (npart * sizeof(int));
  ramses_lvl  = (int    *) malloc (npart * sizeof(int));

  //--- Pos
  for (i = 0; i < ndim; i++)
  {
    SSKIP  fread(&ramses_pos[i][0], sizeof(double), npart, f);  SSKIP
  }
  //--- Vel
  for (i = 0; i < ndim; i++)
  {
    SSKIP  fread(&ramses_vel[i][0], sizeof(double), npart, f);  SSKIP
  }
  //--- Mass
  SSKIP  fread(&ramses_mass[0], sizeof(double), npart, f);  SSKIP

  //--- Id
  SSKIP  fread(&ramses_id[0], sizeof(int), npart, f);  SSKIP

  //--- Level
  SKIP  fread(&ramses_lvl[0], sizeof(int), npart, f);  SKIP

  //--- Birth Epoch
  SKIP  fread(&ramses_age[0], sizeof(double), npart, f);  SKIP

  //--- Metallicity if ((STAR || SINK) && (METAL))
  for (i = 0; i < 11; i++)
  {
    SKIP  fread(&ramses_met[i][0], sizeof(double), npart, f);  SKIP
  }
  fclose (f);


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

  //
  // Convert to human readable units
  //
  double t;
  for (i = 0; i < npart; i++)
  {
    ramses_pos[0][i] *= unit_l;
    ramses_pos[1][i] *= unit_l;
    ramses_pos[2][i] *= unit_l;

    ramses_vel[0][i] *= unit_v;
    ramses_vel[1][i] *= unit_v;
    ramses_vel[2][i] *= unit_v;

    ramses_mass[i]   *= unit_m;

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
  }

//   free (axp_frw);
//   free (hexp_frw);
//   free (tau_frw);
//   free (t_frw);
}



void free_ramses_arrays (void)
{
  int i;

  for (i = 0; i < ndim; i++)
  {
    free (ramses_pos[i]);
    free (ramses_vel[i]);
  }
  for (i = 0; i < 11; i++)
    free (ramses_met[i]);

  free (ramses_pos);
  free (ramses_vel);
  free (ramses_mass);
  free (ramses_id);
  free (ramses_met);
  free (ramses_age);
  free (ramses_lvl);
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
