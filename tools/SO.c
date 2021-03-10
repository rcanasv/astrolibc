/*
 *
 *  \file    SO.c
 *  \brief   Routine to calculate spherical overdensity
 *
 *
 */

#include "SO.h"
#include <time.h>

void get_structure_SO (Catalog * ctlg, Simulation * sim, int * tasks)
{
  // The following are expected to be already in place:
  //     Simulation_init
  //     Catalog_init
  //     Catalog_load_properties
  //     Catalog_fill_SubIDS
  //     Catalog_fill_isolated
  //
  // tasks          'boolean' array with tasks of catalog to be read.
  //                Constructed this way for easy MPI
  //
  // [NOT IMPLEMENTED YET]
  // Overdensity    'int'    specifying overdensity critreia. See allvars.h for macro definition.
  //                Supported criteria:
  //                      - 'R200C'   200x rho_critical      strct.R200c
  //                      - 'R500C'   500x rho_critical      strct.R500c
  //                      - 'R200B'   200x rho_background    strct.R200b
  //                      - 'RBN98'   Bryan & Norman 1998    strct.Rbn98
  //                      - 'ALL'     all of the above
  //
  // NOTE:     - Overdensities are calculated at sim redshift
  //           - Particles will be stored in strct.PSO array. In the case
  //             of 'ALL' all particles on the largest sphere will be returned
  //             the user must use strct.RXXX value for postprocessing
  //

  int  i, j, k, l, n, m;
  int  itasks;
  int  ifiles;
  int  dummyi;

  Particle    * P;
  Particle    * Pbuff;

  Structure   * strct1;
  Structure   * strct2;

  FILE  * f;
  char    fname  [NAME_LENGTH];
  char    buffer [NAME_LENGTH];

  int     numpart;
  int     ninrad;

  int  * files_to_read = NULL;
  int  * loaded_files  = NULL;
  int  * strct_to_get  = NULL;
  int  * npartinfile   = NULL;

  stfExtendedOutput * xtndd;
  int                 ninextended;
  int                 indx;
  int                 id;

  clock_t  start_t, end_t;

  // Useful variables
  double lbox   = sim->Lbox;
  double lbox_2 = lbox / 2.0;

  // Variables for R200 stuff
  double  pi   = acos(-1.0);
  double  fac  = 4.0 * pi / 3.0;
  double  G    = 43009.1;                          // in (kpc/M_sun)*(km/s)^2
  double  H    = sim->cosmology.HubbleParam / 1000.0;  // in (km/s)/kpc
  double  Om_m = sim->cosmology.OmegaM;
  double  Om_l = sim->cosmology.OmegaL;
  double  z    = sim->z;
  double  Om_z = Om_m*pow(1+z,3) / (Om_m*pow(1+z,3) + Om_l);
  double  x    = Om_z - 1;

  double  rhocrit  = 3.0 * H * H / (8.0 * pi * G);
  double  rho500c  = 500.0 * rhocrit;
  double  rho200c  = 200.0 * rhocrit;
  double  rhobn98  = rhocrit * (18*pi*pi + 82*x - 39*x*x);
  double  rho200b  = 200.0 * rhocrit * sim->cosmology.OmegaM;

  printf ("r500c  %e\n", rho500c);
  printf ("r200c  %e\n", rho200c);
  printf ("r200b  %e\n", rho200b);
  printf ("rbn98  %e\n", rhobn98);

  // Allocate arrays
  files_to_read = (int *) malloc ((sim->archive.nfiles) * sizeof(int));
  npartinfile   = (int *) malloc ((sim->archive.nfiles) * sizeof(int));
  strct_to_get  = (int *) malloc ((ctlg->nstruct)       * sizeof(int));

  // Initialize arrays
  for (i = 0; i < sim->archive.nfiles; i++)
  {
    files_to_read[i] = 0;
    npartinfile[i]   = 0;
  }
  for (i = 0; i < sim->archive.nfiles; i++)
    strct_to_get[i] = 0;


  // Loop over tasks
  int itask;
  for (itask = 0; itask < ctlg->archive.nfiles; itask++)
  {
    if (tasks[itask])
    {
      // Tag files that need to be opened
      Catalog_get_files_of_groups (ctlg);
      // TODO: Simuation_tag_neighbour_files (sim, files_to_read);

      // Tag files to read object of interest
      for (i = 1; i <= ctlg->nstruct; i++)
      {
        strct1 = &ctlg->strctProps[i];
        if ((strct1->oTask == itask)  && \
            (strct1->Type == 7)       && \
            (strct1->NumSubs > 0))
        {
          // dummyi contains Central ID
          // dummyd contains total FOF mass
          strct2 = &ctlg->strctProps[strct1->dummyi];
          if (strct2->dummyd >= 10) // HERE MASS IS IN 10^10 Msun units
          {
            strct_to_get[i] = 1;
            for (j = 0; j < strct1->NumFiles; j++)
              files_to_read[strct1->FilesOfGroup[j]] = 1;
          }
        }
      }

      // Get number of particles to allocate array
      for (n = 0, numpart = 0; n < sim->archive.nfiles; n++)
        if (files_to_read[n])
        {
          npartinfile[n] = Simulation_get_npart_ThisFile (sim, n);
          numpart += npartinfile[n];
        }

      // Allocate Particle array memory
      if ((P = (Particle *) malloc (numpart * sizeof(Particle))) == NULL)
      {
        printf ("couldnt allocate memory %ld", numpart*sizeof(Particle));
        exit(0);
      }

      // Load over files to read particles
      for (n = 0, m = 0; n < sim->archive.nfiles; n++)
        if (files_to_read[n])
        {
          // Load particles from simulation
          Simulation_load_particles (sim, i, &Pbuff);
          for (i = 0; i < npartinfile[n]; i++)
          {
            Pbuff[i].dummyi = 0;
            Pbuff[i].StructID = 0;
          }

          // Tag particles host ID
          xtndd = NULL;
          ninextended = 0;
          ninextended = stf_load_extended_output (ctlg, n, &xtndd);
          for (j = 0; j < ninextended; j++)
          {
            id    = xtndd[j].IdStruct;
            indx  = xtndd[j].oIndex;

            if (id > 0)
              Pbuff[indx].StructID = id;
          }
          free (xtndd);

          // Copy Particles
          for (i = 0; i < npartinfile[n]; i++)
            Particle_copy (&Pbuff[i], &P[m++]);
          free (Pbuff);
        }

      if (m != numpart)
        printf ("Something's wrong mismatch on number of particles\n");

      // Loop over centrals
      for (k = 1; k <= ctlg->nstruct; k++)
      {
        strct1 = &ctlg->strctProps[k];
        if (strct_to_get[k] && strct1->inR200 == 0 && strct1->oTask == itask)
        {
          //start_t = clock();
          // Centre of R200 is central galaxy
          strct2 = &ctlg->strctProps[strct1->dummyi];
          ninrad = 0;

          // Tag particles in vicinity
          for (j = 0; j < numpart; j++)
          {
            P[j].Pos[0] -= strct2->Pos[0];
            P[j].Pos[1] -= strct2->Pos[1];
            P[j].Pos[2] -= strct2->Pos[2];
            P[j].dummyi = 0;

            Particle_correct_periodicity (&P[j], lbox_2);

            if ((fabs(P[j].Pos[0]) < 5000.0) && \
                (fabs(P[j].Pos[1]) < 5000.0) && \
                (fabs(P[j].Pos[2]) < 5000.0))
            {
              ninrad++;
              P[j].dummyi = 1;
            }
          }

          // Allocate memory
          Pbuff = (Particle *) malloc (ninrad * sizeof(Particle));

          // Copy particles
          for (j = 0, i = 0; j < numpart; j++)
          {
            if (P[j].dummyi)
            {
              Particle_copy (&P[j], &Pbuff[i]);
              Particle_get_radius (&Pbuff[i]);
              i++;
            }
            P[j].Pos[0] += strct2->Pos[0];
            P[j].Pos[1] += strct2->Pos[1];
            P[j].Pos[2] += strct2->Pos[2];
          }
          //end_t = clock();
          //printf ("%d  took %f seconds\n", k, (end_t - start_t)/(double)CLOCKS_PER_SEC);

          //start_t = clock();
          qsort (Pbuff, ninrad, sizeof(Particle), Particle_rad_compare);
          //end_t = clock();
          //printf ("%d  took %f seconds qsort\n", k, (end_t - start_t)/(double)CLOCKS_PER_SEC);

          //start_t = clock();
          // Starts loop to get SO
          int    icheck = 0;
          double msum   = 0;
          double rad    = 0;
          double rho    = 0;

          strct1->n500c = 0;
          strct1->n200c = 0;
          strct1->n200b = 0;
          strct1->nbn98 = 0;
          for (i = 0; i < ninrad; i++)
          {
            msum += Pbuff[i].Mass;
            rad = Pbuff[i].Radius;
            rho = msum / (fac * rad * rad * rad);

            if (rho < rho500c && strct1->n500c == 0) {strct1->n500c = i;  strct1->R500c = rad;  strct1->M500c = msum;  icheck++;}
            if (rho < rho200c && strct1->n200c == 0) {strct1->n200c = i;  strct1->R200c = rad;  strct1->M200c = msum;  icheck++;}
            if (rho < rho200b && strct1->n200b == 0) {strct1->n200b = i;  strct1->R200b = rad;  strct1->M200b = msum;  icheck++;}
            if (rho < rhobn98 && strct1->nbn98 == 0) {strct1->nbn98 = i;  strct1->Rbn98 = rad;  strct1->Mbn98 = msum;  icheck++;}

            if (icheck == 4)
            {
              strct1->nSO = i;
              break;
            }
          }

          strct1->PSO = (Particle *) malloc (strct1->nSO * sizeof(Particle));
          for (i = 0; i < strct1->nSO; i++)
            Particle_copy (&Pbuff[i], &strct1->PSO[i]);
          free(Pbuff);
          /*
          for (i = 0; i < numpart; i++)
          {
            P[i].Pos[0] += strct2->Pos[0];
            P[i].Pos[1] += strct2->Pos[1];
            P[i].Pos[2] += strct2->Pos[2];
          }
          */
          //end_t = clock();
          //printf ("%d  took %f seconds\n", k, (end_t - start_t)/(double)CLOCKS_PER_SEC);
          printf ("strct  %d  %d  %e  %e  %e  %e\n", k, ninrad, strct1->M500c, strct1->M200c, strct1->M200b, strct1->Mbn98);
        } // if strct_to_get
      } // loop over structures
    } // if tasks[itasks]
  }// Loop over tasks


  // Free memory
  free (files_to_read);
  free (npartinfile);
  free (strct_to_get);
  free (P);
}
