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
  Structure   * strct3;

  FILE  * f;
  char    fname  [NAME_LENGTH];
  char    buffer [NAME_LENGTH];

  int     numpart;
  int     ninrad;
  int     nbuff;

  int  * files_to_read = NULL;
  int  * loaded_files  = NULL;
  int  * strct_to_get  = NULL;
  int  * npartinfile   = NULL;

  stfExtendedOutput * xtndd;
  int                 ninextended;
  int                 indx;
  int                 id;

  int      icheck;
  double   msum;
  double   rad;
  double   rho;

  clock_t  start_t, end_t;

  double  msum_ap;
  double  msum_str;
  double  msum_dif;
  double  ms200c_str, ms200b_str, ms500c_str, msbn98_str;
  double  ms200c_dif, ms200b_dif, ms500c_dif, msbn98_dif;

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
      Pbuff = NULL;
      nbuff = 0;
      for (k = 1; k <= ctlg->nstruct; k++)
      {
        strct1 = &ctlg->strctProps[k];
        if (strct_to_get[k] && strct1->inR200 == 0 && strct1->oTask == itask && strct1->Type == 7 && strct1->NumSubs > 0)
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
          if (ninrad > nbuff)
          {
            nbuff = ninrad;
            if (Pbuff != NULL)
              free (Pbuff);
            Pbuff = (Particle *) malloc (nbuff * sizeof(Particle));
          }

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
          if (i != ninrad)  printf("Something went wrong\n");
          qsort (Pbuff, ninrad, sizeof(Particle), Particle_rad_compare);
          //end_t = clock();
          //printf ("%d  took %f seconds qsort\n", k, (end_t - start_t)/(double)CLOCKS_PER_SEC);

          //start_t = clock();
          // Starts loop to get SO
          icheck = 0;
          msum   = 0;
          rad    = 0;
          rho    = 0;

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

          msum_ap  = 0;  // for apperture acum
          msum_str = 0;  // for structure acum
          msum_dif = 0;  // for diffuse acum
          for (i = 0; i <= strct1->nSO; i++)
          {
            if (Pbuff[i].Type == 4)
            {
              strct3 = &ctlg->strctProps[Pbuff[i].StructID];

              // Spherical appertues e.g. Pillepich
              if (Pbuff[i].StructID == strct1->ID || Pbuff[i].StructID == strct2->ID)
              {
                msum_ap += Pbuff[i].Mass;
                if (Pbuff[i].Radius < 30.0)    strct1->M30  = msum_ap;
                if (Pbuff[i].Radius < 100.0)   strct1->M100 = msum_ap;
              }

              // For IHSC comp
              if (Pbuff[i].StructID > 0 && strct3->Type > 7)
              {
                msum_str += Pbuff[i].Mass;
                if (Pbuff[i].Radius < strct1->R200c)  strct1->ms200c_str = msum_str;
                if (Pbuff[i].Radius < strct1->R200b)  strct1->ms200b_str = msum_str;
                if (Pbuff[i].Radius < strct1->R500c)  strct1->ms500c_str = msum_str;
                if (Pbuff[i].Radius < strct1->Rbn98)  strct1->msbn98_str = msum_str;
              }
              else
              {
                msum_dif += Pbuff[i].Mass;
                if (Pbuff[i].Radius < strct1->R200c)  strct1->ms200c_dif = msum_dif;
                if (Pbuff[i].Radius < strct1->R200b)  strct1->ms200b_dif = msum_dif;
                if (Pbuff[i].Radius < strct1->R500c)  strct1->ms500c_dif = msum_dif;
                if (Pbuff[i].Radius < strct1->Rbn98)  strct1->msbn98_dif = msum_dif;
              }
            }
          }

          strct1->nlowres = 0;
          strct1->tlowres  = 0;
          strct1->rlowres  = 0;
  	  for (i = 0; i < strct1->n200c; i++)
          {
            if (Pbuff[i].Type == 3 || Pbuff[i].Type == 2)
            {
              if (strct1->nlowres == 0)
  	      {
                strct1->tlowres = Pbuff[i].Type;
                strct1->rlowres = Pbuff[i].Radius;
              }
              strct1->nlowres++;
            }
	  }
          //end_t = clock();
          //printf ("%d  took %f seconds\n", k, (end_t - start_t)/(double)CLOCKS_PER_SEC);
          //register double r = Pbuff[ninrad-1].Radius;
          //strct2 = &ctlg->strctProps[strct1->dummyi]; // Central
          //printf ("strct  %d  Rbuff  %e  M200c %e  Mctrl %e  M200b %e R500c %e  R200c %e  Rbn98 %e  R200b %e \n", \
                  k, r, strct1->M200c, strct2->TotMass, strct1->M200b, strct1->R500c, strct1->R200c, strct1->Rbn98, strct1->R200b);
          printf ("strct  %d  %d  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e\n", \
                  k, ninrad, strct1->M500c, strct1->M200c, strct1->M200b, strct1->Mbn98, \
                  strct1->ms200c_str, strct1->ms200b_str, strct1->msbn98_str, \
                  msum_ap, msum_dif, msum_str);
        } // if strct_to_get
      } // loop over structures

      // Free Particle buffer
      if (Pbuff != NULL)
        free(Pbuff);

      // Open file to write
      sprintf (buffer, "%s.ihsc.spho.%d", ctlg->archive.prefix, itask);
      f = fopen (buffer, "w");
      for (k = 1; k <= ctlg->nstruct; k++)
      {
        strct1 = &ctlg->strctProps[k];
        if (strct_to_get[k] && strct1->inR200 == 0 && strct1->oTask == itask && strct1->NumSubs > 0 && strct1->Type == 7)
        {
          strct2 = &ctlg->strctProps[strct1->dummyi]; // Central
          strct3 = &ctlg->strctProps[strct1->SubIDs[strct1->NumSubs-2]]; // Scnd
          fprintf (f, "%e  ", strct2->dummyd);      // Total Stellar Mass
          fprintf (f, "%e  ", strct1->TotMass);     // Mass IHSC
          fprintf (f, "%e  ", strct2->TotMass);     // Mass Central
          fprintf (f, "%e  ", strct3->TotMass);     // Mass Second most massive Gal
          fprintf (f, "%5d ", strct1->NumSubs);     // NumSubs
          fprintf (f, "%5d ", strct2->Central);     // Is Central central?
          fprintf (f, "%5d ", strct1->ID);          // ID IHSC
          fprintf (f, "%5d ", strct2->ID);          // ID Central
          fprintf (f, "%5d ", strct3->ID);          // ID Second most
          fprintf (f, "%e  ", strct1->M30);         // 30 kpc stellar mass no sats
          fprintf (f, "%e  ", strct1->M100);        // 100 kpc stellar mass no sats
          fprintf (f, "%e  ", strct1->R500c);       // R500C
          fprintf (f, "%e  ", strct1->M500c);
          fprintf (f, "%e  ", strct1->ms500c_str);
          fprintf (f, "%e  ", strct1->ms500c_dif);
          fprintf (f, "%e  ", strct1->R200c);       // R200C
          fprintf (f, "%e  ", strct1->M200c);
          fprintf (f, "%e  ", strct1->ms200c_str);
          fprintf (f, "%e  ", strct1->ms200c_dif);
          fprintf (f, "%e  ", strct1->R200b);       // R200B
          fprintf (f, "%e  ", strct1->M200b);
          fprintf (f, "%e  ", strct1->ms200b_str);
          fprintf (f, "%e  ", strct1->ms200b_dif);
          fprintf (f, "%e  ", strct1->Rbn98);       // RBN98
          fprintf (f, "%e  ", strct1->Mbn98);
          fprintf (f, "%e  ", strct1->msbn98_str);
          fprintf (f, "%e  ", strct1->msbn98_dif);
      	  fprintf (f, "%e  ", strct2->Pos[0]);
      	  fprintf (f, "%e  ", strct2->Pos[1]);
      	  fprintf (f, "%e  ", strct2->Pos[2]);
          fprintf (f, "%d  ", strct1->nlowres);
          fprintf (f, "%d  ", strct1->tlowres);
          fprintf (f, "%e  ", strct1->rlowres);
          fprintf (f, "\n");
        }
      }
      fclose (f);
    } // if tasks[itasks]
  }// Loop over tasks


  // Free memory
  free (files_to_read);
  free (npartinfile);
  free (strct_to_get);
  free (P);
}
