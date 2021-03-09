/*
 *
 *  \file    SO.c
 *  \brief   Routine to calculate spherical overdensity
 *
 *
 */

#include "SO.h"

int get_structure_SO (Catalog * ctlg, Simulation * sim, int * tasks)
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
  // Overdensity    'int'    specifying overdensity critreia. See allvars.h for macro definition.
  //                Supported criteria:
  //                      - 'R200C'   200x rho_critical      strct.R200c
  //                      - 'R500C'   500x rho_critical      strct.R500c
  //                      - 'R200B'   200x rho_background    strct.R200b
  //                      - 'RBN98'   Bryan & Norman 1998    strct.Rbn98
  //                      - 'ALL'     all of the above
  //
  // NOTE:     - Overdensities are calculated at sim redshift
  //           - Particles will be stored in strct.SOPart array. In the case
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

  int     npart;
  int     ninrad;

  int  * files_to_read  = NULL;
  int  * loaded_files   = NULL;
  int  * strct_to_get   = NULL;
  int  * npartinfile    = NULL;


  // Useful variables
  double lbox   = sim->Lbox;
  double lbox_2 = lbox / 2.0;

  // Variables for R200 stuff
  double  pi  = acos(-1.0);
  double  fac = 4.0 * pi / 3.0;
  double  G   = 43009.1e-10;                          // in (kpc/M_sun)*(km/s)^2
  double  H   = sim->cosmology.HubbleParam / 1000.0;  // in (km/s)/kpc

  double  rhocrit  = 3.0 * H * H / (8.0 * pi * G);
  double  rho200c  = 200.0 * crit;
  double  rho500c  = 500.0 * crit;
  double  rho200b  = rhocrit;
  double  rhobn98  = rhocrit;


  // Allocate arrays
  files_to_read = (int *) malloc ((sim->archive.nfiles) * sizeof(int));
  npartinfile   = (int *) malloc ((sim->archive.nfiles) * sizeof(int));
  strct_to_get  = (int *) malloc ((ctg->nstruct)        * sizeof(int));

  // Initialize arrays
  for (i = 0; i < sim->archive.nfiles; i++)
  {
    files_to_read[i] = 0;
    npartinfile[i]   = 0;
  }
  for (i = 0; i < sim->archive.nfiles; i++)
    strct_to_get[i] = 0;


  // Loop over tasks
  for (itask = 0; itask < ctg->archive.nfiles; itask++)
  {
    if (tasks[itask])
    {
      // Tag files that need to be opened
      Catalog_get_files_of_groups (ctlg);
      // TODO: Simuation_tag_neighbour_files (sim, files_to_read);

      // Tag files to read object of interest
      for (i = 1; i <= opt.catalog.nstruct; i++)
      {
        strct1 = &opt.catalog.strctProps[i];
        if ((strct1->oTask == itask)  && \
            (strct1->Type == 7)       && \
            (strct1->NumSubs > 0))
        {
          // dummyi contains Central ID
          // dummyd contains total FOF mass
          strct2 = &opt.catalog.strctProps[strct1->dummyi];
          if (strct2->dummyd >= 1e10)
          {
            strct_to_get[i] = 1;
            for (j = 0; j < strct1->NumFiles; j++)
              files_to_read[strct1->FilesOfGroup[j]] = 1;
          }
        }
      }

      // Get number of particles to allocate array
      for (n = 0; n < sim->archive.nfiles; n++)
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
          Simulation_load_particles (sim, i, &Pbuff);
          for (i = 0; i < npartinfile[n]; i++)
          {
            Particle_copy (&Pbuff[i], &P[m++]);
            P.dummyi = 0;
          }
          free (Pbuff);
        }

      if (m != npart)
        printf ("Something's wrong mismatch on number of particles\n");

      // Loop over centrals
      for (k = 1; k <= opt.catalog.nstruct; k++)
      {
        strct1 = &opt.catalog.strctProps[k];
        if (strct_to_get[k] && strct1->inR200 == 0 && strct1->oTask == itask)
        {
          // Centre of R200 is central galaxy
          strct2 = &opt.catalog.strctProps[strct1->dummyi];
          ninrad = 0;

          for (j = 0; j < npart; j++)
          {
            P[j].Pos[0] -= strct2->Pos[0];
            P[j].Pos[1] -= strct2->Pos[1];
            P[j].Pos[2] -= strct2->Pos[2];

            Particle_correct_periodicity (&P[j], lbox_2);

            if ((fabs(P[j].Pos[0]) < 3000.0) && \
                (fabs(P[j].Pos[1]) < 3000.0) && \
                (fabs(P[j].Pos[2]) < 3000.0))
            {
              Particle_get_radius (&P[j]);
              ninrad++;
            }
            else
              P[j].Radius = lbox + j;
          }

          qsort (P, npart, sizeof(Particle), Particle_rad_compare);

          // Starts loop to get SO
          int icheck = 0;
          for (i = 0; i < ninrad; i++)
          {
            msum += P[i].Mass;
            rad = P[i].Radius;
            rho = msum / (fac * rad * rad * rad);

            if (rho > rho500c)  strct1->n500c++;  else {strct1->R500c = rad; strct1->M500c = msun; icheck++;}
            if (rho > rho200c)  strct1->n200c++;  else {strct1->R200c = rad; strct1->M200c = msun; icheck++;}
            if (rho > rho200b)  strct1->n200b++;  else {strct1->R200b = rad; strct1->M200b = msun; icheck++;}
            if (rho > rhobn98)  strct1->nbn98++;  else {strct1->Rbn98 = rad; strct1->Mbn98 = msun; icheck++;}

            if (icheck == 4)
            {
              strct1->nSO = i;
              break;
            }
          }

          strct1->PSO = (Particle *) malloc (strct1->nSO * sizeof(Particle));
          for (i = 0; i < strct1->nSO; i++)
          {
            Particle_copy (&P[i], strct1->PSO[i]);

            P[i].Pos[0] += strct2->Pos[0];
            P[i].Pos[1] += strct2->Pos[1];
            P[i].Pos[2] += strct2->Pos[2];
          }
        } // if strct_to_get
      } // loop over structures
    } // if tasks[itasks]
  }// Loop over tasks


  // Free memory
  free (files_to_read);
  free (npartinfile);
  free (strct_to_get);
  free (P);

  return (0);
}


// --------------------------------------------------- //
//  Parameters
// --------------------------------------------------- //
void calcSO_params (Options * opt)
{
  int   i;
  int   dummy;
  char  buffer   [NAME_LENGTH];
  char  namebuff [NAME_LENGTH];
  char  frmtbuff [NAME_LENGTH];
  char  pathbuff [NAME_LENGTH];
  int   nflsbuff;

  opt->param.file = fopen (opt->param.name, "r");
  if (opt->param.file == NULL)
  {
    printf ("Couldn't open file  %s\n", opt->param.name);
    printf ("Exiting...\n");
    exit (0);
  }

  // Output
  fscanf (opt->param.file, "%s  %s  %s  %d", namebuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->output, namebuff);
  Archive_prefix (&opt->output, namebuff);
  Archive_format (&opt->output, frmtbuff);
  Archive_path   (&opt->output, pathbuff);
  Archive_nfiles (&opt->output, nflsbuff);

  // Catalogues
  fscanf (opt->param.file, "%s  %s  %s  %d", namebuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->catalog.archive, namebuff);
  Archive_prefix (&opt->catalog.archive, namebuff);
  Archive_format (&opt->catalog.archive, frmtbuff);
  Archive_path   (&opt->catalog.archive, pathbuff);
  Archive_nfiles (&opt->catalog.archive, nflsbuff);

  // Simulation
  fscanf (opt->param.file, "%s  %s  %s  %d", namebuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->simulation.archive, namebuff);
  Archive_prefix (&opt->simulation.archive, namebuff);
  Archive_format (&opt->simulation.archive, frmtbuff);
  Archive_path   (&opt->simulation.archive, pathbuff);
  Archive_nfiles (&opt->simulation.archive, nflsbuff);

  // Close
  fclose (opt->param.file);
}

// --------------------------------------------------- //
//  Options
// --------------------------------------------------- //
int calcSO_options (int argc, char ** argv, Options * opt)
{
  int   myopt;
  int   index;
  int   flag = 0;

  extern char * optarg;
  extern int    opterr;
  extern int    optopt;

  struct option lopts[] = {
    {"help",      0, NULL, 'h'},
    {"verbose",   0, NULL, 'v'},
    {"param",     0, NULL, 'p'},
    {"so",        0, NULL, 's'},
    {"extract",   0, NULL, 'x'},
    {0,           0, NULL, 0}
  };

  while ((myopt = getopt_long (argc, argv, "p:ftxvh", lopts, &index)) != -1)
  {
    switch (myopt)
    {
      case 'p':
      	strcpy (opt->param.name, optarg);
        flag++;
        break;

      case 'f':
      	opt->iSO = 1;
      	break;

      case 'x':
      	opt->iExtract = 1;
      	break;

      case 'h':
      	calcSO_usage (0, argv);
        break;

      default:
      	calcSO_usage (1, argv);
    }
  }

  if (flag == 0)
    calcSO_usage (1, argv);
}


// --------------------------------------------------- //
//  Usage
// --------------------------------------------------- //
void calcSO_usage (int opt, char ** argv)
{
  if (opt == 0)
  {
    printf ("                                                                         \n");
    printf ("  calcSO                                                                 \n");
    printf ("                                                                         \n");
    printf ("  Author:            Rodrigo Can\\~as                                    \n");
    printf ("                                                                         \n");
    printf ("  Last edition:      22 - 05 - 2019                                      \n");
    printf ("                                                                         \n");
    printf ("                                                                         \n");
    printf ("  Usage:             %s [Option] [Parameter [argument]] ...\n",      argv[0]);
    printf ("                                                                         \n");
    printf ("  Parameters:                                                            \n");
    printf ("                                                                         \n");
    printf ("                     -i    --input    [string]   Name of input file      \n");
    printf ("                                                                         \n");
    printf ("  Options:                                                               \n");
    printf ("                     -v    --verbose             activate verbose        \n");
    printf ("                     -h    --help                displays description    \n");
    printf ("                                                                         \n");
    exit (0);
  }
  else
  {
    printf ("\t Error:           Some parameters are missing ...\n");
    printf ("\t Usage:           %s [Option] [Option [argument]] ...\n", argv[0]);
    printf ("\t For help try:    %s --help             \n", argv[0]);
    exit (0);
  }
}

/*
// First check that grid is properly loaded
for (k = 0; k < myGrid.nlevelmax; k++)
{
  sprintf (fname, "amr_lvl_%d", k);
  f = fopen(fname, "w");
  for (i = 0; i < myGrid.level[k].num; i++)
  {
    fprintf (f, "%e  ", myGrid.level[k].cell[i].Pos[0]);
    fprintf (f, "%e  ", myGrid.level[k].cell[i].Pos[1]);
    fprintf (f, "%e  ", myGrid.level[k].cell[i].Pos[2]);
    fprintf (f, "\n");
  }
  fclose (f);
}

sprintf (fname, "hydro_all");
f = fopen(fname, "w");
for (k = 9; k < myGrid.nlevelmax; k++)
{
  for (i = 0; i < myGrid.level[k].num; i++)
  {
    for (j = 0; j < 8; j++)
    {
      fprintf (f, "%e  ", myGrid.level[k].cell[i].octPos[j][0]);
      fprintf (f, "%e  ", myGrid.level[k].cell[i].octPos[j][1]);
      fprintf (f, "%e  ", myGrid.level[k].cell[i].octPos[j][2]);
      fprintf (f, "%e  ", myGrid.level[k].cell[i].octRho[j]);
      fprintf (f, "%d  ", myGrid.level[k].cell[i].okOct[j]);
      fprintf (f, "%d  ", k+1);
      fprintf (f, "\n");
    }
  }
}
fclose (f);

*/
