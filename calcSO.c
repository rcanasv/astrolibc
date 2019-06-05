/*
 *
 *  \file    calcSO.c
 *  \brief
 *
 *
 */

#include "base.h"
#include "typedef.h"
#include "archive.h"
#include "catalog.h"
#include "simulation.h"
#include <mpi.h>
#include <time.h>


typedef struct Options
{
  int            iVerbose;
  int            iExtract;
  int            iSO;
  int            nsnap;
  Archive        param;
  Archive        output;
  Catalog        catalog;
  Simulation     simulation;
} Options;


void  calcSO_usage   (int opt,  char ** argv);
int   calcSO_options (int argc, char ** argv, Options * opt);
void  calcSO_params  (Options * opt);


int main (int argc, char ** argv)
{
  int myrank;
  int nprocs;

  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

  int           i, j, k, l, n, m;
  Options       opt;

  Particle    * P;
  Structure   * strct1;
  Structure   * strct2;
  Structure   * strct3;
  Structure   * sorted;
  Structure     tmpstrct;
  Grid        * myGrid;
  gheader       header;

  FILE  * f;
  char    fname  [NAME_LENGTH];
  char    buffer [NAME_LENGTH];
  int     tmpid;
  int     nfiles;
  int     numpart;

  int  * strct_to_get   = NULL;
  int  * files_to_read  = NULL;
  int  * files_of_strct = NULL;

  int  * strct_of_fof = NULL;
  int  * files_of_fof = NULL;

  int       * npartinfile    = NULL;
  Particle  * partbuffer     = NULL;
  Particle ** allPart        = NULL;

  double      dmp_mass;
  double      mratio;

  int    ngaspart;
  int    ndmpart;
  int    nstarpart;
  int    totpart;
  int    nghost;

  int    itasks;
  int    ifiles;

  int     dummyi;

  clock_t  start_t, end_t;


  // Read params
  calcSO_options (argc, argv, &opt);
  calcSO_params  (&opt);


  // Load sim and catalog information
  Simulation_init                 (&opt.simulation);
  Catalog_init                    (&opt.catalog);
  Catalog_load_properties         (&opt.catalog);
  Catalog_fill_SubIDS             (&opt.catalog);
  Catalog_fill_isolated           (&opt.catalog);


  // Tag central galaxy
  for (j = 1; j <= opt.catalog.nstruct; j++)
  {
    strct1 = &opt.catalog.strctProps[j];
    strct1->dummyd = 0.0;
    strct1->dummyi = 0;
    strct1->dummy  = 0;
  }

  for (j = 1; j <= opt.catalog.nstruct; j++)
  {
    strct1 = &opt.catalog.strctProps[j];
    if (strct1->Central == 1 && strct1->HostID > 0)
    {
      strct2 = &opt.catalog.strctProps[strct1->HostID];
      strct1->dummyd = strct1->TotMass + strct2->TotMass;
      strct2->dummyi = strct1->ID;
    }
  }

  for (j = 1; j <= opt.catalog.nstruct; j++)
  {
    strct1 = &opt.catalog.strctProps[j];
    if (strct1->Central != 1 && strct1->Type > 7)
    {
      strct2 = &opt.catalog.strctProps[strct1->HostID];
      strct3 = &opt.catalog.strctProps[strct2->dummyi];
      strct3->dummyd += strct1->TotMass;
    }
  }


  // Allocate arrays
  strct_to_get  = (int *)       malloc ((opt.catalog.nstruct+1)        *sizeof(int));
  files_to_read = (int *)       malloc ((opt.simulation.archive.nfiles)*sizeof(int));
  files_of_fof  = (int *)       malloc ((opt.simulation.archive.nfiles)*sizeof(int));
  npartinfile   = (int *)       malloc ((opt.simulation.archive.nfiles)*sizeof(int));
  myGrid        = (Grid *)      malloc ((opt.simulation.archive.nfiles)*sizeof(Grid));
  allPart       = (Particle **) malloc ((opt.simulation.archive.nfiles)*sizeof(Particle *));

  dmp_mass = 1.0 / (1024.0*1024.0*1024.0) \
             * (opt.simulation.cosmology.OmegaM - opt.simulation.cosmology.OmegaB) \
             /  opt.simulation.cosmology.OmegaM * opt.simulation.unit_m;

  double lbox   = opt.simulation.Lbox;
  double lbox_2 = lbox / 2.0;
  double unit_m = opt.simulation.unit_m;
  double dx;
  double vol;
  int    ninbuffer;

  // Variables for R200 stuff
  double     pi      = acos(-1.0);
  double     fac     = 4.0 * pi / 3.0;
  double     G       = 43009.1e-10;                                    // in (kpc/M_sun)*(km/s)^2
  double     H       = opt.simulation.cosmology.HubbleParam / 1000.0;  // in (km/s)/kpc
  double     crit    = 3.0 * H * H / (8.0 * pi * G);
  double     crit200 = 200.0 * crit;
  double     msum;
  double     rad;
  double     rho;
  int        ninR200;
  double     galpos [3];
  double     prev   [3];
  Particle * partinR200;
  FILE * fout;

  //
  // Loop over tasks
  //
  //for (itasks = 0; itasks < opt.catalog.archive.nfiles; itasks++)
  //for (itasks = myrank*2; itasks < ((myrank+1)*2); itasks++)
  for (itasks = myrank; itasks <= myrank; itasks++)
  {
    start_t = clock();
    printf ("%d  opening catalog file %d\n", myrank, itasks);
    fflush (stdout);
    sprintf (fname, "R200.txt.%d", itasks);
    fout = fopen (fname, "w");

    // Reset array values
    for (i = 0; i <= opt.catalog.nstruct; i++)
      strct_to_get[i] = 0;

    for (i = 0; i < opt.simulation.archive.nfiles; i++)
      files_to_read[i] = 0;


    // Tag files to read of fof of interest
    sprintf (fname, "%s/%s.filesofgroup", opt.catalog.archive.path, opt.catalog.archive.name);
    f = fopen (fname, "r");
    for (i = 1; i <= opt.catalog.nstruct; i++)
    {
      strct1 = &opt.catalog.strctProps[i];

      fgets  (buffer, NAME_LENGTH, f);
      sscanf (buffer, "%d  %d", &tmpid, &nfiles);
      fgets  (buffer, NAME_LENGTH, f);

      get_n_num_from_string (buffer, nfiles, &files_of_strct);

      strct1->NumFiles = nfiles;
      strct1->FilesOfGroup = (int *) malloc (nfiles*sizeof(int));
      for (j = 0; j < nfiles; j++)
         strct1->FilesOfGroup[j] = files_of_strct[j];

       if ((strct1->oTask == itasks) && (strct1->Type == 7) && (strct1->NumSubs > 0))
       {
         strct_to_get[i] = 1;
         for (j = 0; j < nfiles; j++)
           files_to_read[files_of_strct[j]] = 1;
       }
      free (files_of_strct);
    }
    fclose (f);


    // Tag files to read of all
    for (i = 1; i <= opt.catalog.nstruct; i++)
    {
      strct1 = &opt.catalog.strctProps[i];
      if ((strct1->Type > 7) && (strct_to_get[strct1->HostID]))
        for (j = 0; j < strct1->NumFiles; j++)
          files_to_read[strct1->FilesOfGroup[j]] = 1;
    }


    // Load over files to read particles
    for (n = 0; n < opt.simulation.archive.nfiles; n++)
    {
      if(files_to_read[n])
      {
        // Reset particle counters
        ngaspart  = 0;
        ndmpart   = 0;
        nstarpart = 0;
        totpart   = 0;
        nghost    = 0;

        // Load AMR + gas particles
        ramses_amr_init   (&myGrid[n]);
        ramses_amr_load   (&opt.simulation, n, &myGrid[n]);
        ramses_hydro_read (&opt.simulation, n, &myGrid[n]);

        for (k = myGrid[n].nlevelmax-1; k >= 9; k--)
          for (i = 0; i < myGrid[n].level[k].num; i++)
            for (j = 0; j < 8; j++)
              if (myGrid[n].level[k].cell[i].okOct[j])
                ngaspart++;

        // Load dark matter and star particles
        ramses_load_particles (&opt.simulation, n, &partbuffer);
        for (i = 0; i < opt.simulation.npartinfile[n]; i++)
        {
          mratio = fabs(partbuffer[i].Mass/dmp_mass - 1);
          if (mratio < 1e-4 && partbuffer[i].Age == 0)
          {
            partbuffer[i].Type = 1;
            ndmpart++;
          }
          else
          {
            if (partbuffer[i].Age != 0)
            {
              partbuffer[i].Type = 4;
              nstarpart++;
            }
            else
            {
              partbuffer[i].Type = -1;
              nghost++;
            }
          }
        } // Loop over particles in file


        // Add total particles in file(s) with index n
        totpart = ngaspart + ndmpart + nstarpart;
        npartinfile[n] = totpart;
        /*
        printf ("file_index   %8d  ", n);
        printf ("ngaspart      %8d  ", ngaspart);
        printf ("ndmpart      %8d  ", ndmpart);
        printf ("nstarpart    %8d  ", nstarpart);
        printf ("nghost       %8d  ", nghost);
        printf ("nparttot     %8d\n", totpart);
        */
        allPart[n] = (Particle *) malloc (totpart * sizeof(Particle));

        // Copy gas particles to array
        m = 0;
        for (k = myGrid[n].nlevelmax-1; k >= 9; k--)
        {
          for (i = 0; i < myGrid[n].level[k].num; i++)
          {
            for (j = 0; j < 8; j++)
            {
              if (myGrid[n].level[k].cell[i].okOct[j])
              {
                allPart[n][m].Pos[0] = lbox * myGrid[n].level[k].cell[i].octPos[j][0];
                allPart[n][m].Pos[1] = lbox * myGrid[n].level[k].cell[i].octPos[j][1];
                allPart[n][m].Pos[2] = lbox * myGrid[n].level[k].cell[i].octPos[j][2];
                dx  = myGrid[n].level[k].cell[i].dx;
                vol = dx * dx * dx;
                allPart[n][m].Mass = unit_m * vol * myGrid[n].level[k].cell[i].octRho[j];
                allPart[n][m].Id = m;
                allPart[n][m].Type = 2;
                m++;
              }
            }
          }
        }
        // Copy dm + star particles to array
        for (i = 0; i < opt.simulation.npartinfile[n]; i++)
          if ((partbuffer[i].Type == 1) || (partbuffer[i].Type == 4))
            Particle_copy (&partbuffer[i], &allPart[n][m++]);

        // Free memory for next file
        ramses_amr_free (&myGrid[n]);
        free (partbuffer);
      } // If file to read
    } // Loop over files to load particles

    // Calculate R200 for all centrals
    prev[0] = 0.0;
    prev[1] = 0.0;
    prev[2] = 0.0;

printf ("%d  Checkpoint\n", myrank);
fflush (stdout);

    //printf ("IHSC_ID  Ctrl_ID   R200          M200             Rho          Mc/Mvir\n");
    // Loop over FOFs
    for (k = 1; k <= opt.catalog.nstruct; k++)
    {
      strct1 = &opt.catalog.strctProps[k];
      if (strct_to_get[k])
      {
        for (n = 0; n < opt.simulation.archive.nfiles; n++)
          files_of_fof[n] = 0;

        // Tag index of files to look neighbouring particles
        for (i = 1; i <= opt.catalog.nstruct; i++)
        {
          strct2 = &opt.catalog.strctProps[i];
          if (strct2->HostID == strct1->ID || strct2->ID == strct1->ID)
            for (j = 0; j < strct2->NumFiles; j++)
              files_of_fof[strct2->FilesOfGroup[j]] = 1;
        }

        // Centre of R200 is central galaxy
        strct2 = &opt.catalog.strctProps[strct1->dummyi];
        ninbuffer = 0;
        for (n = 0; n < opt.simulation.archive.nfiles; n++)
        {
          if (files_of_fof[n])
          {
            for (j = 0; j < npartinfile[n]; j++)
            {

              allPart[n][j].Pos[0] -= strct2->Pos[0];
              allPart[n][j].Pos[1] -= strct2->Pos[1];
              allPart[n][j].Pos[2] -= strct2->Pos[2];

              while (allPart[n][j].Pos[0] > lbox_2) allPart[n][j].Pos[0] -= lbox;
              while (allPart[n][j].Pos[1] > lbox_2) allPart[n][j].Pos[1] -= lbox;
              while (allPart[n][j].Pos[2] > lbox_2) allPart[n][j].Pos[2] -= lbox;

              while (allPart[n][j].Pos[0] < -lbox_2) allPart[n][j].Pos[0] += lbox;
              while (allPart[n][j].Pos[1] < -lbox_2) allPart[n][j].Pos[1] += lbox;
              while (allPart[n][j].Pos[2] < -lbox_2) allPart[n][j].Pos[2] += lbox;

              if (fabs(allPart[n][j].Pos[0]) < 3000.0)
                if (fabs(allPart[n][j].Pos[1]) < 3000.0)
                  if (fabs(allPart[n][j].Pos[2]) < 3000.0)
                  {
                    allPart[n][j].dummyi = 1;
                    ninbuffer++;
                  }
            }
          }
        }
        partbuffer = (Particle *) malloc (ninbuffer * sizeof(Particle));

        for (n = 0, i = 0; n < opt.simulation.archive.nfiles; n++)
        {
          if (files_of_fof[n])
            for (j = 0; j < npartinfile[n]; j++)
            {
              if (allPart[n][j].dummyi)
              {
                Particle_copy (&allPart[n][j], &partbuffer[i]);
                Particle_get_radius (&partbuffer[i]);
                i++;
              }
              allPart[n][j].Pos[0] += strct2->Pos[0];
              allPart[n][j].Pos[1] += strct2->Pos[1];
              allPart[n][j].Pos[2] += strct2->Pos[2];
              allPart[n][j].dummyi = 0;
            }
        }

        prev[0] = strct2->Pos[0];
        prev[1] = strct2->Pos[1];
        prev[2] = strct2->Pos[2];
        qsort (partbuffer, ninbuffer, sizeof(Particle), Particle_rad_compare);

        //
        for (i = 0, msum = 0, ninR200 = 0; i < ninbuffer; i++)
        {
          msum += partbuffer[i].Mass;
          rad = partbuffer[i].Radius;
          rho = msum / (fac * rad * rad * rad);
          if (rho < crit200)
            break;
          else
            ninR200++;
        }
        fprintf (fout, "%6d   %6d   %e   %e   %e   %e   %d  %d\n", strct1->ID, strct2->ID, rad, msum, rho, log10(strct2->TotMass/msum), ninR200, ninbuffer);
fflush (fout);
        /*
        partinR200 = (Particle *) malloc (ninR200 * sizeof(Particle));
        for (i = 0; i < ninR200; i++)
        {
          Particle_copy (&partbuffer[i], &partinR200[i]);
          //partinR200[i].Type = 1;
        }
        sprintf (opt.output.name, "halo_%03d", strct1->ID);
        gadget_write_snapshot (partinR200, ninR200, &header, &opt.output);
        free (partinR200);
        */
        free (partbuffer);
      }
    }
    fclose (fout);
    end_t = clock();
    printf ("%d  took %f minutes\n", myrank, (end_t - start_t)/CLOCKS_PER_SEC/60.0);
    fflush (stdout);
    for (i = 0; i < opt.simulation.archive.nfiles; i++)
      if (files_to_read[i])
        free(allPart[i]);
  } // Loop over tasks


  // Free memory
  Catalog_free (&opt.catalog);
  free (myGrid);
  free (strct_to_get);
  free (npartinfile);
  free (files_to_read);
  free (allPart);
  MPI_Finalize();
  return 0;



  // --------------------------------------------------- //
  // Diffuse stellar fraction
  //FILE * f;
  //char buffer [NAME_LENGTH];

  if (opt.iSO)
  {
    int     nsat_m09;
    int     nsat_m10;
    int     nsat_m11;
    double  fmass;
    double  radius;
    double  minsat_m08;
    double  minsat_m09;
    double  minsat_m10;
    double  minsat_m11;

      sprintf (buffer, "%s.ihsc", opt.catalog.archive.prefix);
      f = fopen (buffer, "w");
      for (j = 1; j <= opt.catalog.nstruct; j++)
      {
        strct1 = &opt.catalog.strctProps[j];              // IHSC
        strct2 = &opt.catalog.strctProps[strct1->dummyi]; // Central
        if (strct1->Type == 7 && strct1->NumSubs > 0)
        {
          nsat_m09 = 0;
          nsat_m10 = 0;
          nsat_m11 = 0;
          minsat_m08 = 0.0;
          minsat_m09 = 0.0;
          minsat_m10 = 0.0;
          minsat_m11 = 0.0;

          for (k = 1; k < strct1->NumSubs; k++)
          {
            strct3 = &opt.catalog.strctProps[strct1->SubIDs[k]];

            if (strct3->TotMass >= 1e8)
              minsat_m08 += strct3->TotMass;
            if (strct3->TotMass >= 1e9)
              minsat_m09 += strct3->TotMass;
            if (strct3->TotMass >= 1e10)
              minsat_m10 += strct3->TotMass;
            if (strct3->TotMass >= 1e11)
              minsat_m11 += strct3->TotMass;

            if (strct3->TotMass >= 1e9  && strct3->TotMass < 1e10)
              nsat_m09++;
            if (strct3->TotMass >= 1e10 && strct3->TotMass < 1e11)
              nsat_m10++;
            if (strct3->TotMass >= 1e11)
              nsat_m11++;
          }

          /*
          radius = 2.0 * strct2->Rx;
          Structure_calculate_j_r       (strct2, radius);
          Structure_calculate_sigma_v_r (strct2, radius);
          Structure_calculate_sfr       (strct2);
          */

          strct3 = &opt.catalog.strctProps[strct1->SubIDs[strct1->NumSubs-2]];

          fprintf (f, "%e  ", strct2->dummyd);      // Total Stellar Mass
          fprintf (f, "%e  ", strct1->TotMass);     // Mass IHSC
          fprintf (f, "%e  ", strct2->TotMass);     // Mass Central
          fprintf (f, "%e  ", strct3->TotMass);     // Mass Second most massive Gal
          fprintf (f, "%5d ", strct1->NumSubs);     // NumSubs
          fprintf (f, "%5d ", strct2->Central);     // Is Central central?
          fprintf (f, "%5d ", strct1->ID);          // ID IHSC
          fprintf (f, "%5d ", strct2->ID);          // ID Central
          fprintf (f, "%5d ", strct3->ID);          // ID Second most
          fprintf (f, "%5d ", nsat_m09);            // Num sats 1e9  <= M < 1e10
          fprintf (f, "%5d ", nsat_m10);            // Num sats 1e10 <= M < 1e11
          fprintf (f, "%5d ", nsat_m11);            // Num sats 1e11 <= M
          /*
          fprintf (f, "%e ",  strct2->sigma);       // sigma_v(r)
          fprintf (f, "%e ",  strct2->j[3]);        // j(r)
          fprintf (f, "%e ",  radius);              // r
          fprintf (f, "%e ",  strct2->SFR20);       // SFR20
          fprintf (f, "%e ",  strct2->SFR50);       // SFR50
          fprintf (f, "%e ",  strct2->SFR100);      // SFR100
          */
          fprintf (f, "%e ", minsat_m08);            // Mass in sats M > 1e8
          fprintf (f, "%e ", minsat_m09);            // Mass in sats M > 1e9
          fprintf (f, "%e ", minsat_m10);            // Mass in sats M > 1e10
          fprintf (f, "%e ", minsat_m11);            // Mass in sats M > 1e11

          fprintf (f, "%10.5lf ", strct2->Pos[0]);   // Pos
          fprintf (f, "%10.5lf ", strct2->Pos[1]);   // Pos
          fprintf (f, "%10.5lf ", strct2->Pos[2]);   // Pos
          fprintf (f, "\n");
        }
      }
      fclose (f);
    }
  // --------------------------------------------------- //

  /*
  // --------------------------------------------------- //
  // Write snapshots for visualization
  if (opt.iExtract)
  {
    strct_to_get = (int *) malloc ((opt.catalog.nstruct+1)*sizeof(int));
    for (j = 1; j <= opt.catalog.nstruct; j++)
      strct_to_get[j] = 0;

    // Load Particles
    Structure_get_particle_properties (&opt.catalog, &opt.simulation, strct_to_get);

    // Write Gadget Snapshots
    for (i = opt.catalog.nstruct, k = 0; ((k < top)&&(i >=1)); i--)
    {
      printf ("boxsize %e\n", opt.simulation[0].Lbox);
      Structure_correct_periodicity (&tmpstrct, &opt.simulation);
      sprintf (opt.output.name, "%s_%d.gdt_%03d",opt.output.prefix, k, 0);
      gadget_write_snapshot (tmpstrct.Part, m, &header, &opt.output);
    }

    free (strct_to_get);
  }
  // --------------------------------------------------- //
  */

  // Free memory
  Catalog_free (&opt.catalog);

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
