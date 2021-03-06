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
//#include <mpi.h>
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


void  isgal_usage   (int opt,  char ** argv);
int   isgal_options (int argc, char ** argv, Options * opt);
void  isgal_params  (Options * opt);


int main (int argc, char ** argv)
{
  int myrank = 0;
  int nprocs = 1;

  //MPI_Init (&argc, &argv);
  //MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
  //MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

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
  isgal_options (argc, argv, &opt);
  isgal_params  (&opt);


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
    strct1->flg_CorrectedPeriodicity = 0;
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
  double     msum2;
  double     msum3;
  double     rad;
  double     rho;
  int        ninR200;
  int        n3dfof;
  double     galpos [3];
  double     prev   [3];
  Particle * partinR200;
  FILE     * fout;
  int        chkpnt;

  stfExtendedOutput * xtndd;
  int                 ninextended;
  int                 indx;
  int                 id;


  //
  // Loop over tasks
  //
  //for (itasks = 0; itasks < opt.catalog.archive.nfiles; itasks++)
  //for (itasks = myrank*2; itasks < ((myrank+1)*2); itasks++)
  for (itasks = 0; itasks < opt.catalog.archive.nfiles; itasks++)
  {
    chkpnt = 0;
    start_t = clock();
    printf ("%d  opening catalog file %d\n", myrank, itasks);
    fflush (stdout);


    // Reset array values
    for (i = 0; i <= opt.catalog.nstruct; i++)
      strct_to_get[i] = 0;

    for (i = 0; i < opt.simulation.archive.nfiles; i++)
    {
      files_to_read[i] = 0;
      npartinfile[i]   = 0;
    }

    // Tag files to read of fof of interest
    //for (i = 1; i <= 5; i++)
    for (i = 1; i <= opt.catalog.nstruct; i++)
    {
      strct1 = &opt.catalog.strctProps[i];
      //if ((strct1->oTask == itasks) && (strct1->Type == 10))
      if ((strct1->oTask == itasks) && (strct1->Type >= 10))
          strct_to_get[i] = 1;
    }

    Structure_get_particle_properties (&opt.catalog, &opt.simulation, strct_to_get);    
    
    // No need to shift to CM or correct periodicity
    // this is done on the functions

    // Calculate Pos Sigma Tensor
    Structure_calculate_disp_tensor_pos (&opt.catalog, &opt.simulation, strct_to_get);
        
    // Calculate Vel Sigma Tensor
    Structure_calculate_disp_tensor_vel (&opt.catalog, &opt.simulation, strct_to_get);
 
    //
    sprintf (fname, "isgal_subs.txt.%d", itasks);
    fout = fopen (fname, "w");
    for (i = 1; i <= opt.catalog.nstruct; i++)
    {
      strct1 = &opt.catalog.strctProps[i];              
      if (strct_to_get[i])
      {
        strct2 = &opt.catalog.strctProps[strct1->HostID]; // Central
        fprintf (fout, "%5d ", strct1->ID);               // ID gal
        fprintf (fout, "%5d ", strct2->ID);               // ID ihsc
        fprintf (fout, "%e  ", strct1->TotMass);          // Mass IHSC
        fprintf (fout, "%e  ", strct1->Efrac);            // Efrac
        fprintf (fout, "%e  ", strct1->sigmaPosEval[0]);   // Pos XX 
        fprintf (fout, "%e  ", strct1->sigmaPosEval[1]);   // Pos YY 
        fprintf (fout, "%e  ", strct1->sigmaPosEval[2]);   // Pos ZZ 
        fprintf (fout, "%e  ", strct1->sigmaVelEval[0]);   // Vel XX 
        fprintf (fout, "%e  ", strct1->sigmaVelEval[1]);   // Vel YY 
        fprintf (fout, "%e  ", strct1->sigmaVelEval[2]);   // Vel ZZ 
        fprintf (fout, "\n");
      }
    }
    fclose (fout);
    end_t = clock();
    printf ("%d  took %f minutes\n", myrank, (end_t - start_t)/CLOCKS_PER_SEC/60.0);
    fflush (stdout);
  } // Loop over tasks


  // Free memory
  Catalog_free (&opt.catalog);
  free (strct_to_get);
  free (files_to_read);
  return 0;

}


// --------------------------------------------------- //
//  Parameters
// --------------------------------------------------- //
void isgal_params (Options * opt)
{
  int   i;
  int   dummy;
  char  buffer     [NAME_LENGTH];
  char  namebuff   [NAME_LENGTH];
  char  prefixbuff [NAME_LENGTH];
  char  frmtbuff   [NAME_LENGTH];
  char  pathbuff   [NAME_LENGTH];
  int   nflsbuff;

  opt->param.file = fopen (opt->param.name, "r");
  if (opt->param.file == NULL)
  {
    printf ("Couldn't open file  %s\n", opt->param.name);
    printf ("Exiting...\n");
    exit (0);
  }

  // Output
  fscanf (opt->param.file, "%s  %s  %s  %s  %d", namebuff, prefixbuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->output, namebuff);
  Archive_prefix (&opt->output, prefixbuff);
  Archive_format (&opt->output, frmtbuff);
  Archive_path   (&opt->output, pathbuff);
  Archive_nfiles (&opt->output, nflsbuff);

  // Catalogues
  fscanf (opt->param.file, "%s  %s  %s  %s  %d", namebuff, prefixbuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->catalog.archive, namebuff);
  Archive_prefix (&opt->catalog.archive, prefixbuff);
  Archive_format (&opt->catalog.archive, frmtbuff);
  Archive_path   (&opt->catalog.archive, pathbuff);
  Archive_nfiles (&opt->catalog.archive, nflsbuff);

  // Simulation
  fscanf (opt->param.file, "%s  %s  %s  %s  %d", namebuff, prefixbuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->simulation.archive, namebuff);
  Archive_prefix (&opt->simulation.archive, prefixbuff);
  Archive_format (&opt->simulation.archive, frmtbuff);
  Archive_path   (&opt->simulation.archive, pathbuff);
  Archive_nfiles (&opt->simulation.archive, nflsbuff);

  // Close
  fclose (opt->param.file);
}

// --------------------------------------------------- //
//  Options
// --------------------------------------------------- //
int isgal_options (int argc, char ** argv, Options * opt)
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

      case 'h':
      	isgal_usage (0, argv);
        break;

      default:
      	isgal_usage (1, argv);
    }
  }

  if (flag == 0)
    isgal_usage (1, argv);
}


// --------------------------------------------------- //
//  Usage
// --------------------------------------------------- //
void isgal_usage (int opt, char ** argv)
{
  if (opt == 0)
  {
    printf ("                                                                         \n");
    printf ("  isgal                                                                 \n");
    printf ("                                                                         \n");
    printf ("  Author:            Rodrigo Can\\~as                                    \n");
    printf ("                                                                         \n");
    printf ("  Last edition:      04 - 08 - 2019                                      \n");
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
