/*
 *
 *  \file    test.c
 *  \brief
 *
 *
 */


#include "base.h"
#include "typedef.h"
#include "archive.h"
#include "catalog.h"
#include "simulation.h"


typedef struct Options
{
  Archive        param;
  Archive        output;
  Catalog        catalog;
  Simulation     simulation;
} Options;


void  test_usage   (int opt,  char ** argv);
int   test_options (int argc, char ** argv, Options * opt);
void  test_params  (Options * opt);

int main (int argc, char ** argv)
{
  int           i, j, k, l, n, m;
  int           ndim     = 3;
  int           nregions = 8;
  int           coord;
  int           ndiv[ndim];
  int           regionstrctid[nregions];
  int         * strct_to_get;

  double        Lbox;
  double        pos[ndim];
  double        delta[ndim];
  double        regionmaxmass[nregions];
  double        limits [nregions][2][ndim];

  gheader       header;
  Options       opt;
  Structure   * strct;


  test_options (argc, argv, &opt);
  test_params  (&opt);

  Simulation_init         (&opt.simulation);
  Catalog_init            (&opt.catalog);
  Catalog_load_properties (&opt.catalog);

  strct_to_get = (int *) malloc ((opt.catalog.nstruct+1) * sizeof(int));
  Lbox    = opt.simulation.Lbox;

  for (i = 0; i < ndim; i++)
    ndiv[i] = 1;

  coord = 0;
  n = nregions;
  while (n /= 2)
  {
    ndiv[coord]++;
    coord = (coord + 1) % ndim;
  };

  for (i = 0; i < ndim; i++)
    delta[i] = Lbox / (double) ndiv[i];

  for (l = 0; l < nregions; l++)
  {
    m = l;
    for (n = 0; n < ndim; n++)
    {
      limits[l][0][n] = (m % ndiv[n]) * delta[n];
      limits[l][1][n] = limits[l][0][n] + delta[n];
      m /= ndiv[n];
    }
  }


  for (l = 0; l < nregions; l++)
  {
    regionmaxmass[l] = 0.0;
    regionstrctid[l] = 0;
  }

  // To determine region of Particle / Structure
  for (i = n; n <= opt.catalog.nstruct; n++)
  {
    strct = &opt.catalog.strctProps[n];
    if (strct->Type == 7)
    {
      i = (int) (strct->Pos[0] / delta[0]) % ndiv[0];
      j = (int) (strct->Pos[1] / delta[1]) % ndiv[1] ;
      k = (int) (strct->Pos[2] / delta[2]) % ndiv[2] ;
      m = i + (j * ndiv[0]) + (k * ndiv[0] * ndiv[1]);

      if (strct->TotMass > regionmaxmass[m])
      {
        regionmaxmass[m] = strct->TotMass;
        regionstrctid[m] = strct->ID;
      }
    }
  }

  for (l = 0; l <= nregions; l++)
  {
    k = regionstrctid[l];
    strct_to_get[k] = 1;
    for (i = 1; i <= opt.catalog.nstruct; i++)
      if (opt.catalog.strctProps[i].HostID == k)
        strct_to_get[i] = 1;
  }

  Structure_get_particle_properties (&opt.catalog, &opt.simulation, strct_to_get);

  for (i = 1, k = 0; i <= opt.catalog.nstruct; i++)
    if (strct_to_get[i] == 1)
      k += opt.catalog.strctProps[i].NumPart;

  Particle * P = (Particle *) malloc (k * sizeof(Particle));

  for (i = 1, k = 0; i <= opt.catalog.nstruct; i++)
  {
    if (strct_to_get[i] == 1)
      for (j = 0; j < opt.catalog.strctProps[i].NumPart; j++, k++)
        Particle_copy (&opt.catalog.strctProps[i].Part[j], &P[k]);
  }


  // Write header
  header.time           = opt.simulation.a;
  header.redshift       = (1.0 / header.time) - 1.0;
  header.flag_sfr       = 0;
  header.flag_feedback  = 0;
  header.flag_cooling   = 0;
  header.BoxSize        = opt.simulation.Lbox;
  header.Omega0         = opt.simulation.cosmology.OmegaM;
  header.OmegaLambda    = opt.simulation.cosmology.OmegaL;
  header.HubbleParam    = opt.simulation.cosmology.HubbleParam;
  double sqrta = sqrt(header.time);
  for (i = 0; i < k; i++)
  {
    P[i].Pos[0] *= header.HubbleParam;
    P[i].Pos[1] *= header.HubbleParam;
    P[i].Pos[2] *= header.HubbleParam;
    P[i].Vel[0] /= sqrta;
    P[i].Vel[1] /= sqrta;
    P[i].Vel[2] /= sqrta;
    P[i].Mass   *= header.HubbleParam;
  }
  gadget_write_snapshot (P, k, &header, &opt.output);

  free (P);
  free (strct_to_get);
  Catalog_free (&opt.catalog);


  return (0);
}

//
//  Parameters
//
void test_params (Options * opt)
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

  // Simulation
  fscanf (opt->param.file, "%s  %s  %s  %s  %d", namebuff, prefixbuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->simulation.archive, namebuff);
  Archive_prefix (&opt->simulation.archive, prefixbuff);
  Archive_format (&opt->simulation.archive, frmtbuff);
  Archive_path   (&opt->simulation.archive, pathbuff);
  Archive_nfiles (&opt->simulation.archive, nflsbuff);

  // Catalog
  fscanf (opt->param.file, "%s  %s  %s  %s  %d", namebuff, prefixbuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->catalog.archive, namebuff);
  Archive_prefix (&opt->catalog.archive, prefixbuff);
  Archive_format (&opt->catalog.archive, frmtbuff);
  Archive_path   (&opt->catalog.archive, pathbuff);
  Archive_nfiles (&opt->catalog.archive, nflsbuff);

  // Output
  fscanf (opt->param.file, "%s  %s  %s  %s  %d", namebuff, prefixbuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->output, namebuff);
  Archive_prefix (&opt->output, prefixbuff);
  Archive_format (&opt->output, frmtbuff);
  Archive_path   (&opt->output, pathbuff);
  Archive_nfiles (&opt->output, nflsbuff);

  // Close
  fclose (opt->param.file);
}

//
//  Options
//
int test_options (int argc, char ** argv, Options * opt)
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

  while ((myopt = getopt_long (argc, argv, "p:vh", lopts, &index)) != -1)
  {
    switch (myopt)
    {
      case 'p':
      	strcpy (opt->param.name, optarg);
        flag++;
        break;

      case 'h':
      	test_usage (0, argv);
        break;

      default:
      	test_usage (1, argv);
    }
  }

  if (flag == 0)
    test_usage (1, argv);
}


//
//  Usage
//
void test_usage (int opt, char ** argv)
{
  if (opt == 0)
  {
    printf ("                                                                         \n");
    printf ("  test                                                                   \n");
    printf ("                                                                         \n");
    printf ("  Author:            Rodrigo Can\\~as                                    \n");
    printf ("                                                                         \n");
    printf ("  Last edition:      12 - 07 - 2018                                      \n");
    printf ("                                                                         \n");
    printf ("                                                                         \n");
    printf ("  Usage:             %s [Option] [Parameter [argument]] ...\n",      argv[0]);
    printf ("                                                                         \n");
    printf ("  Parameters:                                                            \n");
    printf ("                                                                         \n");
    printf ("                     -p    --param    [string]   Name of param file      \n");
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
