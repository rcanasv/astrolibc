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
  int            verbose;
  int            ilist;
  int            id;
  Archive        list;
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
  int          i, j, k, l;
  Options      opt;


  test_options (argc, argv, &opt);
  test_params  (&opt);


  //
  //  Load catalogs
  //
  Simulation_init                 (&opt.simulation);
  Catalog_init                    (&opt.catalog);
  Catalog_load_properties         (&opt.catalog);
  Catalog_get_particle_properties (&opt.catalog, &opt.simulation);
  Catalog_fill_isolated           (&opt.catalog);


  Structure  ** strct;
  Structure   * strct1;
  Structure   * strct2;
  Structure   * strct3;

  int         * strct_to_get;

  int        ** in3dfof;
  Structure     fof3d;

  Structure     ihsc;
  gheader       header;


  //
  //  Strct_to_get  tells which Centrals + IHSC use
  //
  strct_to_get = (int *) malloc ((opt.catalog.nstruct+1) * sizeof(int));
  for (i = 1; i <= opt.catalog.nstruct; i++)
    strct_to_get[i] = 0;
  strct_to_get[opt.id] = 1; // This should be the FOF ID


  //
  // Get IDs of members of 3DFOF
  //
  for (i = 1; i <= opt.catalog.nstruct; i++)
  {
    strct1 = &opt.catalog.strctProps[i];
    if (strct1->Type > 7 && strct1->HostID == opt.id)
    {
      if (strct1->Central)
        strct_to_get[i] = 1;
      else
        strct_to_get[i] = 2;
    }
  }

  //
  // At this point
  //      1 = IHSC
  //      1 = Central
  //      2 = the rest

  //
  // Merge central with IHSM and fill 3DFOF
  //
  int totpart = 0;
  int n;
  for (i = 1, n = 0; i <= opt.catalog.nstruct; i++)
  {
    strct1 = &opt.catalog.strctProps[i];
    if (strct_to_get[i])
      totpart += strct1->NumPart;
  }

  // Allocate memory for structure
  strct1 = &opt.catalog.strctProps[opt.id];
  strct2 = &opt.catalog.strctProps[0];
  strct2->NumPart = totpart;
  for (i = 0; i < 3; i++)
  {
    strct2->Pos[i] = strct1->Pos[i];
    strct2->Vel[i] = strct1->Vel[i];
  }
  strct2->Part = (Particle *) malloc (totpart * sizeof(Particle));
  strct_to_get[0] = 1;

  // Copy particles to
  for (i = 1, n = 0; i <= opt.catalog.nstruct; i++)
  {

    strct1 = &opt.catalog.strctProps[i];
    if (strct_to_get[i])
    {
      for (k = 0; k < strct1->NumPart; k++)
        Particle_copy (&strct1->Part[k], &strct2->Part[n++]);
    }
  }

  //
  for (i = 0; i <= opt.catalog.nstruct; i++)
  {
    strct1 = &opt.catalog.strctProps[i];
    if (strct_to_get[i] == 1)
      Structure_correct_periodicity (strct1, &opt.simulation);
  }

  //
  //  Calculate Surface density and create files
  //
  FILE      * f;
  double    * r     = NULL;
  double    * Rho   = NULL;
  char        buffer [NAME_LENGTH];
  int         bob;
  int         nbins = 200;
  double      rmin = 0.0;
  double      rmax = 400.0;


  for (i = 0, k = 0; i <= opt.catalog.nstruct; i++)
  {
    strct1 = &opt.catalog.strctProps[i];
    if (strct_to_get[i] == 1)
    {
      // Density profile
      Structure_calculate_spherical_density (strct1, rmin, rmax, nbins, 0, &r, &Rho);
      sprintf (buffer, "rho_density_%07d_%02d.dat", opt.id, k);
      f = fopen (buffer, "w");
      for (j = 0; j < nbins; j++)
        fprintf (f, "%e  %e\n", r[j], Rho[j]);
      fclose (f);

      // Visualization
      sprintf (opt.output.name, "strct_%07d.gdt_%02d", opt.id, k);
      gadget_write_snapshot (strct1->Part, strct1->NumPart, &header, &opt.output);
      k++;
    }
  }

  //
  // Free Memory
  //
  strct1 = &opt.catalog.strctProps[0];
  free (strct1->Part);
  free (strct_to_get);
  Catalog_free (&opt.catalog);

  printf ("All done, memory freed\n");

  return (0);
}


//
//  Parameters
//
void test_params (Options * opt)
{
  int   i;
  int   dummy;
  char  buffer [NAME_LENGTH];

  opt->param.file = fopen (opt->param.name, "r");
  if (opt->param.file == NULL)
  {
    printf ("Couldn't open file  %s\n", opt->param.name);
    printf ("Exiting...\n");
    exit (0);
  }

  // Catalog
  fscanf (opt->param.file, "%s", buffer);  Archive_name   (&opt->catalog.archive, buffer);
                                           Archive_prefix (&opt->catalog.archive, buffer);
  fscanf (opt->param.file, "%s", buffer);  Archive_format (&opt->catalog.archive, buffer);
  fscanf (opt->param.file, "%s", buffer);  Archive_path   (&opt->catalog.archive, buffer);
  fscanf (opt->param.file, "%d", &dummy);  Archive_nfiles (&opt->catalog.archive, dummy);

  // Simulation
  fscanf (opt->param.file, "%s", buffer);  Archive_name   (&opt->simulation.archive, buffer);
                                           Archive_prefix (&opt->simulation.archive, buffer);
  fscanf (opt->param.file, "%s", buffer);  Archive_format (&opt->simulation.archive, buffer);
  fscanf (opt->param.file, "%s", buffer);  Archive_path   (&opt->simulation.archive, buffer);
  fscanf (opt->param.file, "%d", &dummy);  Archive_nfiles (&opt->simulation.archive, dummy);

  // Output
  fscanf (opt->param.file, "%s", buffer);  Archive_name   (&opt->output, buffer);
                                           Archive_prefix (&opt->output, buffer);
  fscanf (opt->param.file, "%s", buffer);  Archive_format (&opt->output, buffer);
  fscanf (opt->param.file, "%s", buffer);  Archive_path   (&opt->output, buffer);
  fscanf (opt->param.file, "%d", &dummy);  Archive_nfiles (&opt->output, dummy);

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
    {"help",    0, NULL, 'h'},
    {"verbose", 0, NULL, 'v'},
    {"param",   1, NULL, 'p'},
    {"id",      2, NULL, 'i'},
    {0,         0, NULL, 0}
  };


  while ((myopt = getopt_long (argc, argv, "i:p:vh", lopts, &index)) != -1)
  {
    switch (myopt)
    {
      case 'p':
      	strcpy (opt->param.name, optarg);
        flag++;
        break;

      case 'i':
        opt->id = atoi(optarg);
        flag++;
        break;

      case 'v':
      	opt->verbose = 1;
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
    printf ("  density.c                                                              \n");
    printf ("                                                                         \n");
    printf ("  Author:            Rodrigo Can\\~as                                    \n");
    printf ("                                                                         \n");
    printf ("  Last edition:      17 - 09 - 2019                                      \n");
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
