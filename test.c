/*
 *
 *  \file    ctlgmatch.c
 *  \brief   This code cross matches structures (galaxies, dmhs, ...)
 *           in the same snapshot identified by two different codes.
 *
 *           Currently only VELOCIraptor - HaloMaker cross match is
 *           available.
 *
 */


#include "astrolibc.h"
#include "allvars.h"
#include "stf.h"
#include "halomaker.h"
#include "ramses.h"
#include "hdf5.h"
#include "hdf5sim.h"


typedef struct Options
{
  int            verbose;
  Archive        param;
  Archive        output;
  Catalog     *  catalog;
  Simulation  *  simulation;
} Options;


void test_usage   (int opt,  char ** argv)
int  test_options (int argc, char ** argv, Options * opt)


int main (int argc, char ** argv)
{
  int          i, j, k;
  Options      opt;
  Structure  * strct;
  Structure  * strct1;
  Structure  * strct2;


  test_options (argc, argv, &opt);
  test_params  (&opt);


  //
  //  Load catalogs
  //
  printf ("simulation init\n");
  Simulation_init (&opt.simulation);

  printf ("BoxSize          %g\n", opt.simulation.Lbox);
  printf ("HubbleParam      %g\n", opt.simulation.cosmology.HubbleParam);
  printf ("Om0              %g\n", opt.simulation.cosmology.OmegaM);
  printf ("OmB              %g\n", opt.simulation.cosmology.OmegaB);
  printf ("OmL              %g\n", opt.simulation.cosmology.OmegaL);
  printf ("time             %g\n", opt.simulation.Time);
  printf ("z                %g\n", opt.simulation.z);
  for (i = 0; i < 6; i++)
    printf ("NumPartThisFile  %u\n", opt.simulation.NpartThisFile[i]);
  for (i = 0; i < 6; i++)
    printf ("NumPartTot       %u\n", opt.simulation.NpartTot[i]);
  for (i = 0; i < 6; i++)
    printf ("Mass             %g\n", opt.simulation.MassTable[i]);

  /*
  printf ("Catalog init\n");
  Catalog_init (&opt.catalog);

  printf ("Catalog load\n");
  Catalog_load (&opt.catalog);

  printf ("Loading particle properties\n");
  Catalog_get_particle_properties (&opt.catalog, &opt.simulation);

  Catalog_free (&opt.catalog);
 */

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

  // Catalogs to load
  opt->catalog    = (Catalog *)    malloc (sizeof(Catalog));
  opt->simulation = (Simulation *) malloc (sizeof(Simulation));

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
    {"input",   0, NULL, 'i'},
    {0,         0, NULL, 0}
  };

  while ((myopt = getopt_long (argc, argv, "i:vh", lopts, &index)) != -1)
  {
    switch (myopt)
    {
      case 'i':
      	strcpy (opt->param.name, optarg);
        flag++;
        break;

      case 'v':
      	opt->verbose = 1;
      	break;

      case 'h':
      	ctlgMatch_usage (0, argv);
        break;

      default:
      	ctlgMatch_usage (1, argv);
    }
  }

  if (flag == 0)
    ctlgMatch_usage (1, argv);

}


//
//  Usage
//
void test_usage (int opt, char ** argv)
{
  if (opt == 0)
  {
    printf ("                                                                         \n");
    printf ("  test.c                                                                 \n");
    printf ("                                                                         \n");
    printf ("  Author:            Rodrigo Can\\~as                                    \n");
    printf ("                                                                         \n");
    printf ("  Last edition:      16 - 03 - 2018                                      \n");
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
