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

#include "base.h"
#include "typedef.h"
#include "archive.h"
#include "catalog.h"
#include "simulation.h"


typedef struct Options
{
  int            verbose;
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
  int          i, j, k;
  Options      opt;


  test_options (argc, argv, &opt);
  test_params  (&opt);


  //
  //  Load catalogs
  //
  //printf ("simulation init\n");
  Simulation_init (&opt.simulation);

  //printf ("Catalog init\n");
  Catalog_init (&opt.catalog);

  //printf ("Catalog load\n");
  Catalog_load_properties (&opt.catalog);

  //printf ("Loading particle properties\n");
  Catalog_get_particle_properties (&opt.catalog, &opt.simulation);

/*
  FILE * f;
  Structure * strct;
  char buffer [NAME_LENGTH];
  for (i = 1; i < opt.catalog.nstruct; i++)
  {
    sprintf (buffer, "galaxy_%03d", i);
    f = fopen (buffer, "w");
    strct = &opt.catalog.strctProps[i];
    for (j = 0; j < strct->NumPart; j++)
    {
      fprintf (f, "%e  ", strct->Part[j].Pos[0]);
      fprintf (f, "%e  ", strct->Part[j].Pos[1]);
      fprintf (f, "%e\n", strct->Part[j].Pos[2]);
    }
    fclose (f);
  }
*/

  Structure * strct;
  Particle  * p;

  double  cmx;
  double  cmy;
  double  cmz;

  double  totmass;
  double  mass100kpc3d;
  double  mass100kpc2d;

  double  halftotmass;
  double  halfmass100kpc3d;
  double  halfmass100kpc2d;

  double  r;
  double  r100kpc3d;
  double  r100kpc2d;

  double r2d;

  FILE * f;
  char buffer [NAME_LENGTH];
  sprintf (buffer, "sizemass_eagle_velociraptor.dat");
  f = fopen (buffer, "w");


  for (i = 1; i < opt.catalog.nstruct; i++)
  {
    strct = &opt.catalog.strctProps[i];

    totmass      = 0.0;
    mass100kpc3d = 0.0;
    mass100kpc2d = 0.0;

    halftotmass      = 0.0;
    halfmass100kpc3d = 0.0;
    halfmass100kpc2d = 0.0;

    r            = 0.0;
    r100kpc3d    = 0.0;
    r100kpc2d    = 0.0;

    if (strct->Type > 7)
    {
      Structure_correct_periodicity       (strct, &opt.simulation);
      Structure_shift_to_centre_of_mass   (strct);
      Structure_get_particle_radius       (strct);

      //
      //  Calculate 3D Sizes
      //
      qsort (strct->Part, strct->NumPart, sizeof(Particle), Particle_rad_compare);

      for (k = 0; k < strct->NumPart; k++)
      {
        // Total Mass
        totmass += strct->Part[k].Mass;

        // Mass Inside 3D 100 kpc
        if (strct->Part[k].Radius < 100.0)
          mass100kpc3d += strct->Part[k].Mass;
      }

      halftotmass      /= 2.0;
      halfmass100kpc3d /= 2.0;
      strct->dummyd     = 0;

      for (k = 0; k < strct->NumPart; k++)
      {
        strct->dummyd += strct->Part[k].Mass;

        if (strct->dummyd < halftotmass)      r         = strct->Part[k].Radius;
        if (strct->dummyd < halfmass100kpc3d) r100kpc3d = strct->Part[k].Radius;
      }

      //
      //  Calculate 2D Sizes
      //
      for (k = 0; k < strct->NumPart; k++)
        strct->Part[k].Pos[2] = 0.0;

      Structure_get_particle_radius (strct);
      qsort (strct->Part, strct->NumPart, sizeof(Particle), Particle_rad_compare);

      for (k = 0; k < strct->NumPart; k++)
      {
        if (strct->Part[k].Radius < 100.0)
          mass100kpc2d += strct->Part[k].Mass;
      }

      halfmass100kpc2d /= 2.0;
      strct->dummyd     = 0;

      for (k = 0; k < strct->NumPart; k++)
      {
        strct->dummyd += strct->Part[k].Mass;
        if (strct->dummyd < halfmass100kpc2d)
         r100kpc2d = strct->Part[k].Radius;
      }
      fprintf (f, "%e  ", totmass);
      fprintf (f, "%e  ", r);
      fprintf (f, "%e  ", mass100kpc3d);
      fprintf (f, "%e  ", r100kpc3d);
      fprintf (f, "%e  ", mass100kpc2d);
      fprintf (f, "%e  ", r100kpc2d);
      fprintf (f, "\n");
    }
  }
  fclose (f);

  //Catalog_free (&opt.catalog);

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
