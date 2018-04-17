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
  int            nstruct;
  int            ilist;
  int          * id;
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
  int          i, j, k;
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


  strct_to_get = (int *) malloc ((opt.catalog.nstruct+1) * sizeof(int));
  for (i = 1; i <= opt.catalog.nstruct; i++)
    strct_to_get[i] = 0;

  for (i = 0; i < opt.nstruct; i++)
    strct_to_get[opt.id[i]] = 1;

  strct = (Structure **) malloc (opt.nstruct*3*sizeof(Structure *));


  //
  // Merge central with IHSM
  //
  int         totpart;
  Structure   tmpstrct;
  int         n;

  for (i = 1, n = 0; i <= opt.catalog.nstruct; i++)
  {
    strct1 = &opt.catalog.strctProps[i];

    if ((strct1->Central == 1) && (strct_to_get[i] == 1) && (strct1->HostID > 0))
    {
      strct2 = &opt.catalog.strctProps[strct1->HostID];
      totpart = strct1->NumPart + strct2->NumPart;
      tmpstrct.NumPart = totpart;
      tmpstrct.Part = (Particle *) malloc (totpart * (sizeof(Particle)));

      for (j = 0, k = 0; j < strct1->NumPart; j++, k++)
        Particle_copy (&strct1->Part[j], &tmpstrct.Part[k]);

      for (j = 0; j < strct2->NumPart; j++, k++)
        Particle_copy (&strct2->Part[j], &tmpstrct.Part[k]);

      strct2->Pos[0] = strct1->Pos[0];
      strct2->Pos[1] = strct1->Pos[1];
      strct2->Pos[2] = strct1->Pos[2];

      tmpstrct.Pos[0] = strct1->Pos[0];
      tmpstrct.Pos[1] = strct1->Pos[1];
      tmpstrct.Pos[2] = strct1->Pos[2];

      strct1->dummyi   = strct1->ID;
      strct2->dummyi   = strct1->ID;
      tmpstrct.dummyi = strct1->ID;

      Structure_correct_periodicity (strct1,    &opt.simulation);
      Structure_correct_periodicity (strct2,    &opt.simulation);
      Structure_correct_periodicity (&tmpstrct, &opt.simulation);

      strct[n++] = strct1;
      strct[n++] = strct2;
      strct[n++] = &tmpstrct;
    }
  }

  //
  //  Calculate Surface density and create files
  //
  FILE      * f;
  double    * r     = NULL;
  double    * Sigma = NULL;
  char        buffer [NAME_LENGTH];
  int         bob;
  int         nbins = 20;
  double      rmin = 0.0;
  double      rmax = 200.0;


  for (i = 0, n = 0; i < opt.nstruct; i++)
  {
    // Central
    strct1 = strct[n++];
    Structure_calculate_surface_density (strct1, NULL, rmin, rmax, nbins, &r, &Sigma);
    sprintf (buffer, "surface_density_%07d_central.dat", strct1->dummyi);
    f = fopen (buffer, "w");
    for (j = 0; j < nbins; j++)
      fprintf (f, "%e  %e\n", r[j], Sigma[j]);
    fclose (f);

    // IHSC
    strct1 = strct[n++];
    Structure_calculate_surface_density (strct1, NULL, rmin, rmax, nbins, &r, &Sigma);
    sprintf (buffer, "surface_density_%07d_ihsc.dat", strct1->dummyi);
    f = fopen (buffer, "w");
    for (j = 0; j < nbins; j++)
      fprintf (f, "%e  %e\n", r[j], Sigma[j]);
    fclose (f);

    // Both
    strct1 = strct[n++];
    Structure_calculate_surface_density (strct1, NULL, rmin, rmax, nbins, &r, &Sigma);
    sprintf (buffer, "surface_density_%07d_both.dat", strct1->dummyi);
    f = fopen (buffer, "w");
    for (j = 0; j < nbins; j++)
      fprintf (f, "%e  %e\n", r[j], Sigma[j]);
    fclose (f);
  }

  //
  // Free Memory
  //
  for (i = 0, n = 0; i < opt.nstruct; i++)
  {
    n = n + 2;
    free (strct[n++]->Part);
  }
  free (strct_to_get);
  free (strct);
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

 if (opt->ilist)
 {
    opt->list.file = fopen (opt->list.name, "r");
    if (opt->list.file == NULL)
    {
      printf ("Couldn't open file  %s\n", opt->list.name);
      printf ("Exiting...\n");
      exit (0);
    }
    fscanf (opt->list.file, "%d", &opt->nstruct);
    opt->id = (int *) malloc (opt->nstruct * sizeof(int));
    for (i = 0; i < opt->nstruct; i++)
      fscanf (opt->list.file, "%d", &opt->id[i]);
    fclose (opt->list.file);
  }
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
    {"list",    2, NULL, 'l'},
    {0,         0, NULL, 0}
  };

  opt->ilist = 0;

  while ((myopt = getopt_long (argc, argv, "i:p:l:vh", lopts, &index)) != -1)
  {
    switch (myopt)
    {
      case 'p':
      	strcpy (opt->param.name, optarg);
        flag++;
        break;

      case 'i':
        opt->nstruct = 1;
      	opt->id = (int *) malloc (opt->nstruct * sizeof(int));
        opt->id[0] = atoi(optarg);
        flag++;
        break;

      case 'l':
        opt->ilist = 1;
        strcpy (opt->list.name, optarg);
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
    printf ("  Last edition:      16 - 04 - 2018                                      \n");
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
