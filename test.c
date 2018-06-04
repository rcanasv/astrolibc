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
  int            nsnap;
  int            ntrees;
  Archive        param;
  Archive        output;
  Catalog      * catalog;
  Archive      * tree;
  Simulation   * simulation;
} Options;


void  test_usage   (int opt,  char ** argv);
int   test_options (int argc, char ** argv, Options * opt);
void  test_params  (Options * opt);


int main (int argc, char ** argv)
{
  int          i, j, k, l;
  Options      opt;
  Structure   * strct1;
  Structure   * strct2;
  Structure   * strct3;


  test_options (argc, argv, &opt);
  test_params  (&opt);


  //
  //  Load catalogs and simulation details
  // for Horizon-AGN only info_XXXX is needed
  //
  // Here I am assuming progenitor trees with a single
  // snapshot connection AND only two snapshots per tree file
  //
  // This means  opt.tree[n] has the progenitors of
  // opt.catalos[n] in opt.catalog[n+1] and so on
  //
  for (i = 0; i < opt.nsnap; i++)
  {
    Simulation_init                 (&opt.simulation[i]);
    Catalog_init                    (&opt.catalog[i]);
    Catalog_load_properties         (&opt.catalog[i]);
    //Catalog_get_particle_properties (&opt.catalog[i], &opt.simulation[i]);
    Catalog_fill_isolated           (&opt.catalog[i]);
    if (i < (opt.nsnap - 1))
      stf_read_treefrog (&opt.tree[i], &opt.catalog[i]);
  }

  //
  // Tag central galaxy and add stellar mass
  //
  for (i = 0; i < opt.nsnap; i++)
  {
    for (j = 1; j <= opt.catalog[i].nstruct; j++)
    {
      strct1 = &opt.catalog[i].strctProps[j];
      if (strct1->Central == 1)
      {
        strct2 = &opt.catalog[i].strctProps[strct1->HostID];
        strct1->dummyd = strct1->TotMass + strct2->TotMass;
        strct2->dummyi = strct1->ID;
      }
    }
  }

  for (i = 0; i < opt.nsnap; i++)
  {
    for (j = 1; j <= opt.catalog[i].nstruct; j++)
    {
      strct1 = &opt.catalog[i].strctProps[j];
      if (strct1->Central != 1 && strct1->Type > 7)
      {
        strct2 = &opt.catalog[i].strctProps[strct1->HostID];
        strct3 = &opt.catalog[i].strctProps[strct2->dummyi];
        strct3->TotMass += strct1->TotMass;
      }
    }
  }

  //
  // Write file with Diffuse stellar mass
  //
  FILE * f;
  char buffer [NAME_LENGTH];
  for (i = 0; i < opt.nsnap; i++)
  {
    sprintf (buffer, "%s.ihsc", opt.catalog[i].archive.prefix);
    f = fopen (buffer, "w");
    for (j = 1; j <= opt.catalog[i].nstruct; j++)
    {
      strct1 = &opt.catalog[i].strctProps[j];
      if (strct1->Type == 7)
      {
        strct2 = &opt.catalog[i].strctProps[strct1->dummyi];
        fprintf (f, "%e  ", strct1->TotMass);
        fprintf (f, "%e  ", strct2->TotMass);
        fprintf (f, "%e  ", strct2->dummyd);
        fprintf (f, "%5d ", strct1->NumSubs);
        fprintf (f, "\n");
      }
    }
    fclose (f);
 }


  //
  // Create `evolutionary tracks'
  //


  //
  // Free catalogues
  //
  for (i = 0; i < opt.nsnap; i++)
    Catalog_free (&opt.catalog[i]);

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

  fscanf (opt->param.file, "%d", &opt->nsnap);
  opt->ntrees = opt->nsnap - 1;

  // Catalogues
  for (i = 0; i < opt->nsnap; i++)
  {
    fscanf (opt->param.file, "%s", buffer);  Archive_name   (&opt->catalog[i].archive, buffer);
                                             Archive_prefix (&opt->catalog[i].archive, buffer);
    fscanf (opt->param.file, "%s", buffer);  Archive_format (&opt->catalog[i].archive, buffer);
    fscanf (opt->param.file, "%s", buffer);  Archive_path   (&opt->catalog[i].archive, buffer);
    fscanf (opt->param.file, "%d", &dummy);  Archive_nfiles (&opt->catalog[i].archive, dummy);
  }

  // Trees
  for (i = 0; i < opt->ntrees; i++)
  {
    fscanf (opt->param.file, "%s", buffer);  Archive_name   (&opt->tree[i], buffer);
                                             Archive_prefix (&opt->tree[i], buffer);
    fscanf (opt->param.file, "%s", buffer);  Archive_format (&opt->tree[i], buffer);
    fscanf (opt->param.file, "%s", buffer);  Archive_path   (&opt->tree[i], buffer);
    fscanf (opt->param.file, "%d", &dummy);  Archive_nfiles (&opt->tree[i], dummy);
  }

  // Simulation
  for (i = 0; i < opt->nsnap; i++)
  {
    fscanf (opt->param.file, "%s", buffer);  Archive_name   (&opt->simulation[i].archive, buffer);
                                             Archive_prefix (&opt->simulation[i].archive, buffer);
    fscanf (opt->param.file, "%s", buffer);  Archive_format (&opt->simulation[i].archive, buffer);
    fscanf (opt->param.file, "%s", buffer);  Archive_path   (&opt->simulation[i].archive, buffer);
    fscanf (opt->param.file, "%d", &dummy);  Archive_nfiles (&opt->simulation[i].archive, dummy);
  }

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
    {0,         0, NULL, 0}
  };

  while ((myopt = getopt_long (argc, argv, "p:vh", lopts, &index)) != -1)
  {
    switch (myopt)
    {
      case 'p':
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
    printf ("  ihsc                                                                   \n");
    printf ("                                                                         \n");
    printf ("  Author:            Rodrigo Can\\~as                                    \n");
    printf ("                                                                         \n");
    printf ("  Last edition:      04 - 06 - 2018                                      \n");
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
