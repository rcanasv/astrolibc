/*
 *
 *  \file    get_structure.c
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


void  get_structure_usage   (int opt,  char ** argv);
int   get_structure_options (int argc, char ** argv, Options * opt);
void  get_structure_params  (Options * opt);


int main (int argc, char ** argv)
{
  int          i, j, k;
  Options      opt;
  Structure *  strct;
  gheader      header;
  int       *  strct_to_get;


  get_structure_options (argc, argv, &opt);
  get_structure_params  (&opt);


  //
  //  Load catalogs
  //
  Simulation_init (&opt.simulation);
  Catalog_init (&opt.catalog);
  Catalog_load_properties (&opt.catalog);

  // 3995614986
  // 4273503266
  // 3840676851

  //
  //  1.  Get Individual structure
  //  2.  Get Several structures a file for each one
  //  3.  Get Several structures into a single file
  //  4.  Get everything inside 3DFOF
  //  2.  Get everything inside 6DFOF
  //

  strct_to_get = (int *) malloc ((opt.catalog.nstruct+1) * sizeof(int));
  for (i = 1; i <= opt.catalog.nstruct; i++)
    strct_to_get[i] = 0;

  for (i = 0; i < opt.nstruct; i++)
  {
//    printf ("%d\n", opt.id[i]);
    strct_to_get[opt.id[i]] = 1;
  }

/*
  for (i = 1; i <= opt.catalog.nstruct; i++)
    if (opt.catalog.strctProps[i].DirectHostID == opt.id[0])
      strct_to_get[i] = 1;
*/


  for (i = 1; i <= opt.catalog.nstruct; i++)
    if (opt.catalog.strctProps[i].HostID == opt.id[0])
      strct_to_get[i] = 1;


  Structure_get_particle_properties (&opt.catalog, &opt.simulation, strct_to_get);



  Particle * P;
  int        numpart = 0;

  for (i = 1; i <= opt.catalog.nstruct; i++)
  {
    if (strct_to_get[i] == 1)
      numpart += opt.catalog.strctProps[i].NumPart;
  }

  P = (Particle *) malloc (numpart * sizeof(Particle));

  k = 0;

  for (i = 1; i <= opt.catalog.nstruct; i++)
  {
    if (strct_to_get[i] == 1)
      for (j = 0; j < opt.catalog.strctProps[i].NumPart; j++, k++)
      {
        Particle_copy (&opt.catalog.strctProps[i].Part[j], &P[k]);
        P[k].Type = 1;
      }
  }

  gadget_write_snapshot (P, numpart, &header, &opt.output);

  free (P);

/*

  for (i = 1, k = 0; i <= opt.catalog.nstruct; i++)
  {
    if (strct_to_get[i] == 1)
    {
      strct = &opt.catalog.strctProps[i];
      for (j = 0; j < strct->NumPart; j++)
        strct->Part[j].Type = 1;
      sprintf (opt.output.name, "%s.gdt_%03d", opt.output.prefix, k);

      gadget_write_snapshot (strct->Part, strct->NumPart, &header, &opt.output);
      k++;
    }
  }
*/


  free (strct_to_get);
  Catalog_free (&opt.catalog);

  return (0);
}


//
//  Parameters
//
void get_structure_params (Options * opt)
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
int get_structure_options (int argc, char ** argv, Options * opt)
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
      	get_structure_usage (0, argv);
        break;

      default:
      	get_structure_usage (1, argv);
    }
  }

  if (flag == 0)
    get_structure_usage (1, argv);

}


//
//  Usage
//
void get_structure_usage (int opt, char ** argv)
{
  if (opt == 0)
  {
    printf ("                                                                         \n");
    printf ("  get_structure.c                                                        \n");
    printf ("                                                                         \n");
    printf ("  Author:            Rodrigo Can\\~as                                    \n");
    printf ("                                                                         \n");
    printf ("  Last edition:      29 - 03 - 2018                                      \n");
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
