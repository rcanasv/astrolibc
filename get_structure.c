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


  test_options (argc, argv, &opt);
  test_params  (&opt);


  //
  //  Load catalogs
  //
  //printf ("simulation init\n");
  Simulation_init (&opt.simulation);
  Catalog_init (&opt.catalog);
  Catalog_load_properties (&opt.catalog);

  Catalog_get_particle_properties (&opt.catalog, &opt.simulation);

  3995614986
  4273503266
  3840676851

  strcpy (stfprefix, argv[1]);
  strcpy (directory, argv[2]);
  strcpy (galprefix, argv[3]);
  strcpy (outprefix, argv[4]);
  strcpy (format,    argv[5]);
  NUMFILES = atoi (argv[6]);
  int ID = atoi (argv[7]);


  if ( (strcmp(format, "gadget") != 0) && (strcmp(format,"ramses") != 0) )
  {
    printf ("Incorrect specified format  %s, exiting...\n", format);
    exit (0);
  }

  int * filestoread = (int *) malloc (NUMFILES        * sizeof(int));
  int * getgal      = (int *) malloc (snap1.nstruct+1 * sizeof(int));

  int acum        = 0;
  int nparttot    = 0;
  int ninextended = 0;
  int dummystruct = 0;
  int dummyindex  = 0;

  double cm [3];
  cm[0] = snap1.strctProps[ID].Pos[0];
  cm[1] = snap1.strctProps[ID].Pos[1];
  cm[2] = snap1.strctProps[ID].Pos[2];

  ID = snap1.strctProps[ID].HostID;

  for (i = 0; i < NUMFILES; i++)
    filestoread[i] = 0;

  for (i = 1; i <= snap1.nstruct; i++)
  {
    getgal[i] = 0;
    if (snap1.strctProps[i].HostID == ID)
      getgal[i] = 1;
  }
  getgal[ID] = 1;

  get_files_to_read_filesofgroup (stfprefix, snap1.nstruct, getgal, filestoread);

  for (i = 1; i <= snap1.nstruct; i++)
    if (getgal[i])
      nparttot += snap1.strctProps[i].NumPart;

  if (!(p = malloc (nparttot * sizeof(struct pgadget))))
  {
    printf ("Cannot allocate memory for particle information for structure %d\n", i);
    printf ("Exiting\n");
    exit (0);
  }

  for (i = 0; i < NUMFILES; i++)
  {
    if (filestoread[i])
    {
      printf ("Reading file %d\n", i);
      ninextended = load_stf_extended_output (stfprefix, i);
      if (ninextended)
      {
        read_ramses_snapshot (directory, galprefix, i);

        for (j = 0; j < ninextended; j++)
        {
          dummystruct = extended_IdStruct[j];
          dummyindex  = extended_oIndex[j];

          if (getgal[dummystruct])
          {
            p[acum].Pos[0] = ramses_pos[0][dummyindex];
            p[acum].Pos[1] = ramses_pos[1][dummyindex];
            p[acum].Pos[2] = ramses_pos[2][dummyindex];
            p[acum].Vel[0] = ramses_vel[0][dummyindex];
            p[acum].Vel[1] = ramses_vel[1][dummyindex];
            p[acum].Vel[2] = ramses_vel[2][dummyindex];
            p[acum].Mass   = ramses_mass  [dummyindex];
            p[acum].Id     = ramses_id    [dummyindex];
            p[acum].Age    = ramses_age   [dummyindex];
            p[acum].Type   = 1;
            acum++;
          }
        }
        free_extended_arrays ();
        free_ramses_arrays ();
      }
    }
  }


  gadget_write_snapshot(strct->Part, strct->NumPart, opt.simulation, outprefix, nfiles);

  free (filestoread);
  free (getgal);

















  //Catalog_free (&opt.catalog);

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
