/*
 *
 *  \file    get_structure.c
 *  \brief
 *
 *
 */


#include "../src/base.h"
#include "../src/typedef.h"
#include "../src/archive.h"
#include "../src/catalog.h"
#include "../src/simulation.h"


typedef struct Options
{
  int            verbose;
  Archive        param;
  Simulation     simulation;
} Options;


void  check_gagetio_usage   (int opt,  char ** argv);
int   check_gagetio_options (int argc, char ** argv, Options * opt);
void  check_gagetio_params  (Options * opt);


int main (int argc, char ** argv)
{
  int          i, j, k;
  int          offset;
  Options      opt;
  gheader      header;
  Particle     *P;

  check_gagetio_options (argc, argv, &opt);
  check_gagetio_params  (&opt);


  // HEADER
  Simulation_init (&opt.simulation);
  printf ("Npart\n");
  for (i = 0; i < 6; i++)
    printf ("   Type %d     %d\n", i, opt.simulation.NpartThisFile[i]);
  printf ("NpartTot\n");
  for (i = 0; i < 6; i++)
    printf ("   Type %d     %d\n", i, opt.simulation.NpartTot[i]);
  printf ("MassTable\n");
  for (i = 0; i < 6; i++)
    printf ("   Type %d     %f\n", i, opt.simulation.MassTable[i]);
  printf("NumFiles    %d\n", opt.simulation.NfilesPerSnapshot);
  printf("Cooling     %d\n", opt.simulation.Cooling);
  printf("Feedback    %d\n", opt.simulation.Feedback);
  printf("SFR         %d\n", opt.simulation.SFR);
  printf("h           %f\n", opt.simulation.h);
  printf("Omega_m     %f\n", opt.simulation.cosmology.OmegaM);
  printf("Omega_m     %f\n", opt.simulation.cosmology.OmegaL);
  printf("a           %f\n", opt.simulation.a);
  printf("z           %f\n", opt.simulation.z);

  // PART DATA
  Simulation_load_particles (&opt.simulation, 0, &P);
  for (k = 0, offset = 0; k < 6; k++)
  {
    printf("PartType %d\n", k);
    for (i = 0; i < 10; i++)
    {
      printf("    ");
      printf("%e  ", P[i+offset].Pos[0]);
      printf("%e  ", P[i+offset].Pos[1]);
      printf("%e  ", P[i+offset].Pos[2]);
      printf("%e  ", P[i+offset].Vel[0]);
      printf("%e  ", P[i+offset].Vel[1]);
      printf("%e  ", P[i+offset].Vel[2]);
      printf("%e  ", P[i+offset].Mass);
      printf("\n");
    }
    offset += opt.simulation.NpartThisFile[k];
  }
  free(P);

  return (0);
}


//
//  Parameters
//
void check_gagetio_params (Options * opt)
{
  int   i;
  int   dummy;
  char  buffer   [NAME_LENGTH];
  char  namebuff [NAME_LENGTH];
  char  prfxbuff [NAME_LENGTH];
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

  // Simulation
  fscanf (opt->param.file, "%s  %s  %s  %s  %d", prfxbuff,namebuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->simulation.archive, namebuff);
  Archive_prefix (&opt->simulation.archive, prfxbuff);
  Archive_format (&opt->simulation.archive, frmtbuff);
  Archive_path   (&opt->simulation.archive, pathbuff);
  Archive_nfiles (&opt->simulation.archive, nflsbuff);

  // Close
  fclose (opt->param.file);
}


//
//  Options
//
int check_gagetio_options (int argc, char ** argv, Options * opt)
{
  int   myopt;
  int   index;
  int   flag = 0;

  extern char * optarg;
  extern int    opterr;
  extern int    optopt;

  struct option lopts[] = {
    {"help",         0, NULL, 'h'},
    {"param",        1, NULL, 'p'},
    {0,              0, NULL,   0}
  };

  while ((myopt = getopt_long (argc, argv, "p:h", lopts, &index)) != -1)
  {
    switch (myopt)
    {
      case 'p':
      	strcpy (opt->param.name, optarg);
        flag++;
        break;

      case 'h':
      	check_gagetio_usage (0, argv);
        break;

      default:
      	check_gagetio_usage (1, argv);
    }
  }

  if (flag == 0)
    check_gagetio_usage (1, argv);

}


//
//  Usage
//
void check_gagetio_usage (int opt, char ** argv)
{
  if (opt == 0)
  {
    printf ("                                                                         \n");
    printf ("  check_gagetio.c                                                        \n");
    printf ("  Author:            Rodrigo Can\\~as                                    \n");
    printf ("  Last edition:      11 - 02 - 2021                                      \n");
    printf ("                                                                         \n");
    printf ("                                                                         \n");
    printf ("  Usage:             %s [Option] [Parameter [argument]] ...\n",      argv[0]);
    printf ("                                                                         \n");
    printf ("  Parameters:                                                            \n");
    printf ("                     -p    --param    [string]   Name of input file      \n");
    printf ("                                                                         \n");
    printf ("  Options:                                                               \n");
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
