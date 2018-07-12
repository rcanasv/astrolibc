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
  Simulation     simulation;
} Options;


void  test_usage   (int opt,  char ** argv);
int   test_options (int argc, char ** argv, Options * opt);
void  test_params  (Options * opt);

int main (int argc, char ** argv)
{
  int           i, j, k, l, n, m;
  Options       opt;

  test_options (argc, argv, &opt);
  test_params  (&opt);

  char    fname  [NAME_LENGTH];
  char    buffer [NAME_LENGTH];

  int     dummy;
  int     dummyi;
  float   dummyf;
  double  dummyd;
  char    dummys [NAME_LENGTH];

  FILE * f;

  Particle * P;

  Simulation * ramses = &opt.simulation;

  //
  //  Read Info file to get Simulation info
  //
  sprintf (fname, "%s/info_%s.txt", ramses->archive.path, ramses->archive.prefix);
  if ((f = fopen (fname, "r")) == NULL)
  {
    printf ("Couldn't open file %s\n", fname);
    exit (0);
  }
  for (i = 0; i < 7; i++)
    fgets (buffer, 100, f);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->Lbox);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->Time);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->a);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->cosmology.HubbleParam);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->cosmology.OmegaM);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->cosmology.OmegaL);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %s ", dummys, dummys, dummys);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->cosmology.OmegaB);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->unit_l);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->unit_d);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->unit_t);
  fclose (f);

  // Length
  ramses->unit_l = ramses->unit_l / 3.08e+21;                                          // in kpc
  // Time
  ramses->unit_v = ramses->unit_l / ramses->unit_t;                                    // in cm / s
  ramses->unit_v = ramses->unit_v / 100000.0;                                          // in km / s
  // Mass
  ramses->unit_m = ramses->unit_d * ramses->unit_l * ramses->unit_l * ramses->unit_l;  // in grams
  ramses->unit_m = ramses->unit_m / 1.989e+33;                                         // in solar masses
  // Box is now in kpc
  ramses->Lbox *= ramses->unit_l;

  printf ("Lbox            %e\n", ramses->Lbox);
  printf ("Time            %e\n", ramses->Time);
  printf ("Scale factor    %e\n", ramses->a);
  printf ("H0              %e\n", ramses->cosmology.HubbleParam);
  printf ("Omega_M         %e\n", ramses->cosmology.OmegaM);
  printf ("Omega_L         %e\n", ramses->cosmology.OmegaL);
  printf ("unit_l          %e\n", ramses->unit_l);
  printf ("unit_d          %e\n", ramses->unit_d);
  printf ("unit_t          %e\n", ramses->unit_t);

  //
  //  Read Particle file to get Simulation info
  //
  if (!strcmp(ramses->archive.format, "ramses_star"))  ramses->format = RAMSES_STAR;
  if (!strcmp(ramses->archive.format, "ramses"))       ramses->format = RAMSES;

  if (ramses->format == RAMSES_STAR) sprintf (fname, "%s/star_%s.out%05d", ramses->archive.name, ramses->archive.prefix, 1);
  if (ramses->format == RAMSES)      sprintf (fname, "%s/part_%s.out%05d", ramses->archive.path, ramses->archive.prefix, 1);

  f = fopen(fname,"r");
  if ((f = fopen (fname, "r")) == NULL)
  {
    printf ("Couldn't open file %s\n", fname);
    exit (0);
  }

  //!--- Header
  if (ramses->format == RAMSES_STAR)
  {
    RMSSSKIP  fread(&ramses->ncpu,     sizeof(int), 1, f);  RMSSSKIP
    RMSSSKIP  fread(&ramses->ndim,     sizeof(int), 1, f);  RMSSSKIP
    RMSSSKIP  fread(&ramses->npart,    sizeof(int), 1, f);  RMSSSKIP
    RMSSSKIP  fread(&ramses->nstarTot, sizeof(int), 1, f);  RMSSSKIP
    // Pos
    RMSSSKIP  printf("%d  ", dummy); fseek(f, dummy, SEEK_CUR); RMSSSKIP  printf("%d\n", dummy);
    RMSSSKIP  printf("%d  ", dummy); fseek(f, dummy, SEEK_CUR); RMSSSKIP  printf("%d\n", dummy);
    RMSSSKIP  printf("%d  ", dummy); fseek(f, dummy, SEEK_CUR); RMSSSKIP  printf("%d\n", dummy);
    // Vel
    RMSSSKIP  printf("%d  ", dummy); fseek(f, dummy, SEEK_CUR); RMSSSKIP  printf("%d\n", dummy);
    RMSSSKIP  printf("%d  ", dummy); fseek(f, dummy, SEEK_CUR); RMSSSKIP  printf("%d\n", dummy);
    RMSSSKIP  printf("%d  ", dummy); fseek(f, dummy, SEEK_CUR); RMSSSKIP  printf("%d\n", dummy);
    // Mass
    RMSSSKIP  printf("%d  ", dummy); fseek(f, dummy, SEEK_CUR); RMSSSKIP  printf("%d\n", dummy);
    // Id
    RMSSSKIP  printf("%d  ", dummy); fseek(f, dummy, SEEK_CUR); RMSSSKIP  printf("%d\n", dummy);
    // Level
    RMSSSKIP  printf("%d  ", dummy); fseek(f, dummy, SEEK_CUR); RMSSSKIP  printf("%d\n", dummy);
    // Birth Epoch
    RMSSSKIP  printf("%d  ", dummy); fseek(f, dummy, SEEK_CUR); RMSSSKIP  printf("%d\n", dummy);
    // Metallicity
    RMSSSKIP  printf("%d  ", dummy); fseek(f, dummy, SEEK_CUR); RMSSSKIP  printf("%d\n", dummy);
    RMSSSKIP  printf("%d  ", dummy); fseek(f, dummy, SEEK_CUR); RMSSSKIP  printf("%d\n", dummy);
    RMSSSKIP  printf("%d  ", dummy); fseek(f, dummy, SEEK_CUR); RMSSSKIP  printf("%d\n", dummy);
    RMSSSKIP  printf("%d  ", dummy); fseek(f, dummy, SEEK_CUR); RMSSSKIP  printf("%d\n", dummy);
    RMSSSKIP  printf("%d  ", dummy); fseek(f, dummy, SEEK_CUR); RMSSSKIP  printf("%d\n", dummy);
    RMSSSKIP  printf("%d  ", dummy); fseek(f, dummy, SEEK_CUR); RMSSSKIP  printf("%d\n", dummy);
    RMSSSKIP  printf("%d  ", dummy); fseek(f, dummy, SEEK_CUR); RMSSSKIP  printf("%d\n", dummy);
    RMSSSKIP  printf("%d  ", dummy); fseek(f, dummy, SEEK_CUR); RMSSSKIP  printf("%d\n", dummy);
    RMSSSKIP  printf("%d  ", dummy); fseek(f, dummy, SEEK_CUR); RMSSSKIP  printf("%d\n", dummy);
    RMSSSKIP  printf("%d  ", dummy); fseek(f, dummy, SEEK_CUR); RMSSSKIP  printf("%d\n", dummy);
    RMSSSKIP  printf("%d  ", dummy); fseek(f, dummy, SEEK_CUR); RMSSSKIP  printf("%d\n", dummy);


    printf ("Ncpu            %d\n", ramses->ncpu);
    printf ("Ndim            %d\n", ramses->ndim);
    printf ("NPart           %d\n", ramses->npart);
    printf ("NstarTot        %d\n", ramses->nstarTot);

  }

  if (ramses->format == RAMSES)
  {
    RMSSSKIP  fread(&ramses->ncpu,     sizeof(int),    1, f);  RMSSSKIP
    RMSSSKIP  fread(&ramses->ndim,     sizeof(int),    1, f);  RMSSSKIP
    RMSSSKIP  fread(&ramses->npart,    sizeof(int),    1, f);  RMSSSKIP
    RMSSSKIP  fread(&ramses->seed[0],  sizeof(int),    4, f);  RMSSSKIP
    RMSSSKIP  fread(&ramses->nstarTot, sizeof(int),    1, f);  RMSSSKIP
    RMSSSKIP  fread(&ramses->mstarTot, sizeof(double), 1, f);  RMSSSKIP
    RMSSSKIP  fread(&ramses->mstarLst, sizeof(double), 1, f);  RMSSSKIP
    RMSSSKIP  fread(&ramses->nsink,    sizeof(int),    1, f);  RMSSSKIP

    printf ("Ncpu            %d\n", ramses->ncpu);
    printf ("Ndim            %d\n", ramses->ndim);
    printf ("NPart           %d\n", ramses->npart);
    for (i = 0; i < 4; i++)
      printf ("LocalSeed[%d]    %d\n", i, ramses->seed[i]);
    printf ("NstarTot        %d\n", ramses->nstarTot);
    printf ("Mstar_tot       %g\n", ramses->mstarTot);
    printf ("Mstar_lost      %g\n", ramses->mstarLst);
    printf ("NumSink         %d\n", ramses->nsink);
 }

/*
  if ((*(part) = (Particle *) malloc (ramses->npart * sizeof(Particle))) == NULL)
  {
    printf ("Couldn't allocate memory for Particle array\n");
    exit(0);
  }

  P = *(part);

  //--- Pos
  for (i = 0; i < ramses->ndim; i++)
  {
    RMSSSKIP
    for (j = 0; j < ramses->npart; j++)
    {
      fread(&dummyd, sizeof(double), 1, f);
      P[j].Pos[i] = dummyd;
    }
    RMSSSKIP
  }

  //--- Vel
  for (i = 0; i < ramses->ndim; i++)
  {
    RMSSSKIP
    for (j = 0; j < ramses->npart; j++)
    {
      fread(&dummyd, sizeof(double), 1, f);
      P[j].Vel[i] = dummyd;
    }
    RMSSSKIP
  }

  //--- Mass
  RMSSSKIP
  for (j = 0; j < ramses->npart; j++)
  {
    fread(&dummyd, sizeof(double), 1, f);
    P[j].Mass = dummyd;
  }
  RMSSSKIP

  //--- Id
  RMSSSKIP
  for (j = 0; j < ramses->npart; j++)
  {
    fread(&dummyi, sizeof(int), 1, f);
    P[j].Id = dummyi;
  }
  RMSSSKIP

  //--- Level
  RMSSSKIP
  for (j = 0; j < ramses->npart; j++)
  {
    fread(&dummyi, sizeof(int), 1, f);
    P[j].Level = dummyi;
  }
  RMSSSKIP

  //--- Birth Epoch
  RMSSSKIP
  for (j = 0; j < ramses->npart; j++)
  {
    fread(&dummyd, sizeof(double), 1, f);
    P[j].Age = dummyd;
  }
  RMSSSKIP

  //--- Metallicity if ((STAR || SINK) && (METAL))
  //
  //  Skip for the moment.
  //
  //for (i = 0; i < 11; i++)
  //{
  //    SKIP  fread(&ramses_met[i][0], sizeof(double), npart, f);  SKIP
  //}
  fclose (f);

  //
  // Convert to human readable units
  //
  for (i = 0; i < ramses->npart; i++)
  {
    P[i].Pos[0] *= ramses->unit_l;
    P[i].Pos[1] *= ramses->unit_l;
    P[i].Pos[2] *= ramses->unit_l;

    P[i].Vel[0] *= ramses->unit_v;
    P[i].Vel[1] *= ramses->unit_v;
    P[i].Vel[2] *= ramses->unit_v;

    P[i].Mass   *= ramses->unit_m;
  }
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

  // Name
  fscanf (opt->param.file, "%s  %s  %s  %s  %d", namebuff, prefixbuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->simulation.archive, namebuff);
  Archive_prefix (&opt->simulation.archive, prefixbuff);
  Archive_format (&opt->simulation.archive, frmtbuff);
  Archive_path   (&opt->simulation.archive, pathbuff);
  Archive_nfiles (&opt->simulation.archive, nflsbuff);

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

  while ((myopt = getopt_long (argc, argv, "p:ftxvh", lopts, &index)) != -1)
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
    printf ("  read star data ramses                                                  \n");
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
