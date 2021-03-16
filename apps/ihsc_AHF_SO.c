/*
 *
 *  \file    ihsc_AHF_SO.c
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
  int            iVerbose;
  int            iFraction;
  int            iExtract;
  Archive        param;
  Archive        output;
  Archive        clean;
  Catalog        stf;
  Catalog        ahf;
  Simulation     sim;
} Options;


void  ihsc_usage   (int opt,  char ** argv);
int   ihsc_options (int argc, char ** argv, Options * opt);
void  ihsc_params  (Options * opt);


int main (int argc, char ** argv)
{
  // 1. Read AHF clean halos ASCII
  // 2. Load AHF Catalog
  // 3. Load AHF Particle list inside
  //    Sort IDs
  // 4. Load Particles from Simulation
  // 5. Load Extended Output
  // 6. Sort Particles by ID

  int           i, j, k, l, n, m;
  Options       opt;
  Structure   * strct1;
  Structure   * strct2;
  Structure   * strct3;
  Structure     tmpstrct;
  int           numpart;
  int         * strct_to_get;
  Particle    * P;
  double        fihsc, mstot;
  int           top;

  ihsc_options (argc, argv, &opt);
  ihsc_params  (&opt);

  Simulation_init         (&opt.sim);
  Catalog_init            (&opt.stf);
  Catalog_load_properties (&opt.stf);
  Catalog_init            (&opt.ahf);
  Catalog_load_properties (&opt.ahf);


  // --------------------------------------------------- //
  //            IHSC SO FROM AHF CLEAN HALOS             //
  // --------------------------------------------------- //
  char         fname   [NAME_LENGTH];
  char         longbuffer [LONG_LENGTH];
  int          ncleans;
  FILE       * f;
  Structure  * strct_clean;

  // 1. Read AHF clean halos ASCII
  sprintf (fname, "%s/%s", opt.clean.path, opt.clean.prefix);
  if ((f = fopen (fname, "r")) == NULL)
  {
    printf ("ERROR: Cannot open file  %s\n", fname);
    exit (0);
  }

  ncleans = 0;
  fgets  (longbuffer, NAME_LENGTH, f);  // Header
  while ( fgets(longbuffer, NAME_LENGTH, f) != NULL) ncleans++;
  rewind (f);
  if ((strct_clean = (Structure *) malloc (ncleans*sizeof(Structure))) == NULL)
  {
    printf ("Cannot allocate memory for Structures\n");
    exit (0);
  }

  fgets  (longbuffer, NAME_LENGTH, f);  // Header
  for (i = 0; i < ncleans; i++)
  {
    fgets  (longbuffer, NAME_LENGTH, f);
    sscanf (longbuffer, "%d  %d  %lf  %lf  %lf  %lf  %lf",              \
                        &strct_clean[i].HostID, &strct_clean[i].ID,     \
                        &strct_clean[i].Mvir,   &strct_clean[i].Rvir,   \
                        &strct_clean[i].Pos[0], &strct_clean[i].Pos[1], \
                        &strct_clean[i].Pos[2]);
  }
  fclose (f);


  for (i = 0; i < 10; i++)
    printf ("%d  %d  %lf  %lf\n", strct_clean[i].HostID, strct_clean[i].ID, strct_clean[i].Mvir, strct_clean[i].Rvir);

  // 2. Load AHF Catalog
  // 3. Load AHF Particle list inside
  //    Sort IDs
  // 4. Load Particles from Simulation
  // 5. Load Extended Output
  // 6. Sort Particles by ID

  // --------------------------------------------------- //



  // --------------------------------------------------- //
  //                    Free catalogues                  //
  // --------------------------------------------------- //
  Catalog_free (&opt.stf);
  Catalog_free (&opt.ahf);
  // --------------------------------------------------- //


  return (0);
}



//  Parameters
void ihsc_params (Options * opt)
{
  int   i;
  int   dummy;
  char  buffer     [NAME_LENGTH];
  char  prefixbuff [NAME_LENGTH];
  char  namebuff   [NAME_LENGTH];
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

  // Output
  fscanf (opt->param.file, "%s  %s  %s  %s  %d", prefixbuff, namebuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->output, namebuff);
  Archive_prefix (&opt->output, prefixbuff);
  Archive_format (&opt->output, frmtbuff);
  Archive_path   (&opt->output, pathbuff);
  Archive_nfiles (&opt->output, nflsbuff);

  // Catalog STF
  fscanf (opt->param.file, "%s  %s  %s  %s  %d", prefixbuff, namebuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->stf.archive, namebuff);
  Archive_prefix (&opt->stf.archive, prefixbuff);
  Archive_format (&opt->stf.archive, frmtbuff);
  Archive_path   (&opt->stf.archive, pathbuff);
  Archive_nfiles (&opt->stf.archive, nflsbuff);

  // Catalog  AHF
  fscanf (opt->param.file, "%s  %s  %s  %s  %d", prefixbuff, namebuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->ahf.archive, namebuff);
  Archive_prefix (&opt->ahf.archive, prefixbuff);
  Archive_format (&opt->ahf.archive, frmtbuff);
  Archive_path   (&opt->ahf.archive, pathbuff);
  Archive_nfiles (&opt->ahf.archive, nflsbuff);

  // Simulation
  fscanf (opt->param.file, "%s  %s  %s  %s  %d", prefixbuff, namebuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->sim.archive, namebuff);
  Archive_prefix (&opt->sim.archive, prefixbuff);
  Archive_format (&opt->sim.archive, frmtbuff);
  Archive_path   (&opt->sim.archive, pathbuff);
  Archive_nfiles (&opt->sim.archive, nflsbuff);

  // Close
  fclose (opt->param.file);
}


//
//  Options
//
int ihsc_options (int argc, char ** argv, Options * opt)
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
    {"extract",   0, NULL, 'x'},
    {0,           0, NULL, 0}
  };

  opt->iVerbose  = 0;
  opt->iExtract  = 0;

  while ((myopt = getopt_long (argc, argv, "p:ftxvhs", lopts, &index)) != -1)
  {
    switch (myopt)
    {
      case 'p':
      	strcpy (opt->param.name, optarg);
        flag++;
        break;

      case 'x':
      	opt->iExtract = 1;
      	break;

      case 'h':
      	ihsc_usage (0, argv);
        break;

      default:
      	ihsc_usage (1, argv);
    }
  }

  if (flag == 0)
    ihsc_usage (1, argv);
}


//
//  Usage
//
void ihsc_usage (int opt, char ** argv)
{
  if (opt == 0)
  {
    printf ("                                                                         \n");
    printf ("  ihsc_AHF_SO                                                            \n");
    printf ("                                                                         \n");
    printf ("  Author:            Rodrigo Can\\~as                                    \n");
    printf ("                                                                         \n");
    printf ("  Last edition:      16 - March - 2021                                      \n");
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
