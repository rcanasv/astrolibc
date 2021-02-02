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
  int            iVerbose;
  int            iFraction;
  int            iTrack;
  int            iExtract;
  int            nsnap;
  int            ntrees;
  Archive        param;
  Archive        output;
  Catalog      * catalog;
  Archive      * tree;
  Simulation   * simulation;
  int            top;
  double         mass_min;
  double         mass_max;
  double         frac_min;
  double         frac_max;
} Options;


void  test_usage   (int opt,  char ** argv);
int   test_options (int argc, char ** argv, Options * opt);
void  test_params  (Options * opt);


int main (int argc, char ** argv)
{
  int           i, j, k, l, n, m;
  Options       opt;
  Structure   * strct1;
  Structure   * strct2;
  Structure   * strct3;
  Structure   * sorted;
  Structure     tmpstrct;
  int           numpart;
  int        ** strct_to_get;
  Particle    * P;
  gheader       header;
  double        fihsc, mstot;
  double        mass_min, mass_max;
  double        frac_min, frac_max;
  int           top;

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
    if (opt.iTrack)
      if (i < (opt.nsnap - 1))
        stf_read_treefrog (&opt.tree[i], &opt.catalog[i]);
  }

  FILE * f;
  char buffer [NAME_LENGTH];

  sprintf (buffer, "%s", opt.output.prefix);
  f = fopen (buffer, "w");
  for (j = 1; j <= opt.catalog[0].nstruct; j++)
  {
    strct1 = &opt.catalog[0].strctProps[j];
printf ("%d  %d\n", strct1->ID, strct1->NumMatch);
    if (strct1->Type > 7 && strct1->NumMatch)
    {
      strct2 = &opt.catalog[1].strctProps[strct1->MatchIDs[0]];  
      fprintf (f, "%d  ", strct1->ID);       // stf  ID
      fprintf (f, "%d  ", strct2->ID);       // hmkr ID
      fprintf (f, "%e  ", strct1->TotMass);  // stf  Mass
      fprintf (f, "%e  ", strct2->TotMass);  // hmkr Mass
      fprintf (f, "\n");
    }
  }
  fclose (f);
  // --------------------------------------------------- //



  // --------------------------------------------------- //
  //
  // Free catalogues
  //
  for (i = 0; i < opt.nsnap; i++)
    Catalog_free (&opt.catalog[i]);
  free (opt.catalog);
  free (opt.simulation);
  if (opt.iTrack)
    free (opt.tree);
  // --------------------------------------------------- //


  return (0);
}




//
//  Parameters
//
void test_params (Options * opt)
{
  int   i;
  int   dummy;
  char  buffer   [NAME_LENGTH];
  char  prefixbuff [NAME_LENGTH];
  char  namebuff [NAME_LENGTH];
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

  // Output
  fscanf (opt->param.file, "%s  %s  %s  %d", namebuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->output, namebuff);
  Archive_prefix (&opt->output, namebuff);
  Archive_format (&opt->output, frmtbuff);
  Archive_path   (&opt->output, pathbuff);
  Archive_nfiles (&opt->output, nflsbuff);


  // Update number of trees
  fscanf (opt->param.file, "%d", &opt->nsnap);
  opt->ntrees = opt->nsnap - 1;

  opt->catalog    = (Catalog    *) malloc (opt->nsnap * sizeof(Catalog));
  opt->simulation = (Simulation *) malloc (opt->nsnap * sizeof(Simulation));
  if (opt->iTrack)
    opt->tree     = (Archive    *) malloc (opt->nsnap * sizeof(Archive));


  // Catalogues
  for (i = 0; i < opt->nsnap; i++)
  {
    fscanf (opt->param.file, "%s  %s  %s  %s  %d", prefixbuff, namebuff, frmtbuff, pathbuff, &nflsbuff);
    Archive_name   (&opt->catalog[i].archive, namebuff);
    Archive_prefix (&opt->catalog[i].archive, prefixbuff);
    Archive_format (&opt->catalog[i].archive, frmtbuff);
    Archive_path   (&opt->catalog[i].archive, pathbuff);
    Archive_nfiles (&opt->catalog[i].archive, nflsbuff);
  }


  // Simulation
  for (i = 0; i < opt->nsnap; i++)
  {
    fscanf (opt->param.file, "%s  %s  %s  %s  %d", prefixbuff, namebuff, frmtbuff, pathbuff, &nflsbuff);
    Archive_name   (&opt->simulation[i].archive, namebuff);
    Archive_prefix (&opt->simulation[i].archive, prefixbuff);
    Archive_format (&opt->simulation[i].archive, frmtbuff);
    Archive_path   (&opt->simulation[i].archive, pathbuff);
    Archive_nfiles (&opt->simulation[i].archive, nflsbuff);
  }

  // Trees
  if (opt->iTrack)
  {
    // Trees
    for (i = 0; i < opt->ntrees; i++)
    {
      fscanf (opt->param.file, "%s  %s  %s  %s  %d", prefixbuff, namebuff, frmtbuff, pathbuff, &nflsbuff);
      Archive_name   (&opt->tree[i], namebuff);
      Archive_prefix (&opt->tree[i], prefixbuff);
      Archive_format (&opt->tree[i], frmtbuff);
      Archive_path   (&opt->tree[i], pathbuff);
      Archive_nfiles (&opt->tree[i], nflsbuff);
    }
  }


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
    {"fraction",  0, NULL, 'f'},
    {"track",     0, NULL, 't'},
    {"extract",   0, NULL, 'x'},
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

      case 'f':
      	opt->iFraction = 1;
      	break;

      case 't':
      	opt->iTrack = 1;
      	break;

      case 'x':
      	opt->iExtract = 1;
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
    printf ("  Last edition:      25 - 06 - 2018                                      \n");
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

/*
void ihsc_prog_tree (Catalog * ctlgs, int pihscid, int plevel, int maxlvls, char * mainbuff, char * brnchbuff, FILE * file)
{
  int          i, j;
  int          level;
  char         lvl_mainbuff   [3000];
  char         lvl_brnchbuff  [3000];
  char         strctbuff      [3000];
  Structure  * ctrl;
  Structure  * ihsc;
  Structure  * strct;

  level = plevel + 1;

  ihsc = &ctlgs[plevel].strctProps[pihscid];

  if (plevel == 0)
  {
    sprintf (mainbuff,   "% 7d  % 7d  ", ihsc->ID, ihsc->ID);
    sprintf (brnchbuff,  "% 7d  % 7d  ", ihsc->ID, -1);
  }

  sprintf (lvl_mainbuff,  "%s", mainbuff);
  sprintf (lvl_brnchbuff, "%s", brnchbuff);

  for (i = 0, j = 0; (i < ihsc->NumMatch && j < 4); i++)
  {
    strct = &ctlgs[level].strctProps[ihsc->MatchIDs[i]];
    ctrl  = &ctlgs[level].strctProps[strct->dummyi];

    if ((strct->Type == 7) && (ctrl->dummyd > 1e8))
    {
      if (j == 0)
        sprintf (strctbuff, "%s% 7d  ", mainbuff, strct->ID);
      else
        sprintf (strctbuff, "%s% 7d  ", brnchbuff, strct->ID);

      sprintf (lvl_brnchbuff, "%s% 7d  ", brnchbuff, -1);

      if (level != maxlvls)
        ihsc_prog_tree (ctlgs, strct->ID, level, maxlvls, strctbuff, lvl_brnchbuff, file);
      else
        fprintf (file, "%s\n", strctbuff);

      j++;
    }
  }

  if (ihsc->NumMatch == 0)
  {
    for (i = 0; i < (maxlvls-level); i++)
    {
      sprintf (lvl_mainbuff, "% 7d  ", 0);
      strcat  (mainbuff, lvl_mainbuff);
    }
    fprintf (file, "%s\n", mainbuff);
  }
}
*/
