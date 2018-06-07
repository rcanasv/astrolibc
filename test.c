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
      if (strct1->Central == 1 && strct1->HostID > 0)
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
  // Diffuse stellar fraction
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


  int top = 5;
  /*
  double top_mass[top];
  int    top_id[top];
  for (i = 0; i < top; i++)
    massive[i] = 0.0;
  double tmp;

  for (i = 1; i <= &opt.catalog[0].nstruct; i++)
  {
    strct1 = &opt.catalog[0].strctProps[i];
    if (strct1->Type == 7)
    {
      j = 0;
      do
      {
        if (strct->TotMass > massive[j])
        {
          for (k = top-2; k > j; k--)
          {
            massive[k+1] = massive[k];
          }
        }
      } while(j < top);
    }
  }
  */


  //
  // Create `evolutionary tracks'
  //

  FILE * f1;
  FILE * f2;
  FILE * f3;
  char   buffer1 [NAME_LENGTH];
  char   buffer2 [NAME_LENGTH];
  char   buffer3 [NAME_LENGTH];

  sprintf (buffer1, "%s.track_mihsc", opt.output.prefix);
  sprintf (buffer2, "%s.track_mctrl", opt.output.prefix);
  sprintf (buffer3, "%s.track_mstot", opt.output.prefix);

  f1 = fopen (buffer1, "w");
  f2 = fopen (buffer2, "w");
  f3 = fopen (buffer3, "w");

  Structure * ihsc;
  Structure * ctrl;
  Structure * ihscp;
  Structure * ctrlp;

//  for (i = 1; i <= opt.catalog[0].nstruct; i++)
  for (i = 1; i <= top; i++)
  {
    ihsc = &opt.catalog[0].strctProps[i];
    if (ihsc->Type == 7)
    {
      ctrl = &opt.catalog[0].strctProps[ihsc->dummyi];

      fprintf (f1, "%e  ", ihsc->TotMass);
      fprintf (f2, "%e  ", ctrl->TotMass);
      fprintf (f3, "%e  ", ctrl->dummyd);

      for (j = 1; j < opt.ntrees; j++)
      {
        ihscp = &opt.catalog[j].strctProps[ihsc->MatchIDs[0]];
        ctrlp = &opt.catalog[j].strctProps[ctrl->MatchIDs[0]];

        fprintf (f1, "%e  ", ihscp->TotMass);
        fprintf (f2, "%e  ", ctrlp->TotMass);
        fprintf (f3, "%e  ", ctrlp->dummyd);

        ihsc = ihscp;
        ctrl = ctrlp;
      }

      fprintf (f1, "\n");
      fprintf (f2, "\n");
      fprintf (f3, "\n");
    }
  }
  fclose (f1);
  fclose (f2);
  fclose (f3);


  //
  // Free catalogues
  //
  for (i = 0; i < opt.nsnap; i++)
    Catalog_free (&opt.catalog[i]);

  free (opt.catalog);
  free (opt.tree);
  free (opt.simulation);

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

  fscanf (opt->param.file, "%d", &opt->nsnap);
  opt->ntrees = opt->nsnap - 1;

  opt->catalog    = (Catalog    *) malloc (opt->nsnap * sizeof(Catalog));
  opt->tree       = (Archive    *) malloc (opt->nsnap * sizeof(Archive));
  opt->simulation = (Simulation *) malloc (opt->nsnap * sizeof(Simulation));


  // Catalogues
  for (i = 0; i < opt->nsnap; i++)
  {
    fscanf (opt->param.file, "%s  %s  %s  %d", namebuff, frmtbuff, pathbuff, &nflsbuff);
    Archive_name   (&opt->catalog[i].archive, namebuff);
    Archive_prefix (&opt->catalog[i].archive, namebuff);
    Archive_format (&opt->catalog[i].archive, frmtbuff);
    Archive_path   (&opt->catalog[i].archive, pathbuff);
    Archive_nfiles (&opt->catalog[i].archive, nflsbuff);
  }

  // Trees
  for (i = 0; i < opt->ntrees; i++)
  {
    fscanf (opt->param.file, "%s  %s  %s  %d", namebuff, frmtbuff, pathbuff, &nflsbuff);
    Archive_name   (&opt->tree[i], namebuff);
    Archive_prefix (&opt->tree[i], namebuff);
    Archive_format (&opt->tree[i], frmtbuff);
    Archive_path   (&opt->tree[i], pathbuff);
    Archive_nfiles (&opt->tree[i], nflsbuff);
  }

  // Simulation
  for (i = 0; i < opt->nsnap; i++)
  {
    fscanf (opt->param.file, "%s  %s  %s  %d", namebuff, frmtbuff, pathbuff, &nflsbuff);
    Archive_name   (&opt->simulation[i].archive, namebuff);
    Archive_prefix (&opt->simulation[i].archive, namebuff);
    Archive_format (&opt->simulation[i].archive, frmtbuff);
    Archive_path   (&opt->simulation[i].archive, pathbuff);
    Archive_nfiles (&opt->simulation[i].archive, nflsbuff);
  }

  // Output
  fscanf (opt->param.file, "%s  %s  %s  %d", namebuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->output, namebuff);
  Archive_prefix (&opt->output, namebuff);
  Archive_format (&opt->output, frmtbuff);
  Archive_path   (&opt->output, pathbuff);
  Archive_nfiles (&opt->output, nflsbuff);

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
