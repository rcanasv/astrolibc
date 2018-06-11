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
  int          i, j, k, l, n, m;
  Options      opt;
  Structure   * strct1;
  Structure   * strct2;
  Structure   * strct3;
  Structure   * sorted;


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
    Catalog_fill_SubIDS             (&opt.catalog[i]);
    Catalog_fill_isolated           (&opt.catalog[i]);
    if (i < (opt.nsnap - 1))
      stf_read_treefrog (&opt.tree[i], &opt.catalog[i]);
  }


  // --------------------------------------------------- //


  //
  // Tag central galaxy and add stellar mass
  //
  for (i = 0; i < opt.nsnap; i++)
  {
    for (j = 1; j <= opt.catalog[i].nstruct; j++)
    {
      strct1 = &opt.catalog[i].strctProps[j];
      strct1->dummyd = 0.0;
    }
  }

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


  // --------------------------------------------------- //


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
        fprintf (f, "%5d ", strct2->Central);
        fprintf (f, "%5d ", strct2->ID);
        fprintf (f, "\n");
      }
    }
    fclose (f);
  }

  int top = 5;
  sorted = (Structure *) malloc ((opt.catalog[0].nstruct+1) * sizeof(Structure));
  memcpy(sorted, &opt.catalog[0].strctProps[0], (opt.catalog[0].nstruct+1) * sizeof(Structure));
  qsort (&sorted[1], opt.catalog[0].nstruct, sizeof(Structure), Structure_dummyd_compare);


  // --------------------------------------------------- //


  //
  // Create `evolutionary tracks'
  //
  FILE * f1;
  FILE * f2;
  FILE * f3;
  FILE * f4;
  char   buffer1 [NAME_LENGTH];
  char   buffer2 [NAME_LENGTH];
  char   buffer3 [NAME_LENGTH];
  char   buffer4 [NAME_LENGTH];

  sprintf (buffer1, "%s.track_mihsc",  opt.output.prefix);
  sprintf (buffer2, "%s.track_mctrl",  opt.output.prefix);
  sprintf (buffer3, "%s.track_mstot",  opt.output.prefix);
  sprintf (buffer4, "%s.track_idctrl", opt.output.prefix);

  f1 = fopen (buffer1, "w");
  f2 = fopen (buffer2, "w");
  f3 = fopen (buffer3, "w");
  f4 = fopen (buffer4, "w");

  Structure * ihsc;
  Structure * ctrl;
  Structure * ihscp;
  Structure * ctrlp;

  //  for (i = 1; i <= opt.catalog[0].nstruct; i++)
  for (i = opt.catalog[0].nstruct, k = 0; k < top; i--, k++)
  {
    ctrl = &opt.catalog[0].strctProps[sorted[i].ID];
    ihsc = &opt.catalog[0].strctProps[ctrl->HostID];

    fprintf (f1, "%e  ", ihsc->TotMass);
    fprintf (f2, "%e  ", ctrl->TotMass);
    fprintf (f3, "%e  ", ctrl->dummyd);
    fprintf (f4, "%5d  %5d   |   ", ctrl->ID, ihsc->dummyi);

    for (j = 1; j < opt.ntrees; j++)
    {
//      ihscp = &opt.catalog[j].strctProps[ihsc->MatchIDs[0]];
//      ctrlp = &opt.catalog[j].strctProps[ihscp->dummyi];
      ctrlp = &opt.catalog[j].strctProps[ctrl->MatchIDs[0]];
      ihscp = &opt.catalog[j].strctProps[ctrlp->HostID];

      fprintf (f1, "%e  ", ihscp->TotMass);
      fprintf (f2, "%e  ", ctrlp->TotMass);
      fprintf (f3, "%e  ", ctrlp->dummyd);
      fprintf (f4, "%5d  %5d   |   ", ctrlp->ID, ihscp->dummyi);

      ihsc = ihscp;
      ctrl = ctrlp;
    }

    fprintf (f1, "\n");
    fprintf (f2, "\n");
    fprintf (f3, "\n");
    fprintf (f4, "\n");
  }
  fclose (f1);
  fclose (f2);
  fclose (f3);
  fclose (f4);


  // --------------------------------------------------- //


  //
  // Write snapshots for visualization
  //
  int          numpart;
  int       ** strct_to_get;
  Particle   * P;
  gheader      header;

  strct_to_get = (int **) malloc (opt.nsnap * (sizeof(int *)));
  for (i = 0; i < opt.nsnap; i++)
  {
    strct_to_get[i] = (int *) malloc (opt.catalog[i].nstruct+1);
    for (j = 1; j <= opt.catalog[i].nstruct; j++)
      strct_to_get[i][j] = 0;
  }

  //
  // Tag galaxies to be extracted.
  // For visualization purposes label is used
  // also as type of particles
  //
  //   1 - central
  //   2 - ihsc
  //   3 - rest
  //
  for (i = opt.catalog[0].nstruct, k = 0; k < top; i--, k++)
  {
    numpart = 0;
    ctrl = &opt.catalog[0].strctProps[sorted[i].ID];
    ihsc = &opt.catalog[0].strctProps[ctrl->HostID];

    strct_to_get[0][ihsc->ID] = 2;
    numpart += ihsc->NumPart;

    for (n = 0; n < ihsc->NumSubs; n++)
    {
      strct1 = &opt.catalog[0].strctProps[ihsc->SubIDs[n]];
      strct_to_get[0][strct1->ID] = 3;
      numpart += strct1->NumPart;
    }

    numpart += ctrl->NumPart;
    strct_to_get[0][ctrl->ID] = numpart;

    for (j = 1; j < opt.ntrees; j++)
    {
//      ihscp = &opt.catalog[j].strctProps[ihsc->MatchIDs[0]];
//      ctrlp = &opt.catalog[j].strctProps[ihscp->dummyi];
      ctrlp = &opt.catalog[j].strctProps[ctrl->MatchIDs[0]];
      ihscp = &opt.catalog[j].strctProps[ctrlp->HostID];

      strct_to_get[j][ihscp->ID] = 2;
      numpart += ihscp->NumPart;

      for (n = 0; n < ihscp->NumSubs; n++)
      {
        strct1 = &opt.catalog[j].strctProps[ihscp->SubIDs[n]];
        strct_to_get[j][strct1->ID] = 3;
        numpart += strct1->NumPart;
      }

      numpart += ctrlp->NumPart;
      strct_to_get[j][ctrlp->ID] = numpart;

      ihsc = ihscp;
      ctrl = ctrlp;
    }
  }


  // Load Particles
  for (i = 0; i < opt.nsnap; i++)
    Structure_get_particle_properties (&opt.catalog[i], &opt.simulation[i], strct_to_get[i]);


  // Write Gadget Snapshots
  for (i = opt.catalog[0].nstruct, k = 0; k < top; i--, k++)
  {
    m = 0;
    numpart = strct_to_get[0][sorted[i].ID];
    P = (Particle *) malloc (numpart * sizeof(Particle));

    ctrl = &opt.catalog[0].strctProps[sorted[i].ID];
    ihsc = &opt.catalog[0].strctProps[ctrl->HostID];

    // Tag IHSC as 2
    for (n = 0; n < ihsc->NumPart; n++, m++)
    {
      Particle_copy (&ihsc->Part[n], &P[m]);
      P[m].Type = 2;
    }

    // Tag Satellites as 3
    for (n = 0; n < ihsc->NumSubs; n++)
    {
      strct1 = &opt.catalog[0].strctProps[ihsc->SubIDs[n]];
      for (l = 0; l < strct1->NumPart; l++, m++)
      {
        Particle_copy (&strct1->Part[l], &P[m]);
        P[m].Type = 3;
      }
    }

    // Tag Central as 1
    for (n = 0; n < ctrl->NumPart; n++, m++)
    {
      Particle_copy (&ctrl->Part[l], &P[m]);
      P[m].Type = 1;
    }

    if (m == numpart)
    {
      sprintf (opt.output.name, "%s_%d.gdt_%03d",opt.output.prefix, k, 0);
      gadget_write_snapshot (P, m, &header, &opt.output);
    }
    else
      printf ("Error in total number of particles\n");
    free (P);


    for (j = 1; j < opt.ntrees; j++)
    {
//      ihscp = &opt.catalog[j].strctProps[ihsc->MatchIDs[0]];
//      ctrlp = &opt.catalog[j].strctProps[ihscp->dummyi];
      ctrlp = &opt.catalog[j].strctProps[ctrl->MatchIDs[0]];
      ihscp = &opt.catalog[j].strctProps[ctrlp->HostID];

      m = 0;
      numpart = strct_to_get[j][ctrlp->ID];
      P = (Particle *) malloc (numpart * sizeof(Particle));

      // Tag IHSC as 2
      for (n = 0; n < ihscp->NumPart; n++, m++)
      {
        Particle_copy (&ihscp->Part[n], &P[m]);
        P[m].Type = 2;
      }

      // Tag Satellites as 3
      for (n = 0; n < ihscp->NumSubs; n++)
      {
        strct1 = &opt.catalog[j].strctProps[ihscp->SubIDs[n]];
        for (l = 0; l < strct1->NumPart; l++, m++)
        {
          Particle_copy (&strct1->Part[l], &P[m]);
          P[m].Type = 3;
        }
      }

      // Tag Central as 1
      for (n = 0; n < ctrlp->NumPart; n++, m++)
      {
        Particle_copy (&ctrlp->Part[l], &P[m]);
        P[m].Type = 1;
      }


      if (m == numpart)
      {
        sprintf (opt.output.name, "%s_%d.gdt_%03d",opt.output.prefix, k, j);
        gadget_write_snapshot (P, m, &header, &opt.output);
      }
      else
        printf ("Error in total number of particles\n");
      free (P);


      ihsc = ihscp;
      ctrl = ctrlp;
    }
  }


  //
  // Free catalogues
  //
  for (i = 0; i < opt.nsnap; i++)
    Catalog_free (&opt.catalog[i]);
  for (i = 0; i < opt.nsnap; i++)
    free (strct_to_get[i]);
  free (strct_to_get);
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
