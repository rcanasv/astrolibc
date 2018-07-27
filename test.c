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

void ihsc_prog_tree (Catalog * ctlgs, int pid, int plevel, int maxlvls, char * mainbuff, char * brnchbuff, FILE * file, int branch);

int main (int argc, char ** argv)
{
  int           i, j, k, l, n, m;
  Options       opt;
  Structure   * strct1;
  Structure   * strct2;
  Structure   * strct3;
  Structure   * sorted;
  double        fihsc;
  double        mstot;
  double        mass_min;
  double        mass_max;
  double        frac_min;
  double        frac_max;
  int           top;

  test_options (argc, argv, &opt);
  test_params  (&opt);

//  opt.nsnap = 3;
//  opt.ntrees = opt.nsnap-1;

  for (i = 0; i < opt.nsnap; i++)
  {
    Simulation_init                 (&opt.simulation[i]);
    Catalog_init                    (&opt.catalog[i]);
    Catalog_load_properties         (&opt.catalog[i]);
    Catalog_fill_SubIDS             (&opt.catalog[i]);
    Catalog_fill_isolated           (&opt.catalog[i]);

    if (opt.iTrack)
      if (i < opt.ntrees)
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
      strct1->dummyi = 0;
      strct1->dummy  = 0;
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
        strct3->dummyd += strct1->TotMass;
      }
    }
  }

  sorted = (Structure *) malloc ((opt.catalog[0].nstruct+1) * sizeof(Structure));
  memcpy (sorted, &opt.catalog[0].strctProps[0], (opt.catalog[0].nstruct+1) * sizeof(Structure));
  qsort (&sorted[1], opt.catalog[0].nstruct, sizeof(Structure), Structure_dummyd_compare);

  int ID = sorted[opt.catalog[0].nstruct].ID;
  printf ("%d\n", ID);

  // --------------------------------------------------- //
  //
  // Write snapshots for visualization
  //
  int       ** strct_to_get;
  Particle   * P;
  gheader      header;
  Structure    tmpstrct;
  Structure  * ctrl;
  Structure  * ctrlp;
  Structure  * sat;
  Structure  * satp;
  Structure  * satd;

  for (i = 1; i <= opt.catalog[0].nstruct; i++)
  {
    strct1 = &opt.catalog[0].strctProps[i];
    for (j = 1; j < opt.nsnap; j++)
    {
      if (strct1->NumMatch)
      {
        strct2 = &opt.catalog[j].strctProps[strct1->MatchIDs[0]];
        strct2->dummy = strct1->ID;
        strct1 = strct2;
      }
      else
        break;
    }
  }

  //
  // Struct to get array
  //
  strct_to_get = (int **) malloc (opt.nsnap * (sizeof(int *)));
  for (i = 0; i < opt.nsnap; i++)
  {
    strct_to_get[i] = (int *) malloc ((opt.catalog[i].nstruct+1)*sizeof(int));
    for (j = 1; j <= opt.catalog[i].nstruct; j++)
      strct_to_get[i][j] = 0;
  }

  //
  // Tag structures to get
  //
  ctrl = &opt.catalog[0].strctProps[ID];
  strct_to_get[0][ctrl->ID] = 1;
  for (i = 1; i < opt.nsnap; i++)
  {
    if (ctrl->NumMatch)
    {
      ctrlp = &opt.catalog[i].strctProps[ctrl->MatchIDs[0]];
      strct_to_get[i][ctrlp->ID] = 1;
      ctrl = ctrlp;
    }
    else
      break;
  }

  strct1 = &opt.catalog[0].strctProps[ID];
  for (i = 0; i < strct1->NumSubs; i++)
  {
    sat = &opt.catalog[0].strctProps[strct1->SubIDs[i]];
    printf("sat id %d\n", sat->ID);
    strct_to_get[0][sat->ID] = 1;
    sat->Pos[0] = strct1->Pos[0];
    sat->Pos[1] = strct1->Pos[1];
    sat->Pos[2] = strct1->Pos[2];
    ctrl = strct1;
    for (j = 1; j < opt.nsnap; j++)
    {
      if (sat->NumMatch)
      {
        ctrlp = &opt.catalog[j].strctProps[ctrl->MatchIDs[0]];
        satp = &opt.catalog[j].strctProps[sat->MatchIDs[0]];
        strct_to_get[j][satp->ID] = 1;
        satp->Pos[0] = ctrlp->Pos[0];
        satp->Pos[1] = ctrlp->Pos[1];
        satp->Pos[2] = ctrlp->Pos[2];
        sat = satp;
        ctrl = ctrlp;
      }
      else
        break;
    }
  }

  //
  // Load Particles
  //
  for (i = 0; i < opt.nsnap; i++)
    Structure_get_particle_properties (&opt.catalog[i], &opt.simulation[i], strct_to_get[i]);

  //
  //  Write gadget_files
  //
  ctrl = &opt.catalog[0].strctProps[ID];
  Structure_correct_periodicity (ctrl, &opt.simulation[0]);
  sprintf (opt.output.name, "%s_central.gdt_%03d",opt.output.prefix, 0);
  gadget_write_snapshot (ctrl->Part, ctrl->NumPart, &header, &opt.output);
  for (i = 1; i < opt.nsnap; i++)
  {
    ctrlp = &opt.catalog[i].strctProps[ctrl->MatchIDs[0]];
    Structure_correct_periodicity (ctrlp, &opt.simulation[i]);
    sprintf (opt.output.name, "%s_central.gdt_%03d",opt.output.prefix, i);
    gadget_write_snapshot (ctrlp->Part, ctrlp->NumPart, &header, &opt.output);
    ctrl = ctrlp;
  }

printf ("HERE\n");
  ctrl = &opt.catalog[0].strctProps[ID];
  for (i = 0; i < ctrl->NumSubs; i++)
  {
    sat = &opt.catalog[0].strctProps[ctrl->SubIDs[i]];
    printf ("sat  %03d ID  %06d  NumPart %09d\n", i, sat->ID, sat->NumPart);
    Structure_correct_periodicity (sat, &opt.simulation[0]);
    sprintf (opt.output.name, "%s_sat_%02d.gdt_%03d",opt.output.prefix, i, 0);
    gadget_write_snapshot (sat->Part, sat->NumPart, &header, &opt.output);
    for (j = 1; j < opt.nsnap; j++)
    {
      if (sat->NumMatch)
      {
        satp = &opt.catalog[j].strctProps[sat->MatchIDs[0]];
        Structure_correct_periodicity (satp, &opt.simulation[j]);
        sprintf (opt.output.name, "%s_sat_%02d.gdt_%03d",opt.output.prefix, i, j);
        printf ("satp %03d ID  %06d  NumPart %09d\n", i, satp->ID, satp->NumPart);
        gadget_write_snapshot (satp->Part, satp->NumPart, &header, &opt.output);
        sat = satp;
      }
      else
        break;
    }
  }


  //
  // Free catalogues
  //
  for (i = 0; i < opt.nsnap; i++)
    free (strct_to_get[i]);
  free (strct_to_get);
  for (i = 0; i < opt.nsnap; i++)
    Catalog_free (&opt.catalog[i]);
  free (opt.catalog);
  free (opt.simulation);
  if (opt.iTrack)
    free (opt.tree);

  return (0);
}

//
//  IHSC Progenitor tree
//
void ihsc_prog_tree (Catalog * ctlgs, int pid, int plevel, int maxlvls, char * mainbuff, char * brnchbuff, FILE * file, int branch)
{
  int          i, j;
  int          id;
  int          type;
  int          num;
  int          level;
  char         lvl_mainbuff   [3000];
  char         lvl_brnchbuff  [3000];
  char         strctbuff      [3000];
  Structure  * ctrl;
  Structure  * ihsc;
  Structure  * strct;
  Structure  * ctrlp;
  Structure  * ctrlpihsc;
  Structure  * ihscp;
  Structure  * ihscpctrl;

  level = plevel + 1;

  if (branch == 0)
  {
    ctrl = &ctlgs[plevel].strctProps[pid];
    ihsc = &ctlgs[plevel].strctProps[ctrl->HostID];
  }
  else
  {
    ihsc = &ctlgs[plevel].strctProps[pid];
  }

  if (plevel == 0)
  {
    sprintf (mainbuff,   "% 7d  % 7d  ", ihsc->ID, ihsc->ID);
    sprintf (brnchbuff,  "% 7d  % 7d  ", ihsc->ID, -1);
  }

  sprintf (lvl_mainbuff,  "%s", mainbuff);
  sprintf (lvl_brnchbuff, "%s", brnchbuff);

  for (i = 0, j = 0; (i < ihsc->NumMatch && j < 4); i++)
  {
    if (branch == 0)
    {
      ctrlp     = &ctlgs[level].strctProps[ctrl->MatchIDs[0]];
      ctrlpihsc = &ctlgs[level].strctProps[ctrlp->HostID];
      ihscp     = &ctlgs[level].strctProps[ctrlp->HostID];
      ihscpctrl = &ctlgs[level].strctProps[ihscp->dummyi];
      strct     = &ctlgs[level].strctProps[ctrl->MatchIDs[0]];
      id        = ctrlp->ID;
      type      = 10;
      i--;
    }
    else
    {
      ihscp     = &ctlgs[level].strctProps[ihsc->MatchIDs[i]];
      ihscpctrl = &ctlgs[level].strctProps[ihscp->dummyi];
      strct     = &ctlgs[level].strctProps[ihsc->MatchIDs[i]];
      id        = ihscp->ID;
      type      = 7;
      if (ihscp->ID == ctrlpihsc->ID)
        continue;
    }

//    if (!(strct->Type%type) && (ihscpctrl->dummyd > 1e8))
    if (!(strct->Type%type) && (ihscpctrl->dummyd > 1e10))
    {
      if (j == 0)
        sprintf (strctbuff, "%s% 7d  ", mainbuff,  ihscp->ID);
      else
        sprintf (strctbuff, "%s% 7d  ", brnchbuff, ihscp->ID);

      sprintf (lvl_brnchbuff, "%s% 7d  ", brnchbuff, -1);

      if (level != maxlvls)
        ihsc_prog_tree (ctlgs, id, level, maxlvls, strctbuff, lvl_brnchbuff, file, branch);
      else
        fprintf (file, "%s%d\n", strctbuff,branch);

      j++;
      branch++;
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
    branch++;
  }
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

  // Output
  fscanf (opt->param.file, "%s  %s  %s  %s  %d", prfxbuff,namebuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->output, namebuff);
  Archive_prefix (&opt->output, prfxbuff);
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
    fscanf (opt->param.file, "%s  %s  %s  %s  %d", prfxbuff,namebuff, frmtbuff, pathbuff, &nflsbuff);
    Archive_name   (&opt->catalog[i].archive, namebuff);
    Archive_prefix (&opt->catalog[i].archive, prfxbuff);
    Archive_format (&opt->catalog[i].archive, frmtbuff);
    Archive_path   (&opt->catalog[i].archive, pathbuff);
    Archive_nfiles (&opt->catalog[i].archive, nflsbuff);
  }


  // Simulation
  for (i = 0; i < opt->nsnap; i++)
  {
    fscanf (opt->param.file, "%s  %s  %s  %s  %d", prfxbuff,namebuff, frmtbuff, pathbuff, &nflsbuff);
    Archive_name   (&opt->simulation[i].archive, namebuff);
    Archive_prefix (&opt->simulation[i].archive, prfxbuff);
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
      fscanf (opt->param.file, "%s  %s  %s  %s  %d", prfxbuff, namebuff, frmtbuff, pathbuff, &nflsbuff);
      Archive_name   (&opt->tree[i], namebuff);
      Archive_prefix (&opt->tree[i], prfxbuff);
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
