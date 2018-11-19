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

//  opt.nsnap = 3;
//  opt.ntrees = opt.nsnap-1;

  for (i = 0; i < opt.nsnap; i++)
  {
    //Simulation_init                 (&opt.simulation[i]);
    Catalog_init                    (&opt.catalog[i]);
    Catalog_load_properties         (&opt.catalog[i]);
    //Catalog_get_particle_properties (&opt.catalog[i], &opt.simulation[i]);
    Catalog_fill_SubIDS             (&opt.catalog[i]);
    Catalog_fill_isolated           (&opt.catalog[i]);

    if (opt.iTrack)
      if (i < (opt.nsnap - 1))
        stf_read_treefrog (&opt.tree[i], &opt.catalog[i]);
  }

  mass_min = pow(10.0, opt.mass_min);
  mass_max = pow(10.0, opt.mass_max);
  frac_min = pow(10.0, opt.frac_min);
  frac_max = pow(10.0, opt.frac_max);
  top      = opt.top;

  printf ("frac min  %e\n", frac_min);
  printf ("frac max  %e\n", frac_max);
  printf ("mass min  %e\n", mass_min);
  printf ("mass max  %e\n", mass_max);


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

  // --------------------------------------------------- //


  // --------------------------------------------------- //
  //
  // Diffuse stellar fraction
  //
  FILE * f;
  char buffer [NAME_LENGTH];

  if (opt.iFraction)
  {
    int  nsat_m09;
    int  nsat_m10;
    int  nsat_m11;

    for (i = 0; i < opt.nsnap; i++)
    {
      sprintf (buffer, "%s.ihsc", opt.catalog[i].archive.prefix);
      f = fopen (buffer, "w");
      for (j = 1; j <= opt.catalog[i].nstruct; j++)
      {
        strct1 = &opt.catalog[i].strctProps[j];
        if (strct1->Type == 7 && strct1->NumSubs > 0)
        {
          nsat_m09 = 0;
          nsat_m10 = 0;
          nsat_m11 = 0;

          for (k = 1; k < strct1->NumSubs; k++)
          {
            strct3 = &opt.catalog[i].strctProps[strct1->SubIDs[k]];
            if (strct3->TotMass >= 1e9  && strct3->TotMass < 1e10)
              nsat_m09++;
            if (strct3->TotMass >= 1e10 && strct3->TotMass < 1e11)
              nsat_m10++;
            if (strct3->TotMass >= 1e11)
              nsat_m11++;
          }

          strct2 = &opt.catalog[i].strctProps[strct1->dummyi];
          strct3 = &opt.catalog[i].strctProps[strct1->SubIDs[strct1->NumSubs-2]];
          fprintf (f, "%e  ", strct2->dummyd);    // Total Stellar Mass
          fprintf (f, "%e  ", strct1->TotMass);   // Mass IHSC
          fprintf (f, "%e  ", strct2->TotMass);   // Mass Central
          fprintf (f, "%e  ", strct3->TotMass);   // Mass Second most massive Gal
          fprintf (f, "%5d ", strct1->NumSubs);   // NumSubs
          fprintf (f, "%5d ", strct2->Central);   // Is Central central?
          fprintf (f, "%5d ", strct1->ID);        // ID IHSC
          fprintf (f, "%5d ", strct2->ID);        // ID Central
          fprintf (f, "%5d ", strct3->ID);        // ID Second most
          fprintf (f, "%5d ", nsat_m09);          // Num sats 1e9  <= M < 1e10
          fprintf (f, "%5d ", nsat_m10);          // Num sats 1e10 <= M < 1e11
          fprintf (f, "%5d ", nsat_m11);          // Num sats 1e11 <= M
          fprintf (f, "\n");
        }
      }
      fclose (f);
    }
  }
  // --------------------------------------------------- //


  // --------------------------------------------------- //
  //
  // Create `evolutionary tracks'
  //
  FILE * f1;
  FILE * f2;
  FILE * f3;
  FILE * f4;
//  FILE * f5;
//  FILE * f6;
//  FILE * f7;
//  FILE * f8;
//  FILE * f9;
//  FILE * f10;
//  FILE * f11;

  char   buffer1  [NAME_LENGTH];
  char   buffer2  [NAME_LENGTH];
  char   buffer3  [NAME_LENGTH];
  char   buffer4  [NAME_LENGTH];
//  char   buffer5  [NAME_LENGTH];
//  char   buffer6  [NAME_LENGTH];
//  char   buffer7  [NAME_LENGTH];
//  char   buffer8  [NAME_LENGTH];
//  char   buffer9  [NAME_LENGTH];
//  char   buffer10 [NAME_LENGTH];
//  char   buffer11 [NAME_LENGTH];
//  char   buffer12 [NAME_LENGTH];
//  char   buffer13 [NAME_LENGTH];
//  char   buffer14 [NAME_LENGTH];

  Structure * ihsc;
  Structure * ctrl;
  Structure * ihscp;
  Structure * ctrlp;
  Structure * ihscpctrl;

  if (opt.iTrack)
  {
    sprintf (buffer1,  "%s.track_mihsc",       opt.output.prefix);
    sprintf (buffer2,  "%s.track_mctrl",       opt.output.prefix);
    sprintf (buffer3,  "%s.track_mstot",       opt.output.prefix);
    sprintf (buffer4,  "%s.track_idctrl",      opt.output.prefix);
//    sprintf (buffer5,  "%s.track_mmsub",       opt.output.prefix);
//    sprintf (buffer6,  "%s.track_subs_m09",    opt.output.prefix);
//    sprintf (buffer7,  "%s.track_subs_m10",    opt.output.prefix);
//    sprintf (buffer8,  "%s.track_subs_m11",    opt.output.prefix);
//    sprintf (buffer9,  "%s.track_gmerger_m09", opt.output.prefix);
//    sprintf (buffer10, "%s.track_gmerger_m10", opt.output.prefix);
//    sprintf (buffer11, "%s.track_gmerger_m11", opt.output.prefix);
//    sprintf (buffer12, "%s.track_smerger_m09", opt.output.prefix);
//    sprintf (buffer13, "%s.track_smerger_m10", opt.output.prefix);
//    sprintf (buffer14, "%s.track_smerger_m11", opt.output.prefix);

    f1  = fopen (buffer1,  "w");
    f2  = fopen (buffer2,  "w");
    f3  = fopen (buffer3,  "w");
    f4  = fopen (buffer4,  "w");
//    f5  = fopen (buffer5,  "w");
//    f6  = fopen (buffer6,  "w");
//    f7  = fopen (buffer7,  "w");
//    f8  = fopen (buffer8,  "w");
//    f9  = fopen (buffer9,  "w");
//    f10 = fopen (buffer10, "w");
//    f11 = fopen (buffer11, "w");
//    f12 = fopen (buffer12, "w");
//    f13 = fopen (buffer13, "w");
//    f14 = fopen (buffer14, "w");

    for (i = opt.catalog[0].nstruct, k = 0; ((k < top)&&(i >=1)); i--)
    {
      ctrl = &opt.catalog[0].strctProps[sorted[i].ID];
      ihsc = &opt.catalog[0].strctProps[ctrl->HostID];

      fihsc = ihsc->TotMass/ctrl->dummyd;
      mstot = ctrl->dummyd;

      if ((fihsc > frac_min) && \
          (fihsc < frac_max) && \
          (mstot > mass_min) && \
          (mstot < mass_max))
      {
        fprintf (f1, "%e  ",  ihsc->TotMass);
        fprintf (f2, "%e  ",  ctrl->TotMass);
        fprintf (f3, "%e  ",  ctrl->dummyd);
        fprintf (f4, "%7d  ", ctrl->ID);

        for (j = 1; j < opt.nsnap; j++)
        {
          ctrlp     = &opt.catalog[j].strctProps[ctrl->MatchIDs[0]];
          ihscp     = &opt.catalog[j].strctProps[ctrlp->HostID];
          ihscpctrl = &opt.catalog[j].strctProps[ihscp->dummyi];

          fprintf (f1, "%e  ", ihscp->TotMass);
          fprintf (f2, "%e  ", ctrlp->TotMass);
          fprintf (f3, "%e  ", ihscpctrl->dummyd);
          if (ctrlp->ID == ihscpctrl->ID)
            fprintf (f4, "%7d  ", ctrlp->ID);
          else
            fprintf (f4, "%7d  ", ihscpctrl->ID*(-1));

          ihsc = ihscp;
          ctrl = ctrlp;
        }

        fprintf (f1, "\n");
        fprintf (f2, "\n");
        fprintf (f3, "\n");
        fprintf (f4, "\n");

        k++;
      }
    }
    fclose (f1);
    fclose (f2);
    fclose (f3);
    fclose (f4);
  }
  printf ("%d\n", k);
  // --------------------------------------------------- //


  // --------------------------------------------------- //
  //
  // Create System merger trees
  //
  char    mainbuff  [3000];
  char    brnchbuff [3000];
  FILE  * ftree;
  if (opt.iTrack)
  {
    sprintf (buffer1, "%s.track_ihsc_tree", opt.output.prefix);
    ftree  = fopen (buffer1,  "w");
    for (i = opt.catalog[0].nstruct, k = 0; ((k < top)&&(i >=1)); i--)
    {
      ctrl = &opt.catalog[0].strctProps[sorted[i].ID];
      ihsc = &opt.catalog[0].strctProps[ctrl->HostID];

      fihsc = ihsc->TotMass/ctrl->dummyd;
      mstot = ctrl->dummyd;

      if ((fihsc > frac_min) && \
          (fihsc < frac_max) && \
          (mstot > mass_min) && \
          (mstot < mass_max))
      {
        //ihsc_prog_tree (opt.catalog, ihsc->ID, 0, opt.nsnap, mainbuff, brnchbuff, ftree);
        ihsc_prog_tree (opt.catalog, ctrl->ID, 0, opt.nsnap, mainbuff, brnchbuff, ftree, 0);
        k++;
      }
    }
    fclose (ftree);
  }
  // --------------------------------------------------- //


  // --------------------------------------------------- //
  //
  // Write snapshots for visualization
  //
  int          numpart;
  int       ** strct_to_get;
  Particle   * P;
  gheader      header;
  Structure    tmpstrct;

  if (opt.iExtract)
  {
    strct_to_get = (int **) malloc (opt.nsnap * (sizeof(int *)));
    for (i = 0; i < opt.nsnap; i++)
    {
      strct_to_get[i] = (int *) malloc ((opt.catalog[i].nstruct+1)*sizeof(int));
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
    for (i = opt.catalog[0].nstruct, k = 0; ((k < top)&&(i >=1)); i--)
    {
      ctrl      = &opt.catalog[0].strctProps[sorted[i].ID];
      ihsc      = &opt.catalog[0].strctProps[ctrl->HostID];

      fihsc = ihsc->TotMass/ctrl->dummyd;
      mstot = ctrl->dummyd;

      if ((fihsc > frac_min) && \
          (fihsc < frac_max) && \
          (mstot > mass_min) && \
          (mstot < mass_max))
      {
        numpart = 0;
        ctrl = &opt.catalog[0].strctProps[sorted[i].ID];
        ihsc = &opt.catalog[0].strctProps[ctrl->HostID];

        strct_to_get[0][ihsc->ID] = 2;
        numpart += ihsc->NumPart;

        for (n = 0; n < ihsc->NumSubs; n++)
        {
          strct1 = &opt.catalog[0].strctProps[ihsc->SubIDs[n]];
          if (strct1->ID != ctrl->ID)
          {
            strct_to_get[0][strct1->ID] = 3;
            numpart += strct1->NumPart;
          }
        }

        numpart += ctrl->NumPart;
        strct_to_get[0][ctrl->ID] = numpart;

        for (j = 1; j < opt.ntrees; j++)
        {
          numpart = 0;

          ctrlp = &opt.catalog[j].strctProps[ctrl->MatchIDs[0]];
          ihscp = &opt.catalog[j].strctProps[ctrlp->HostID];

          strct_to_get[j][ihscp->ID] = 2;
          numpart += ihscp->NumPart;

          for (n = 0; n < ihscp->NumSubs; n++)
          {
            strct1 = &opt.catalog[j].strctProps[ihscp->SubIDs[n]];
            if (strct1->ID != ctrlp->ID)
            {
              strct_to_get[j][strct1->ID] = 3;
              numpart += strct1->NumPart;
            }
          }

          numpart += ctrlp->NumPart;
          strct_to_get[j][ctrlp->ID] = numpart;

          ihsc = ihscp;
          ctrl = ctrlp;
        }
        k++;
      }
    }

    // Load Particles
    for (i = 0; i < opt.nsnap; i++)
      Structure_get_particle_properties (&opt.catalog[i], &opt.simulation[i], strct_to_get[i]);

    // Write Gadget Snapshots
    for (i = opt.catalog[0].nstruct, k = 0; ((k < top)&&(i >=1)); i--)
    {
      ctrl      = &opt.catalog[0].strctProps[sorted[i].ID];
      ihsc      = &opt.catalog[0].strctProps[ctrl->HostID];

      fihsc = ihsc->TotMass/ctrl->dummyd;
      mstot = ctrl->dummyd;

      if ((fihsc > frac_min) && \
          (fihsc < frac_max) && \
          (mstot > mass_min) && \
          (mstot < mass_max))
      {
        m = 0;
        numpart = strct_to_get[0][sorted[i].ID];
        tmpstrct.Part = (Particle *) malloc (numpart * sizeof(Particle));
        tmpstrct.NumPart = numpart;

        ctrl = &opt.catalog[0].strctProps[sorted[i].ID];
        ihsc = &opt.catalog[0].strctProps[ctrl->HostID];

        // Tag IHSC as 2
        for (n = 0; n < ihsc->NumPart; n++, m++)
        {
          Particle_copy (&ihsc->Part[n], &tmpstrct.Part[m]);
          tmpstrct.Part[m].Type = 2;
        }

        // Tag Satellites as 3
        for (n = 0; n < ihsc->NumSubs; n++)
        {
          strct1 = &opt.catalog[0].strctProps[ihsc->SubIDs[n]];
          if (strct1->ID != ctrl->ID)
          {
            for (l = 0; l < strct1->NumPart; l++, m++)
            {
              Particle_copy (&strct1->Part[l], &tmpstrct.Part[m]);
              tmpstrct.Part[m].Type = 3;
            }
          }
        }

        // Tag Central as 1
        for (n = 0; n < ctrl->NumPart; n++, m++)
        {
          Particle_copy (&ctrl->Part[n], &tmpstrct.Part[m]);
          tmpstrct.Part[m].Type = 1;
        }

        tmpstrct.Pos[0] = ctrl->Pos[0];
        tmpstrct.Pos[1] = ctrl->Pos[1];
        tmpstrct.Pos[2] = ctrl->Pos[2];
        printf ("boxsize %e\n", opt.simulation[0].Lbox);
        Structure_correct_periodicity (&tmpstrct, &opt.simulation[0]);

        if (m == numpart)
        {
          sprintf (opt.output.name, "%s_%d.gdt_%03d",opt.output.prefix, k, 0);
          gadget_write_snapshot (tmpstrct.Part, m, &header, &opt.output);
        }
        else
          printf ("Error in total number of particles\n");
        free (tmpstrct.Part);

        for (j = 1; j < opt.ntrees; j++)
        {
          ctrlp = &opt.catalog[j].strctProps[ctrl->MatchIDs[0]];
          ihscp = &opt.catalog[j].strctProps[ctrlp->HostID];

          m = 0;
          numpart = strct_to_get[j][ctrlp->ID];
          tmpstrct.Part = (Particle *) malloc (numpart * sizeof(Particle));
          tmpstrct.NumPart = numpart;

          // Tag IHSC as 2
          for (n = 0; n < ihscp->NumPart; n++, m++)
          {
            Particle_copy (&ihscp->Part[n], &tmpstrct.Part[m]);
            tmpstrct.Part[m].Type = 2;
          }

          // Tag Satellites as 3
          for (n = 0; n < ihscp->NumSubs; n++)
          {
            strct1 = &opt.catalog[j].strctProps[ihscp->SubIDs[n]];
            if (strct1->ID != ctrlp->ID)
            {
              if (strct1->Central == 0)
              {
                for (l = 0; l < strct1->NumPart; l++, m++)
                {
                  Particle_copy (&strct1->Part[l], &tmpstrct.Part[m]);
                  tmpstrct.Part[m].Type = 3;
                }
              }
              else
              {
                for (l = 0; l < strct1->NumPart; l++, m++)
                {
                  Particle_copy (&strct1->Part[l], &tmpstrct.Part[m]);
                  tmpstrct.Part[m].Type = 4;
                }
              }
            }
          }

          // Tag Central as 1
          for (n = 0; n < ctrlp->NumPart; n++, m++)
          {
            Particle_copy (&ctrlp->Part[n], &tmpstrct.Part[m]);
            tmpstrct.Part[m].Type = 1;
          }

          tmpstrct.Pos[0] = ctrlp->Pos[0];
          tmpstrct.Pos[1] = ctrlp->Pos[1];
          tmpstrct.Pos[2] = ctrlp->Pos[2];
          printf ("boxsize %e\n", opt.simulation[j].Lbox);
          Structure_correct_periodicity (&tmpstrct, &opt.simulation[j]);

          if (m == numpart)
          {
            sprintf (opt.output.name, "%s_%d.gdt_%03d",opt.output.prefix, k, j);
            gadget_write_snapshot (tmpstrct.Part, m, &header, &opt.output);
          }
          else
            printf ("Error in total number of particles\n");
          free (tmpstrct.Part);

          ihsc = ihscp;
          ctrl = ctrlp;
        }
       k++;
      }
    }

    for (i = 0; i < opt.nsnap; i++)
      free (strct_to_get[i]);
    free (strct_to_get);
  }
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

  // Limits
  //   - Mass
  fscanf (opt->param.file, "%lf  %lf", &opt->mass_min, &opt->mass_max);
  //   - Fraction
  fscanf (opt->param.file, "%lf  %lf", &opt->frac_min, &opt->frac_max);
  //   - Number
  fscanf (opt->param.file, "%d", &opt->top);


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
    fscanf (opt->param.file, "%s  %s  %s  %d", namebuff, frmtbuff, pathbuff, &nflsbuff);
    Archive_name   (&opt->catalog[i].archive, namebuff);
    Archive_prefix (&opt->catalog[i].archive, namebuff);
    Archive_format (&opt->catalog[i].archive, frmtbuff);
    Archive_path   (&opt->catalog[i].archive, pathbuff);
    Archive_nfiles (&opt->catalog[i].archive, nflsbuff);
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

  // Trees
  if (opt->iTrack)
  {
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
