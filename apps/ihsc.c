/*
 *
 *  \file    test.c
 *  \brief
 *
 *
 */


#include "../src/base.h"
#include "../src/typedef.h"
#include "../src/archive.h"
#include "../src/catalog.h"
#include "../src/simulation.h"
#include "../tools/SO.h"


typedef struct Options
{
  int            iVerbose;
  int            iFraction;
  int            iTrack;
  int            iExtract;
  int            iSO;
  int            iSO_AHF;
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


void  ihsc_usage   (int opt,  char ** argv);
int   ihsc_options (int argc, char ** argv, Options * opt);
void  ihsc_params  (Options * opt);
void  ihsc_prog_tree (Catalog * ctlgs, int pid, int plevel, int maxlvls, char * mainbuff, char * brnchbuff, FILE * file, int branch);


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

  ihsc_options (argc, argv, &opt);
  ihsc_params  (&opt);

  printf ("opt.iSO  %d\n", opt.iSO);
  printf ("opt.iTrack %d\n", opt.iTrack);
  printf ("opt.iFraction %d\n", opt.iFraction);
  printf ("opt.iExtract %d\n", opt.iExtract);

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
    if (opt.iExtract || opt.iSO)
      Simulation_init                 (&opt.simulation[i]);

    Catalog_init                    (&opt.catalog[i]);
    Catalog_load_properties         (&opt.catalog[i]);
    Catalog_fill_SubIDS             (&opt.catalog[i]);
    Catalog_fill_isolated           (&opt.catalog[i]);

    if (opt.iTrack && i < (opt.nsnap - 1))
        stf_read_treefrog (&opt.tree[i], &opt.catalog[i]);
  }

  mass_min = pow(10.0, opt.mass_min);
  mass_max = pow(10.0, opt.mass_max);
  frac_min = pow(10.0, opt.frac_min);
  frac_max = pow(10.0, opt.frac_max);
  top      = opt.top;

  if (opt.iVerbose)
  {
    printf ("frac min  %e\n", frac_min);
    printf ("frac max  %e\n", frac_max);
    printf ("mass min  %e\n", mass_min);
    printf ("mass max  %e\n", mass_max);
  }


  // Tag central galaxy and add stellar mass
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


  //sorted = (Structure *) malloc ((opt.catalog[0].nstruct+1) * sizeof(Structure));
  //memcpy (sorted, &opt.catalog[0].strctProps[0], (opt.catalog[0].nstruct+1) * sizeof(Structure));
  //qsort (&sorted[1], opt.catalog[0].nstruct, sizeof(Structure), Structure_dummyd_compare);


  /*
  int * strct2get;
  for (i = 0; i < opt.nsnap; i++)
  {
    strct2get = (int *) malloc ((opt.catalog[i].nstruct+1) * sizeof(int));
    for (j = 1; j <= opt.catalog[i].nstruct; j++)
    {
      strct1 = &opt.catalog[i].strctProps[j];
      if (strct1->Central == 1)
        strct2get[j] = 1;
      else
        strct2get[j] = 0;
    }
    Structure_get_particle_properties (&opt.catalog[i], &opt.simulation[i], strct2get);
    Structure_calculate_fmass_radius (&opt.catalog[i], &opt.simulation[i], strct2get, 0.50);
    if (opt.simulation[i].format == RAMSES || opt.simulation[i].format == RAMSES_STAR)
      ramses_structure_calculate_star_age (&opt.simulation[i], &opt.catalog[i], strct2get);
    free (strct2get);
  }
  */


  // --------------------------------------------------- //
  //              IHSC Mass fraction FOF                 //
  // --------------------------------------------------- //
  FILE * f;
  char buffer [NAME_LENGTH];
  if (opt.iFraction)
  {
    int     nsat_m08;
    int     nsat_m09;
    int     nsat_m10;
    int     nsat_m11;
    double  fmass;
    double  radius;
    double  minsat_m08;
    double  minsat_m09;
    double  minsat_m10;
    double  minsat_m11;

    for (i = 0; i < opt.nsnap; i++)
    {
      sprintf (buffer, "%s.ihsc", opt.catalog[i].archive.prefix);
      //sprintf (buffer, "test.ihsc");
      f = fopen (buffer, "w");
      for (j = 1; j <= opt.catalog[i].nstruct; j++)
      {
        strct1 = &opt.catalog[i].strctProps[j];              // IHSC
        strct2 = &opt.catalog[i].strctProps[strct1->dummyi]; // Central
        if (strct1->Type == 7 && strct1->NumSubs > 0)
        {
          nsat_m08 = 0;
          nsat_m09 = 0;
          nsat_m10 = 0;
          nsat_m11 = 0;
          minsat_m08 = 0.0;
          minsat_m09 = 0.0;
          minsat_m10 = 0.0;
          minsat_m11 = 0.0;

          for (k = 0; k < (strct1->NumSubs-1); k++)
          {
            strct3 = &opt.catalog[i].strctProps[strct1->SubIDs[k]];

            if (strct3->TotMass >= 1e8)
              minsat_m08 += strct3->TotMass;
            if (strct3->TotMass >= 1e9)
              minsat_m09 += strct3->TotMass;
            if (strct3->TotMass >= 1e10)
              minsat_m10 += strct3->TotMass;
            if (strct3->TotMass >= 1e11)
              minsat_m11 += strct3->TotMass;

            if (strct3->TotMass >= 1e8  && strct3->TotMass < 1e9)
              nsat_m08++;
            if (strct3->TotMass >= 1e9  && strct3->TotMass < 1e10)
              nsat_m09++;
            if (strct3->TotMass >= 1e10 && strct3->TotMass < 1e11)
              nsat_m10++;
            if (strct3->TotMass >= 1e11)
              nsat_m11++;

            /*
            fmass = strct3->TotMass/strct2->TotMass;
            if (fmass >= 0.001 && fmass < 0.05)
            {
              nsat_f0p001_0p05++;
              minsat_f0p001_0p05 += strct3->TotMass;
            }
            if (fmass >= 0.05 && fmass < 0.3)
            {
              nsat_f0p05_0p30++;
              minsat_f0p05_0p30 += strct3->TotMass;
            }
            if (fmass >= 0.3)
            {
              nsat_f0p30_1p0++;
              minsat_f0p30_1p0 += strct3->TotMass;
            }
            */
          }

          /*
          radius = 2.0 * strct2->Rx;
          Structure_calculate_j_r       (strct2, radius);
          Structure_calculate_sigma_v_r (strct2, radius);
          Structure_calculate_sfr       (strct2);
          */

          strct3 = &opt.catalog[i].strctProps[strct1->SubIDs[strct1->NumSubs-2]];

          fprintf (f, "%e  ", strct2->dummyd);      // Total Stellar Mass
          fprintf (f, "%e  ", strct1->TotMass);     // Mass IHSC
          fprintf (f, "%e  ", strct2->TotMass);     // Mass Central
          fprintf (f, "%e  ", strct3->TotMass);     // Mass Second most massive Gal
          fprintf (f, "%5d ", strct1->NumSubs);     // NumSubs
          fprintf (f, "%5d ", strct2->Central);     // Is Central central?
          fprintf (f, "%5d ", strct1->ID);          // ID IHSC
          fprintf (f, "%5d ", strct2->ID);          // ID Central
          fprintf (f, "%5d ", strct3->ID);          // ID Second most
          fprintf (f, "%5d ", nsat_m08);            // Num sats 1e8  <= M < 1e9
          fprintf (f, "%5d ", nsat_m09);            // Num sats 1e9  <= M < 1e10
          fprintf (f, "%5d ", nsat_m10);            // Num sats 1e10 <= M < 1e11
          fprintf (f, "%5d ", nsat_m11);            // Num sats 1e11 <= M

          /*
          fprintf (f, "%5d ", nsat_f0p001_0p05);    // Num sats 0.001  <= f < 0.05
          fprintf (f, "%5d ", nsat_f0p05_0p30);     // Num sats 0.05   <= f < 0.3
          fprintf (f, "%5d ", nsat_f0p30_1p0);      // Num sats 0.3    <= f
          fprintf (f, "%e ",  minsat_f0p001_0p05);  // Min sats 0.001  <= f < 0.05
          fprintf (f, "%e ",  minsat_f0p05_0p30);   // Min sats 0.05   <= f < 0.3
          fprintf (f, "%e ",  minsat_f0p30_1p0);    // Min sats 0.3    <= f
          */

          /*
          fprintf (f, "%e ",  strct2->sigma);       // sigma_v(r)
          fprintf (f, "%e ",  strct2->j[3]);        // j(r)
          fprintf (f, "%e ",  radius);              // r
          fprintf (f, "%e ",  strct2->SFR20);       // SFR20
          fprintf (f, "%e ",  strct2->SFR50);       // SFR50
          fprintf (f, "%e ",  strct2->SFR100);      // SFR100
          */
          fprintf (f, "%e ", minsat_m08);            // Mass in sats M > 1e8
          fprintf (f, "%e ", minsat_m09);            // Mass in sats M > 1e9
          fprintf (f, "%e ", minsat_m10);            // Mass in sats M > 1e10
          fprintf (f, "%e ", minsat_m11);            // Mass in sats M > 1e11

          fprintf (f, "%10.5lf ", strct2->Pos[0]);   // Pos
          fprintf (f, "%10.5lf ", strct2->Pos[1]);   // Pos
          fprintf (f, "%10.5lf ", strct2->Pos[2]);   // Pos
          fprintf (f, "\n");
        }
      }
      fclose (f);

      // --------------------------------------------------- //
      //                     Satellites                      //
      // --------------------------------------------------- //
      sprintf (buffer, "%s.sats", opt.catalog[i].archive.prefix);
      f = fopen (buffer, "w");
      for (j = 1; j <= 1; j++)
      {
        strct1 = &opt.catalog[i].strctProps[j];              // IHSC
        strct2 = &opt.catalog[i].strctProps[strct1->dummyi]; // Central
        if (strct1->Type == 7 && strct1->NumSubs > 0)
        {
          for (k = 0; k < (strct1->NumSubs-1); k++)
          {
            strct3 = &opt.catalog[i].strctProps[strct1->SubIDs[k]];
            fprintf (f, "%5d  ", strct3->ID);         // Total Stellar Mass
            fprintf (f, "%e  ", strct3->TotMass);     // Mass IHSC
            fprintf (f, "%e  ", strct3->Pos[0]-strct2->Pos[0]);     // Delta x
            fprintf (f, "%e  ", strct3->Pos[1]-strct2->Pos[1]);     // Delta y
            fprintf (f, "%e  ", strct3->Pos[2]-strct2->Pos[2]);     // Delta z
            fprintf (f, "%e  ", strct3->Vel[0]-strct2->Vel[0]);     // Vx
            fprintf (f, "%e  ", strct3->Vel[1]-strct2->Vel[1]);     // Vy
            fprintf (f, "%e  ", strct3->Vel[2]-strct2->Vel[2]);     // Vz
            fprintf (f, "\n");
          }
        }
      }
      fclose (f);
    }
  }
  // --------------------------------------------------- //



  // --------------------------------------------------- //
  //                    IHSC SO                          //
  // --------------------------------------------------- //
  if (opt.iSO)
  {
    Particle * Pbuff;
    int  msum_ap;
    int  msum_str;
    int  msum_dif;

    double  ms200c_str, ms200b_str, ms500c_str, msbn98_str;
    double  ms200c_dif, ms200b_dif, ms500c_dif, msbn98_dif;

    for (i = 0; i < opt.nsnap; i++)
    {
      // Get SO information
      int * SO_tasks = (int *) malloc (opt.catalog[i].archive.nfiles * sizeof(int));
      for (j = 0; j < opt.catalog[i].archive.nfiles; j++)
        SO_tasks[j] = 1;
      get_structure_SO (&opt.catalog[i], &opt.simulation[i], SO_tasks);

      /*
      // Open file to write
      sprintf (buffer, "%s.ihsc.so", opt.catalog[i].archive.prefix);
      f = fopen (buffer, "w");

      // Loop over structures
      for (j = 1; j <= opt.catalog[i].nstruct; j++)
      {
        strct1 = &opt.catalog[i].strctProps[j];              // IHSC
        strct2 = &opt.catalog[i].strctProps[strct1->dummyi]; // Central
        Pbuff = strct1->PSO;

        for (k = 0; k < strct1->nSO; k++)
        {
          msum_ap  = 0;  // for apperture acum
          msum_str = 0;  // for structure acum
          msum_dif = 0;  // for diffuse acum
          if (Pbuff[k].Type == 4)
          {
            strct3 = &opt.catalog[i].strctProps[Pbuff[k].StructID];

            // Spherical appertues e.g. Pillepich
            if (Pbuff[k].StructID == strct1->ID || Pbuff[k].StructID == strct2->ID)
            {
              msum_ap += Pbuff[k].Mass;
              if (Pbuff[k].Radius < 30.0)    strct1->M30  = msum_ap;
              if (Pbuff[k].Radius < 100.0)   strct1->M100 = msum_ap;
            }

            // For IHSC comp
            if (Pbuff[k].StructID > 0 && strct3->Type > 7)
            {
              msum_str += Pbuff[k].Mass;
              if (Pbuff[k].Radius < strct1->R200c)  ms200c_str = msum_str;
              if (Pbuff[k].Radius < strct1->R200b)  ms200b_str = msum_str;
              if (Pbuff[k].Radius < strct1->R500c)  ms500c_str = msum_str;
              if (Pbuff[k].Radius < strct1->Rbn98)  msbn98_str = msum_str;
            }
            else
            {
              msum_dif += Pbuff[k].Mass;
              if (Pbuff[k].Radius < strct1->R200c)  ms200c_dif = msum_dif;
              if (Pbuff[k].Radius < strct1->R200b)  ms200b_dif = msum_dif;
              if (Pbuff[k].Radius < strct1->R500c)  ms500c_dif = msum_dif;
              if (Pbuff[k].Radius < strct1->Rbn98)  msbn98_dif = msum_dif;
            }
          }
        }
      }

      // Write properties
      for (k = 1; k <= opt.catalog[i].nstruct; k++)
      {
        strct1 = &opt.catalog[i].strctProps[k];              // IHSC
        if (strct1->Type == 7 && strct1->NumSubs > 0 && strct1->nSO > 0)
        {
          strct2 = &opt.catalog[i].strctProps[strct1->dummyi]; // Central
          strct3 = &opt.catalog[i].strctProps[strct1->SubIDs[strct1->NumSubs-2]]; // Scnd

          fprintf (f, "%e  ", strct2->dummyd);      // Total Stellar Mass
          fprintf (f, "%e  ", strct1->TotMass);     // Mass IHSC
          fprintf (f, "%e  ", strct2->TotMass);     // Mass Central
          fprintf (f, "%e  ", strct3->TotMass);     // Mass Second most massive Gal
          fprintf (f, "%5d ", strct1->NumSubs);     // NumSubs
          fprintf (f, "%5d ", strct2->Central);     // Is Central central?
          fprintf (f, "%5d ", strct1->ID);          // ID IHSC
          fprintf (f, "%5d ", strct2->ID);          // ID Central
          fprintf (f, "%5d ", strct3->ID);          // ID Second most
          fprintf (f, "%e  ", strct1->M30);         // 30 kpc stellar mass no sats
          fprintf (f, "%e  ", strct1->M100);        // 100 kpc stellar mass no sats
          fprintf (f, "%e  ", strct1->R500c);       // R500C
          fprintf (f, "%e  ", strct1->M500c);
          fprintf (f, "%e  ", ms500c_str);
          fprintf (f, "%e  ", ms500c_dif);
          fprintf (f, "%e  ", strct1->R200c);       // R200C
          fprintf (f, "%e  ", strct1->M200c);
          fprintf (f, "%e  ", ms200c_str);
          fprintf (f, "%e  ", ms200c_dif);
          fprintf (f, "%e  ", strct1->R200b);       // R200B
          fprintf (f, "%e  ", strct1->M200b);
          fprintf (f, "%e  ", ms200b_str);
          fprintf (f, "%e  ", ms200b_dif);
          fprintf (f, "%e  ", strct1->Rbn98);       // RBN98
          fprintf (f, "%e  ", strct1->Mbn98);
          fprintf (f, "%e  ", msbn98_str);
          fprintf (f, "%e  ", msbn98_dif);
          fprintf (f, "\n");
        }
      }
      fclose (f);
      printf("HERE IN SO 7\n");
      */
    }
  }
  // --------------------------------------------------- //

  // --------------------------------------------------- //
  //            IHSC SO FROM AHF CLEAN HALOS             //
  // --------------------------------------------------- //
  if (opt.iSO_AHF)
  {
    // 1. Read AHF clean halos ASCII
    // 2. Load AHF Catalog
    // 3. Load AHF Particle list inside
    //    Sort IDs
    // 4. Load Particles from Simulation
    // 5. Load Extended Output
    // 6. Sort Particles by ID
  }
  // --------------------------------------------------- //



  // --------------------------------------------------- //
  //            Create `evolutionary tracks'             //
  // --------------------------------------------------- //
  if (opt.iTrack)
  {
    FILE * f1;
    FILE * f2;
    FILE * f3;
    FILE * f4;
    FILE * f5;

    char   buffer1  [NAME_LENGTH];
    char   buffer2  [NAME_LENGTH];
    char   buffer3  [NAME_LENGTH];
    char   buffer4  [NAME_LENGTH];
    char   buffer5  [NAME_LENGTH];

    Structure * ihsc;
    Structure * ctrl;
    Structure * ihscp;
    Structure * ctrlp;
    Structure * ihscpctrl;

    int ok = 0;

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
        printf ("%d  %d\n", k, sorted[i].ID);
        fflush (stdout);
        for (j = 1; j < opt.nsnap; j++)
        {

          printf ("%d  %d\n", ctrl->NumMatch, ctrl->iMatch);
          fflush (stdout);
          if (ctrl->NumMatch)
            ctrlp     = &opt.catalog[j].strctProps[ctrl->MatchIDs[0]];
          else
            break;


          if (ctrlp->HostID != ctrlp->ID)
            ihscp     = &opt.catalog[j].strctProps[ctrlp->HostID];
          else
            break;
          ihscpctrl = &opt.catalog[j].strctProps[ihscp->dummyi];

          ihsc = ihscp;
          ctrl = ctrlp;
        }

        ctrl = &opt.catalog[0].strctProps[sorted[i].ID];
        if (j == opt.nsnap)
        {
          ctrl->dummyi = 1;
          k++;
        }
        else
        {
          ctrl->dummyi = 0;
        }

      }
    }

    sprintf (buffer1,  "%s.track_mihsc",       opt.output.prefix);
    sprintf (buffer2,  "%s.track_mctrl",       opt.output.prefix);
    sprintf (buffer3,  "%s.track_mstot",       opt.output.prefix);
    sprintf (buffer4,  "%s.track_idctrl",      opt.output.prefix);
    sprintf (buffer5,  "%s.track_idihsc",      opt.output.prefix);

    f1  = fopen (buffer1,  "w");
    f2  = fopen (buffer2,  "w");
    f3  = fopen (buffer3,  "w");
    f4  = fopen (buffer4,  "w");
    f5  = fopen (buffer5,  "w");

    for (i = opt.catalog[0].nstruct, k = 0; ((k < top)&&(i >=1)); i--)
    {
      ok = 1;
      ctrl = &opt.catalog[0].strctProps[sorted[i].ID];
      ihsc = &opt.catalog[0].strctProps[ctrl->HostID];

      fihsc = ihsc->TotMass/ctrl->dummyd;
      mstot = ctrl->dummyd;

      if (ctrl->dummyi)
      {
        printf ("%d  %d\n", k, sorted[i].ID);
        fflush (stdout);

        fprintf (f1, "%e  ",  ihsc->TotMass);
        fprintf (f2, "%e  ",  ctrl->TotMass);
        fprintf (f3, "%e  ",  ctrl->dummyd);
        fprintf (f4, "%7d  ", ctrl->ID);
        fprintf (f5, "%7d  ", ihsc->ID);


        for (j = 1; j < opt.nsnap; j++)
        {
          ctrlp     = &opt.catalog[j].strctProps[ctrl->MatchIDs[0]];
          ihscp     = &opt.catalog[j].strctProps[ctrlp->HostID];
          ihscpctrl = &opt.catalog[j].strctProps[ihscp->dummyi];

          fprintf (f1, "%e  ", ihscp->TotMass);
          fprintf (f5, "%7d  ", ihscp->ID);
          fprintf (f2, "%e  ", ctrlp->TotMass);
          fprintf (f3, "%e  ", ihscpctrl->dummyd);
          if (ctrlp->ID == ihscpctrl->ID)
            fprintf (f4, "%7d  ", ctrlp->ID);
          else
            fprintf (f4, "%7d  ", ihscpctrl->ID*(-1));

          ihsc = ihscp;
          ctrl = ctrlp;
        }


        fprintf (f1, "\n"); fflush (f1);
        fprintf (f2, "\n"); fflush (f2);
        fprintf (f3, "\n"); fflush (f3);
        fprintf (f4, "\n"); fflush (f4);
        fprintf (f5, "\n"); fflush (f5);
        k++;
      }
    }
    fclose (f1);
    fclose (f2);
    fclose (f3);
    fclose (f4);
    fclose (f5);
  }
  // --------------------------------------------------- //

  /*
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

      if (ctrl->dummyi)
      {
        //ihsc_prog_tree (opt.catalog, ihsc->ID, 0, opt.nsnap, mainbuff, brnchbuff, ftree);
        ihsc_prog_tree (opt.catalog, ctrl->ID, 0, opt.nsnap, mainbuff, brnchbuff, ftree, 0);
        k++;
      }
    }
    fclose (ftree);
  }
  // --------------------------------------------------- //
  */


  // --------------------------------------------------- //
  //          Write snapshots for visualization          //
  // --------------------------------------------------- //
  if (opt.iExtract)
  {
    Structure * ihsc;
    Structure * ctrl;
    Structure * ihscp;
    Structure * ctrlp;
    Structure * ihscpctrl;

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

      if (ctrl->dummyi)
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
            strct_to_get[0][strct1->ID] = 1;
            numpart += strct1->NumPart;
          }
        }

        numpart += ctrl->NumPart;
        strct_to_get[0][ctrl->ID] = numpart;

        for (j = 1; j < opt.nsnap; j++)
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
              strct_to_get[j][strct1->ID] = 1;
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

    printf ("HERE2\n");

    for (i = 1; i <= opt.catalog[0].nstruct; i++)
      printf ("%d\n", strct_to_get[0][i]);
    // Load Particles
    for (i = 0; i < opt.nsnap; i++)
      Structure_get_particle_properties (&opt.catalog[i], &opt.simulation[i], strct_to_get[i]);

      printf ("HERE2\n");
    // Write Gadget Snapshots
    for (i = opt.catalog[0].nstruct, k = 0; ((k < top)&&(i >=1)); i--)
    {
      ctrl      = &opt.catalog[0].strctProps[sorted[i].ID];
      ihsc      = &opt.catalog[0].strctProps[ctrl->HostID];

      fihsc = ihsc->TotMass/ctrl->dummyd;
      mstot = ctrl->dummyd;

      printf ("HERE2\n");
      if (ctrl->dummyi)
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
              tmpstrct.Part[m].Type = 1;
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

        for (j = 1; j < opt.nsnap; j++)
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
                  tmpstrct.Part[m].Type = 1;
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
          tmpstrct.flg_CorrectedPeriodicity = 0;
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
  //                    Free catalogues                  //
  // --------------------------------------------------- //
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

  printf ("branch  %d  ihsc_id  %d  type   %d  plvel  %d\n", branch, ihsc->ID, ihsc->Type, plevel);
  fflush(stdout);

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
void ihsc_params (Options * opt)
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
  fscanf (opt->param.file, "%s  %s  %s  %s  %d", prefixbuff, namebuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->output, namebuff);
  Archive_prefix (&opt->output, prefixbuff);
  Archive_format (&opt->output, frmtbuff);
  Archive_path   (&opt->output, pathbuff);
  Archive_nfiles (&opt->output, nflsbuff);


  // Limits
  //   - Mass
  //   - Fraction
  //   - Number
  fscanf (opt->param.file, "%lf  %lf", &opt->mass_min, &opt->mass_max);
  fscanf (opt->param.file, "%lf  %lf", &opt->frac_min, &opt->frac_max);
  fscanf (opt->param.file, "%d", &opt->top);


  // Update number of trees
  fscanf (opt->param.file, "%d", &opt->nsnap);
  opt->ntrees = opt->nsnap - 1;


  // Allocate memory
  opt->catalog    = (Catalog    *) malloc (opt->nsnap * sizeof(Catalog));
  if (opt->iExtract || opt->iSO)
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
  if (opt->iExtract || opt->iSO)
  {
    for (i = 0; i < opt->nsnap; i++)
    {
      fscanf (opt->param.file, "%s  %s  %s  %s  %d", prefixbuff, namebuff, frmtbuff, pathbuff, &nflsbuff);
      Archive_name   (&opt->simulation[i].archive, namebuff);
      Archive_prefix (&opt->simulation[i].archive, prefixbuff);
      Archive_format (&opt->simulation[i].archive, frmtbuff);
      Archive_path   (&opt->simulation[i].archive, pathbuff);
      Archive_nfiles (&opt->simulation[i].archive, nflsbuff);
    }
  }

  // Trees
  if (opt->iTrack)
  {
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
    {"fraction",  0, NULL, 'f'},
    {"track",     0, NULL, 't'},
    {"extract",   0, NULL, 'x'},
    {"so",        0, NULL, 's'},
    {"soAHF",     0, NULL, 'a'},
    {0,           0, NULL, 0}
  };

  opt->iVerbose  = 0;
  opt->iFraction = 0;
  opt->iTrack    = 0;
  opt->iExtract  = 0;
  opt->iSO       = 0;
  opt->iSO_AHF   = 0;

  while ((myopt = getopt_long (argc, argv, "p:ftxvhs", lopts, &index)) != -1)
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

      case 'v':
      	opt->iVerbose = 1;
      	break;

      case 'x':
      	opt->iExtract = 1;
      	break;

      case 's':
      	opt->iSO = 1;
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
    printf ("  ihsc                                                                   \n");
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
