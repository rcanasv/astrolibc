/*
 *
 *  \file    calcSO.c
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
  int            iExtract;
  int            iSO;
  int            nsnap;
  Archive        param;
  Archive        output;
  Catalog        catalog;
  Simulation     simulation;
} Options;


void  calcSO_usage   (int opt,  char ** argv);
int   calcSO_options (int argc, char ** argv, Options * opt);
void  calcSO_params  (Options * opt);


int main (int argc, char ** argv)
{
  int           i, j, k, l, n, m;
  Options       opt;

  Particle    * P;
  Structure   * strct1;
  Structure   * strct2;
  Structure   * strct3;
  Structure   * sorted;
  Structure     tmpstrct;
  Grid          myGrid;
  gheader       header;

  FILE  * f;
  char    fname  [NAME_LENGTH];
  char    buffer [NAME_LENGTH];
  int     tmpid;
  int     nfiles;
  int     numpart;
  int   * strct_to_get;
  int   * files_to_read;
  int   * files_of_strct = NULL;

  calcSO_options (argc, argv, &opt);
  calcSO_params  (&opt);

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

  Simulation_init                 (&opt.simulation);
  Catalog_init                    (&opt.catalog);
  Catalog_load_properties         (&opt.catalog);
  Catalog_fill_SubIDS             (&opt.catalog);
  Catalog_fill_isolated           (&opt.catalog);

  // --------------------------------------------------- //
  // Tag central galaxy
  for (i = 0; i < opt.nsnap; i++)
  {
    for (j = 1; j <= opt.catalog.nstruct; j++)
    {
      strct1 = &opt.catalog.strctProps[j];
      strct1->dummyd = 0.0;
      strct1->dummyi = 0;
      strct1->dummy  = 0;
    }
  }

  for (j = 1; j <= opt.catalog.nstruct; j++)
  {
    strct1 = &opt.catalog.strctProps[j];
    if (strct1->Central == 1 && strct1->HostID > 0)
    {
      strct2 = &opt.catalog.strctProps[strct1->HostID];
      strct1->dummyd = strct1->TotMass + strct2->TotMass;
      strct2->dummyi = strct1->ID;
    }
  }

  for (j = 1; j <= opt.catalog.nstruct; j++)
  {
    strct1 = &opt.catalog.strctProps[j];
    if (strct1->Central != 1 && strct1->Type > 7)
    {
      strct2 = &opt.catalog.strctProps[strct1->HostID];
      strct3 = &opt.catalog.strctProps[strct2->dummyi];
      strct3->dummyd += strct1->TotMass;
    }
  }

  strct_to_get  = (int *) malloc ((opt.catalog.nstruct+1)*sizeof(int));
  for (i = 1; i <= opt.catalog.nstruct; i++)
    strct_to_get[i] = 0;

  sprintf (fname, "%s/%s.filesofgroup", opt.catalog.archive.path, opt.catalog.archive.name);
  f = fopen (fname, "r");

  int idtoget;

  for (i = 1; i <= opt.catalog.nstruct; i++)
  {
    strct_to_get[i] = 0;

    fgets  (buffer, NAME_LENGTH, f);
    sscanf (buffer, "%d  %d", &tmpid, &nfiles);
    fgets  (buffer, NAME_LENGTH, f);

    get_n_num_from_string (buffer, nfiles, &files_of_strct);
    if ((nfiles == 1) && (files_of_strct[0] == 0))
    {
      strct_to_get[i] = 1;
      idtoget = i;
    }
    free (files_of_strct);
  }
  fclose (f);


  double galpos [3];
  strct1 = &opt.catalog.strctProps[idtoget];
  strct2 = &opt.catalog.strctProps[strct1->dummyi];
  galpos[0] = strct2->Pos[0];
  galpos[1] = strct2->Pos[1];
  galpos[2] = strct2->Pos[2];

  ramses_amr_init                 (&myGrid);
  ramses_amr_load                 (&opt.simulation, 0, &myGrid);
  ramses_hydro_read               (&opt.simulation, 0, &myGrid);

  int ngaspart  = 0;
  int ndmpart   = 0;
  int nstarpart = 0;

  for (k = myGrid.nlevelmax-1; k >= 9; k--)
  {
    for (i = 0; i < myGrid.level[k].num; i++)
    {
      for (j = 0; j < 8; j++)
      {
        if (myGrid.level[k].cell[i].okOct[j])
          ngaspart++;
      }
    }
  }

  //
  // Load all particles in file(s)
  //


  //
  // Allocate memory for DM, Stars and Gas
  //
  Particle * gasPart = (Particle *) malloc (ngaspart * sizeof(Particle));
  n = 0;
  double lbox = opt.simulation.Lbox;
  double dx, vol;
  double unit_m = opt.simulation.unit_m;

  for (k = myGrid.nlevelmax-1; k >= 9; k--)
  {
    for (i = 0; i < myGrid.level[k].num; i++)
    {
      for (j = 0; j < 8; j++)
      {
        if (myGrid.level[k].cell[i].okOct[j])
        {
          gasPart[n].Pos[0] = lbox * myGrid.level[k].cell[i].octPos[j][0];
          gasPart[n].Pos[1] = lbox * myGrid.level[k].cell[i].octPos[j][1];
          gasPart[n].Pos[2] = lbox * myGrid.level[k].cell[i].octPos[j][2];
          dx  = myGrid.level[k].cell[i].dx;
          vol = dx * dx * dx;
          gasPart[n].Mass = unit_m * vol * myGrid.level[k].cell[i].octRho[j];
          gasPart[n].Id = n;
          gasPart[n].Type = 1;
          n++;
        }
      }
    }
  }

  for (i = 0; i < ngaspart; i++)
  {
    gasPart[i].Pos[0] -= galpos[0];
    gasPart[i].Pos[1] -= galpos[1];
    gasPart[i].Pos[2] -= galpos[2];
    Particle_get_radius (&gasPart[i]);
  }
  qsort (gasPart, ngaspart, sizeof(Particle), Particle_rad_compare);

  double msum = 0.0;
  double pi  = acos(-1.0);
  double fac = 4.0 * pi / 3.0;
  double G = 43009.1e-10;                // in (kpc/M_sun)*(km/s)^2
  double H = opt.simulation.cosmology.HubbleParam / 1000.0; // in (km/s)/kpc
  double crit    = 3.0 * H * H / (8.0 * pi * G);
  double crit200 = 200.0 * crit;
  double rad, rho;
  int  ninR200;

   printf ("Rho_crit  %e\n", crit);

  for (i = 0; i < ngaspart; i++)
  {
    msum += gasPart[i].Mass;
    rad = gasPart[i].Radius;
    rho = msum / (fac * rad * rad * rad);
    if (rho < crit200)
      break;
    else
      gasPart[i].dummyi = 1;
  }
  printf ("R200  %e  M200  %e  Rho  %e\n", rad, msum, rho);

  for (i = 0, ninR200 = 0; i < ngaspart; i++)
    if (gasPart[i].dummyi == 1)
      ninR200++;

  Particle * partinR200;
  partinR200 = (Particle *) malloc (ninR200 * sizeof(Particle));
  for (i = 0, j = 0; i < ngaspart; i++)
    if (gasPart[i].dummyi == 1)
      Particle_copy (&gasPart[i], &partinR200[j++]);

  //FILE * ff = fopen("gaspart", "w");
  //for (i = 0; i < ngaspart; i++)
  //  fprintf (ff, "%e  %e  %e  %e\n", gasPart[i].Pos[0], gasPart[i].Pos[1], gasPart[i].Pos[2], gasPart[i].Mass);
  //fclose (ff);
  gadget_write_snapshot (partinR200, ninR200, &header, &opt.output);

  free (gasPart);
  free (partinR200);

  ramses_amr_free (&myGrid);
  return 0;


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
  // Free memory
  Catalog_free (&opt.catalog);
  ramses_amr_free (&myGrid);
  return 0;

  // --------------------------------------------------- //
  // Diffuse stellar fraction
  //FILE * f;
  //char buffer [NAME_LENGTH];

  if (opt.iSO)
  {
    int     nsat_m09;
    int     nsat_m10;
    int     nsat_m11;
    double  fmass;
    double  radius;
    double  minsat_m08;
    double  minsat_m09;
    double  minsat_m10;
    double  minsat_m11;

      sprintf (buffer, "%s.ihsc", opt.catalog.archive.prefix);
      f = fopen (buffer, "w");
      for (j = 1; j <= opt.catalog.nstruct; j++)
      {
        strct1 = &opt.catalog.strctProps[j];              // IHSC
        strct2 = &opt.catalog.strctProps[strct1->dummyi]; // Central
        if (strct1->Type == 7 && strct1->NumSubs > 0)
        {
          nsat_m09 = 0;
          nsat_m10 = 0;
          nsat_m11 = 0;
          minsat_m08 = 0.0;
          minsat_m09 = 0.0;
          minsat_m10 = 0.0;
          minsat_m11 = 0.0;

          for (k = 1; k < strct1->NumSubs; k++)
          {
            strct3 = &opt.catalog.strctProps[strct1->SubIDs[k]];

            if (strct3->TotMass >= 1e8)
              minsat_m08 += strct3->TotMass;
            if (strct3->TotMass >= 1e9)
              minsat_m09 += strct3->TotMass;
            if (strct3->TotMass >= 1e10)
              minsat_m10 += strct3->TotMass;
            if (strct3->TotMass >= 1e11)
              minsat_m11 += strct3->TotMass;

            if (strct3->TotMass >= 1e9  && strct3->TotMass < 1e10)
              nsat_m09++;
            if (strct3->TotMass >= 1e10 && strct3->TotMass < 1e11)
              nsat_m10++;
            if (strct3->TotMass >= 1e11)
              nsat_m11++;
          }

          /*
          radius = 2.0 * strct2->Rx;
          Structure_calculate_j_r       (strct2, radius);
          Structure_calculate_sigma_v_r (strct2, radius);
          Structure_calculate_sfr       (strct2);
          */

          strct3 = &opt.catalog.strctProps[strct1->SubIDs[strct1->NumSubs-2]];

          fprintf (f, "%e  ", strct2->dummyd);      // Total Stellar Mass
          fprintf (f, "%e  ", strct1->TotMass);     // Mass IHSC
          fprintf (f, "%e  ", strct2->TotMass);     // Mass Central
          fprintf (f, "%e  ", strct3->TotMass);     // Mass Second most massive Gal
          fprintf (f, "%5d ", strct1->NumSubs);     // NumSubs
          fprintf (f, "%5d ", strct2->Central);     // Is Central central?
          fprintf (f, "%5d ", strct1->ID);          // ID IHSC
          fprintf (f, "%5d ", strct2->ID);          // ID Central
          fprintf (f, "%5d ", strct3->ID);          // ID Second most
          fprintf (f, "%5d ", nsat_m09);            // Num sats 1e9  <= M < 1e10
          fprintf (f, "%5d ", nsat_m10);            // Num sats 1e10 <= M < 1e11
          fprintf (f, "%5d ", nsat_m11);            // Num sats 1e11 <= M
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
    }
  // --------------------------------------------------- //

  /*
  // --------------------------------------------------- //
  // Write snapshots for visualization
  if (opt.iExtract)
  {
    strct_to_get = (int *) malloc ((opt.catalog.nstruct+1)*sizeof(int));
    for (j = 1; j <= opt.catalog.nstruct; j++)
      strct_to_get[j] = 0;

    // Load Particles
    Structure_get_particle_properties (&opt.catalog, &opt.simulation, strct_to_get);

    // Write Gadget Snapshots
    for (i = opt.catalog.nstruct, k = 0; ((k < top)&&(i >=1)); i--)
    {
      printf ("boxsize %e\n", opt.simulation[0].Lbox);
      Structure_correct_periodicity (&tmpstrct, &opt.simulation);
      sprintf (opt.output.name, "%s_%d.gdt_%03d",opt.output.prefix, k, 0);
      gadget_write_snapshot (tmpstrct.Part, m, &header, &opt.output);
    }

    free (strct_to_get);
  }
  // --------------------------------------------------- //
  */

  // Free memory
  Catalog_free (&opt.catalog);

  return (0);
}


// --------------------------------------------------- //
//  Parameters
// --------------------------------------------------- //
void calcSO_params (Options * opt)
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

  // Catalogues
  fscanf (opt->param.file, "%s  %s  %s  %d", namebuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->catalog.archive, namebuff);
  Archive_prefix (&opt->catalog.archive, namebuff);
  Archive_format (&opt->catalog.archive, frmtbuff);
  Archive_path   (&opt->catalog.archive, pathbuff);
  Archive_nfiles (&opt->catalog.archive, nflsbuff);

  // Simulation
  fscanf (opt->param.file, "%s  %s  %s  %d", namebuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->simulation.archive, namebuff);
  Archive_prefix (&opt->simulation.archive, namebuff);
  Archive_format (&opt->simulation.archive, frmtbuff);
  Archive_path   (&opt->simulation.archive, pathbuff);
  Archive_nfiles (&opt->simulation.archive, nflsbuff);

  // Close
  fclose (opt->param.file);
}

// --------------------------------------------------- //
//  Options
// --------------------------------------------------- //
int calcSO_options (int argc, char ** argv, Options * opt)
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
    {"so",        0, NULL, 's'},
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
      	opt->iSO = 1;
      	break;

      case 'x':
      	opt->iExtract = 1;
      	break;

      case 'h':
      	calcSO_usage (0, argv);
        break;

      default:
      	calcSO_usage (1, argv);
    }
  }

  if (flag == 0)
    calcSO_usage (1, argv);
}


// --------------------------------------------------- //
//  Usage
// --------------------------------------------------- //
void calcSO_usage (int opt, char ** argv)
{
  if (opt == 0)
  {
    printf ("                                                                         \n");
    printf ("  calcSO                                                                 \n");
    printf ("                                                                         \n");
    printf ("  Author:            Rodrigo Can\\~as                                    \n");
    printf ("                                                                         \n");
    printf ("  Last edition:      22 - 05 - 2019                                      \n");
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
// First check that grid is properly loaded
for (k = 0; k < myGrid.nlevelmax; k++)
{
  sprintf (fname, "amr_lvl_%d", k);
  f = fopen(fname, "w");
  for (i = 0; i < myGrid.level[k].num; i++)
  {
    fprintf (f, "%e  ", myGrid.level[k].cell[i].Pos[0]);
    fprintf (f, "%e  ", myGrid.level[k].cell[i].Pos[1]);
    fprintf (f, "%e  ", myGrid.level[k].cell[i].Pos[2]);
    fprintf (f, "\n");
  }
  fclose (f);
}

sprintf (fname, "hydro_all");
f = fopen(fname, "w");
for (k = 9; k < myGrid.nlevelmax; k++)
{
  for (i = 0; i < myGrid.level[k].num; i++)
  {
    for (j = 0; j < 8; j++)
    {
      fprintf (f, "%e  ", myGrid.level[k].cell[i].octPos[j][0]);
      fprintf (f, "%e  ", myGrid.level[k].cell[i].octPos[j][1]);
      fprintf (f, "%e  ", myGrid.level[k].cell[i].octPos[j][2]);
      fprintf (f, "%e  ", myGrid.level[k].cell[i].octRho[j]);
      fprintf (f, "%d  ", myGrid.level[k].cell[i].okOct[j]);
      fprintf (f, "%d  ", k+1);
      fprintf (f, "\n");
    }
  }
}
fclose (f);

*/
