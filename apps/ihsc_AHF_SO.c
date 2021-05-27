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
  int         iVerbose;
  int         iFraction;
  int         iExtract;
  int         iIncludeFOF;
  int         iGetSO;
  int         iGetICL;
  int         iGals;
  int         iAsciiStars;
  int         iProps;
  int         region;
  int         rho;
  Archive     param;
  Archive     output;
  Archive     clean;
  Archive     tree;
  Catalog     stf;
  Catalog     ahf;
  Simulation  sim;
} Options;


void  ihsc_usage   (int opt,  char ** argv);
int   ihsc_options (int argc, char ** argv, Options * opt);
void  ihsc_params  (Options * opt);


int main (int argc, char ** argv)
{
  // 1. Read AHF clean halos ASCII
  // 2. Load Catalogues
  //      - STF
  //      - AHF
  // 3. Load AHF Particle list inside and sort IDs
  // 4. Load Simulation Particles and Extended Output
  // 5. Sort Particles by ID
  // 6. Loop over clean structures and copy particles
  //      - Write Gadget snapshot (if desired)
  // 7. Write IHSC mass fractions

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


  // --------------------------------------------------- //
  //            IHSC SO FROM AHF CLEAN HALOS             //
  // --------------------------------------------------- //
  char         fname    [LONG_LENGTH];
  char         buffer   [LONG_LENGTH];
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
  fgets (buffer, NAME_LENGTH, f);  // Header
  while (fgets(buffer, NAME_LENGTH, f) != NULL) ncleans++;
  rewind (f);
  if ((strct_clean = (Structure *) malloc (ncleans*sizeof(Structure))) == NULL)
  {
    printf ("Cannot allocate memory for Structures\n");
    exit (0);
  }

  fgets  (buffer, NAME_LENGTH, f);  // Header
  for (i = 0; i < ncleans; i++)
  {
    fgets  (buffer, NAME_LENGTH, f);
    sscanf (buffer, "%ld  %ld  %lf  %lf  %lf  %lf  %lf",                \
                        &strct_clean[i].HostID, &strct_clean[i].ID,     \
                        &strct_clean[i].Mvir,   &strct_clean[i].Rvir,   \
                        &strct_clean[i].Pos[0], &strct_clean[i].Pos[1], \
                        &strct_clean[i].Pos[2]);
  }
  fclose (f);


  // 2. Load Catalogues
  //      - STF
  //      - AHF
  Catalog_init            (&opt.stf);
  Catalog_load_properties (&opt.stf);
  Catalog_init            (&opt.ahf);
  Catalog_load_properties (&opt.ahf);
  stf_read_treefrog       (&opt.tree, &opt.stf);


  // 3. Load AHF Particle list inside and sort IDs
  ahf_catalog_get_particle_list (&opt.ahf);
  for (i = 1; i <= opt.ahf.nstruct; i++)
  {
    strct1 = &opt.ahf.strctProps[i];
    qsort (&strct1->PIDs[0], strct1->NumPart, sizeof(long), long_compare);
  }


  // 4. Load Simultion Partilces and Extended Outpu
  int  ninextended;
  int  id, indx;
  stfExtendedOutput * xtndd;

  Simulation_init (&opt.sim);                              // Initialize sim
  numpart = Simulation_get_npart_ThisFile (&opt.sim, 0);   // Get numpart
  Simulation_load_particles (&opt.sim, 0, &P);             // Load particles into P

  for (i = 0; i < numpart; i++)
  {
    P[i].dummyi = 0;
    P[i].StructID = 0;
  }


  // Tag particles host ID
  xtndd = NULL;
  ninextended = 0;
  ninextended = stf_load_extended_output (&opt.stf, 0, &xtndd);

  int starOffset = 0;
  for (k = 0; k < 4; k++)
    starOffset += opt.sim.NpartThisFile[k];

  for (j = 0; j < ninextended; j++)
  {
    id    = xtndd[j].IdStruct;
    indx  = xtndd[j].oIndex;

    if (opt.sim.format == GIZMO_SIMBA)
      indx += starOffset;

    if (id > 0)
      P[indx].StructID = id;
    P[indx].indx = indx;
  }
  free (xtndd);


  // 5. Sort Particles by ID
  qsort(&P[0], numpart, sizeof(Particle), Particle_id_compare);


  // 6. Loop over clean structures and copy particles
  //      - Write Gadget snapshot (if desired)
  gheader  header;
  Archive  output;
  for (i = 0, n = 1; i < ncleans; i++)
  {
    if (strct_clean[i].HostID == opt.region)
    {
      strct1 = &opt.ahf.strctProps[strct_clean[i].ID];
      strct1->PSO = (Particle *) malloc (strct1->NumPart*sizeof(Particle));
      strct1->nSO = strct1->NumPart;

      for (k = 0, j = 0; k < numpart; k++)
      {
        if (P[k].Id == strct1->PIDs[j])
        {
          Particle_copy (&P[k], &strct1->PSO[j]);
          j++;
        }
      }

      for (j = 1; j <= opt.stf.nstruct; j++)
        opt.stf.strctProps[j].dummyi = 0;

      if (opt.iGetSO)
      {
        sprintf (output.name, "%s.ihsc_AHF_%03dc.gdt_%03d", opt.stf.archive.prefix, opt.rho, n);
        gadget_write_snapshot (&strct1->PSO[0], strct1->NumPart, &header, &output);
      }

      strct1->ms200c_str = 0;
      strct1->ms200c_dif = 0;
      register int ndif = 0;


      for (j = 0; j < strct1->nSO; j++)
      {
        //strct1->PSO[j].Pos[0] -= strct1->Pos[0];
        //strct1->PSO[j].Pos[1] -= strct1->Pos[1];
        //strct1->PSO[j].Pos[2] -= strct1->Pos[2];
        strct1->PSO[j].Radius = 1;

        if (strct1->PSO[j].Type == 4)
        {
          //Particle_get_radius(&strct1->PSO[j]);
          //strct1->PSO[j].Radius = 1;
          strct2 = &opt.stf.strctProps[strct1->PSO[j].StructID]; // Cross catalog check

          // For IHSC comp
          //if (strct1->PSO[j].StructID > 0 && strct2->Type > 7 || (strct2->Type == 7 && strct2->NumSubs == 0))
          if (strct2->Type > 7 || (strct2->Type == 7 && strct2->NumSubs == 0)*opt.iIncludeFOF)
          {
            strct1->ms200c_str += strct1->PSO[j].Mass;
            strct2->dummyi = 1;
          }
          else
          {
            strct1->ms200c_dif += strct1->PSO[j].Mass;
            strct1->PSO[j].Radius = 0.;  // Tag these for writing ihsc output
            ndif++;
          }
        }
      }

      // All ICL particles will be at the front of the array
      qsort (&strct1->PSO[0], strct1->nSO, sizeof(Particle), Particle_rad_compare);

      // Write Gadget snapshot of ICL
      if (opt.iGetICL)
      {
        //sprintf (output.name, "%s.ihsc_AHF_%03dc_icl.gdt_%03d", opt.stf.archive.prefix, opt.rho, n-1);
        sprintf (output.name, "%s.ihsc_AHF_%03dc_icl.gdt_%03d", opt.stf.archive.prefix, opt.rho, n-1);
        if (opt.iIncludeFOF)
          sprintf (output.name, "%s.ihsc_AHF_%03dc_3DFOFs_icl.gdt_%03d", opt.stf.archive.prefix, opt.rho, n-1);
        gadget_write_snapshot (&strct1->PSO[0], ndif, &header, &output);
      }

      // Write particle properties of ICL
      if (opt.iAsciiStars)
      {
        sprintf (output.name, "%s.ihsc_AHF_%03dc.stars_icl_%03d", opt.stf.archive.prefix, opt.rho, n-1);
        if (opt.iIncludeFOF)
          sprintf (output.name, "%s.ihsc_AHF_%03dc_3DFOFs.stars_icl_%03d", opt.stf.archive.prefix, opt.rho, n-1);
        FILE * fff = fopen (output.name, "w");
        for (j = 0; j < strct1->nSO; j++)
          if (strct1->PSO[j].Type == 4 || strct1->PSO[j].Type == 5)
          {
            fprintf (fff, "%e  ", strct1->PSO[j].Pos[0]);
            fprintf (fff, "%e  ", strct1->PSO[j].Pos[1]);
            fprintf (fff, "%e  ", strct1->PSO[j].Pos[2]);
            fprintf (fff, "%e  ", strct1->PSO[j].Vel[0]);
            fprintf (fff, "%e  ", strct1->PSO[j].Vel[1]);
            fprintf (fff, "%e  ", strct1->PSO[j].Vel[2]);
            fprintf (fff, "%e  ", strct1->PSO[j].Age);
            fprintf (fff, "%d  ",  (int)strct1->PSO[j].Radius);
            fprintf (fff, "%e  ", strct1->PSO[j].Mass);
            fprintf (fff, "%d  ", strct1->PSO[j].Type);
            fprintf (fff, "%ld ", strct1->PSO[j].Id);
            fprintf (fff, "\n");
          }
        fclose (fff);
      }

      // Now print Galaxies info for GSMF
      if (opt.iGals)
      {
        sprintf (output.name, "%s.ihsc_AHF_%03dc.gals_%03d", opt.stf.archive.prefix, opt.rho, n-1);
        if (opt.iIncludeFOF)
          sprintf (output.name, "%s.ihsc_AHF_%03dc_3DFOFs.gals_%03d", opt.stf.archive.prefix, opt.rho, n-1);
        FILE * ff = fopen (output.name, "w");
        for (j = 1; j <= opt.stf.nstruct; j++)
        {
          strct2 = &opt.stf.strctProps[j];
          if (strct2->dummyi == 1)
          {
            fprintf (ff, "%ld  ", strct2->ID);
            fprintf (ff, "%d  ",  strct2->Type);
            fprintf (ff, "%e  ",  strct2->TotMass);
            fprintf (ff, "%d  ",  strct2->NumPart);
            fprintf (ff, "%d  ",  strct2->NumSubs);
            fprintf (ff, "%e  ",  strct2->Pos[0]);
            fprintf (ff, "%e  ",  strct2->Pos[1]);
            fprintf (ff, "%e  ",  strct2->Pos[2]);
            printf ("%d  %d\n", ctrl->NumMatch, ctrl->iMatch);
            if (strct2->NumMatch)
              fprintf (ff, "%ld  ", strct2->MatchIDs[0]);
            else
              fprintf (ff, "%ld  ", 0);
            fprintf (ff, "\n");
          }
        }
        fclose (ff);
      }
    n++;
    } // Region IF block
  } // End of Ncleans loop


  // 7. Write IHSC mass fractions
  if (opt.iProps)
  {
    sprintf (buffer, "%s.ihsc_AHF_%03dc", opt.stf.archive.prefix, opt.rho);
    if (opt.iIncludeFOF)
      sprintf (buffer, "%s.ihsc_AHF_%03dc_3DFOFs", opt.stf.archive.prefix, opt.rho);
    f = fopen (buffer, "w");
    for (i = 0; i < ncleans; i++)
    {
      if (strct_clean[i].HostID == opt.region)
      {
        strct1 = &opt.ahf.strctProps[strct_clean[i].ID];
        fprintf (f, "%15ld  ", strct1->ID);
        fprintf (f, "%d     ", strct1->NumPart);
        fprintf (f, "%e     ", strct1->Mvir);
        fprintf (f, "%e     ", strct1->Rvir);
        fprintf (f, "%e     ", strct1->Pos[0]);
        fprintf (f, "%e     ", strct1->Pos[1]);
        fprintf (f, "%e     ", strct1->Pos[2]);
        fprintf (f, "%e     ", strct1->Vel[1]);
        fprintf (f, "%e     ", strct1->Vel[2]);
        fprintf (f, "%e     ", strct1->Rmax);
        fprintf (f, "%e     ", strct1->mbpOffset);
        fprintf (f, "%e     ", strct1->comOffset);
        fprintf (f, "%e     ", strct1->Vmax);
        fprintf (f, "%e     ", strct1->Vdisp);
        fprintf (f, "%e     ", strct1->Ekin);
        fprintf (f, "%e     ", strct1->Epot);
        fprintf (f, "%e     ", strct1->cNFW);
        fprintf (f, "%e     ", strct1->ms200c_str + strct1->ms200c_dif);
        fprintf (f, "%e     ", strct1->ms200c_dif);
        fprintf (f, "%e     ", strct1->ms200c_str);
        fprintf (f, "\n");
      }
    }
    fclose (f);
  }
  // --------------------------------------------------- //



  // --------------------------------------------------- //
  //                    Free catalogues                  //
  // --------------------------------------------------- //
  free (P);
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

  // Region
  fscanf (opt->param.file, "%d", &opt->region);

  // Overdensity
  fscanf (opt->param.file, "%d", &opt->rho);

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

  // Clean list
  fscanf (opt->param.file, "%s  %s  %s  %s  %d", prefixbuff, namebuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->clean, namebuff);
  Archive_prefix (&opt->clean, prefixbuff);
  Archive_format (&opt->clean, frmtbuff);
  Archive_path   (&opt->clean, pathbuff);
  Archive_nfiles (&opt->clean, nflsbuff);

  // Simulation
  fscanf (opt->param.file, "%s  %s  %s  %s  %d", prefixbuff, namebuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->sim.archive, namebuff);
  Archive_prefix (&opt->sim.archive, prefixbuff);
  Archive_format (&opt->sim.archive, frmtbuff);
  Archive_path   (&opt->sim.archive, pathbuff);
  Archive_nfiles (&opt->sim.archive, nflsbuff);

  // Tree
  fscanf (opt->param.file, "%s  %s  %s  %s  %d", prefixbuff, namebuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->tree.archive, namebuff);
  Archive_prefix (&opt->tree.archive, prefixbuff);
  Archive_format (&opt->tree.archive, frmtbuff);
  Archive_path   (&opt->tree.archive, pathbuff);
  Archive_nfiles (&opt->tree.archive, nflsbuff);

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
    {"help",         0, NULL, 'h'},
    {"verbose",      0, NULL, 'v'},
    {"param",        0, NULL, 'P'},
    {"extract",      0, NULL, 'x'},
    {"include-fofs", 0, NULL, 'f'},
    {"get-so",       0, NULL, 's'},
    {"get-icl",      0, NULL, 'i'},
    {"gals-info",    0, NULL, 'g'},
    {"properties",   0, NULL, 'p'},
    {"ascii-stars",  0, NULL, 'a'},
    {0,              0, NULL, 0}
  };

  opt->iVerbose    = 0;
  opt->iExtract    = 0;
  opt->iIncludeFOF = 0;
  opt->iGetSO      = 0;
  opt->iGetICL     = 0;
  opt->iGals       = 0;
  opt->iProps      = 0;
  opt->iAsciiStars = 0;

  while ((myopt = getopt_long (argc, argv, "P:hvxfsigpa", lopts, &index)) != -1)
  {
    switch (myopt)
    {
      case 'P':
      	strcpy (opt->param.name, optarg);
        flag++;
        break;

      case 'f':
      	opt->iIncludeFOF = 1;
      	break;

      case 'x':
      	opt->iExtract = 1;
      	break;

      case 'v':
      	opt->iVerbose = 1;
      	break;

      case 's':
      	opt->iGetSO = 1;
      	break;

      case 'i':
      	opt->iGetICL = 1;
      	break;

      case 'g':
      	opt->iGals = 1;
      	break;

      case 'p':
      	opt->iProps = 1;
      	break;

      case 'a':
      	opt->iAsciiStars = 1;
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
