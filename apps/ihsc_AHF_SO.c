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
  int            region;
  int            rho;
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

printf ("READ CLEAN\n");

  // 2. Load Catalogues
  //      - STF
  //      - AHF
  Catalog_init            (&opt.stf);
  Catalog_load_properties (&opt.stf);
  Catalog_init            (&opt.ahf);
  Catalog_load_properties (&opt.ahf);


for (i = 1; i <= 10; i++)
  printf ("%ld  %d  %e\n", opt.ahf.strctProps[i].ID, opt.ahf.strctProps[i].NumPart,opt.ahf.strctProps[i].Mvir);

for (i = 1; i <= 10; i++)
  printf ("%ld  %d  %e\n", opt.stf.strctProps[i].ID, opt.stf.strctProps[i].NumPart,opt.stf.strctProps[i].TotMass);

  // 3. Load AHF Particle list inside and sort IDs
  ahf_catalog_get_particle_list (&opt.ahf);
  for (i = 1; i <= opt.ahf.nstruct; i++)
  {
    strct1 = &opt.ahf.strctProps[i];
    qsort (&strct1->PIDs[0], strct1->NumPart, sizeof(long), long_compare);
  }
exit(0);
printf ("CATALOGS LOADED\n");

  // 4. Load Simultion Partilces and Extended Outpu
  int  ninextended;
  int  id, indx;
  stfExtendedOutput * xtndd;

  Simulation_init (&opt.sim);                              // Initialize sim

printf ("SImulation initialized\n");

  numpart = Simulation_get_npart_ThisFile (&opt.sim, 0);   // Get numpart

printf ("NumPart %d\n", numpart);

  Simulation_load_particles (&opt.sim, 0, &P);             // Load particles into P
  for (i = 0; i < numpart; i++)
  {
    P[i].dummyi = 0;
    P[i].StructID = 0;
  }

printf ("SIMULATION LOADED\n");

  // Tag particles host ID
  xtndd = NULL;
  ninextended = 0;
  ninextended = stf_load_extended_output (&opt.stf, 0, &xtndd);
  for (j = 0; j < ninextended; j++)
  {
    id    = xtndd[j].IdStruct;
    indx  = xtndd[j].oIndex;

    if (id > 0)
      P[indx].StructID = id;
    P[indx].indx = indx;
  }
  free (xtndd);


  // 5. Sort Particles by ID
  qsort(&P[0], numpart, sizeof(Particle), Particle_id_compare);

printf ("PARTICLES SORTED\n");

  // 6. Loop over clean structures and copy particles
  //      - Write Gadget snapshot (if desired)
  gheader  header;
  Archive  output;
  for (i = 0, n = 0; i < ncleans; i++)
  {
    if (strct_clean[i].HostID == opt.region)
    {
      strct1 = &opt.ahf.strctProps[strct_clean[i].ID];
printf ("region  %d  line %ld strct  %ld  numpart  %d\n", opt.region, strct_clean[i].ID, strct1->ID, strct1->NumPart);  
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

      sprintf (output.name, "%s.ihsc_AHF_%03dc.gdt_%03d", opt.stf.archive.prefix, opt.rho, n++);
      gadget_write_snapshot (&strct1->PSO[0], strct1->NumPart, &header, &output);

      strct1->ms200c_str = 0;
      strct1->ms200c_dif = 0;
      for (j = 0; j < strct1->nSO; j++)
      {
	strct1->PSO[j].Pos[0] -= strct1->Pos[0];
	strct1->PSO[j].Pos[1] -= strct1->Pos[1];
	strct1->PSO[j].Pos[2] -= strct1->Pos[2];

        if (strct1->PSO[j].Type == 4)
        {
          Particle_get_radius(&strct1->PSO[j]);
          strct2 = &opt.stf.strctProps[strct1->PSO[j].StructID]; // Cross catalog check

          // For IHSC comp
          if (strct1->PSO[j].StructID > 0 && strct2->Type > 7)
            strct1->ms200c_str += strct1->PSO[j].Mass;
          else
            strct1->ms200c_dif += strct1->PSO[j].Mass;
        }
      }
    }
  }


  // 7. Write IHSC mass fractions
  sprintf (buffer, "%s.ihsc_AHF_%03dc", opt.stf.archive.prefix, opt.rho);
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
