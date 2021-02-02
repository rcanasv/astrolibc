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
  Archive        param;
  Archive        output;
  Catalog        catalog;
  Simulation     simulation;
} Options;


void  test_usage   (int opt,  char ** argv);
int   test_options (int argc, char ** argv, Options * opt);
void  test_params  (Options * opt);


int main (int argc, char ** argv)
{
  int          i, j, k;
  Options      opt;


  test_options (argc, argv, &opt);
  test_params  (&opt);


  //
  //  Load catalogs
  //
  //printf ("simulation init\n");
  Simulation_init (&opt.simulation);

  //printf ("Catalog init\n");
  Catalog_init (&opt.catalog);

  //printf ("Catalog load\n");
  Catalog_load_properties (&opt.catalog);

  //printf ("Loading particle properties\n");
  Catalog_get_particle_properties (&opt.catalog, &opt.simulation);


  Structure * strct;
  Structure * strct1;
  Structure * strct2;
  Particle  * p;

  printf ("Tagging isolated galaxies\n");

  //
  //  Tag isolated galaxies in VELOCIraptor
  //
  int * isolated   = (int *) malloc ((opt.catalog.nstruct+1)*sizeof(int));
  int * looselyint = (int *) malloc ((opt.catalog.nstruct+1)*sizeof(int));
  int * highlyint  = (int *) malloc ((opt.catalog.nstruct+1)*sizeof(int));
  int * central    = (int *) malloc ((opt.catalog.nstruct+1)*sizeof(int));

  for (i = 1; i <= opt.catalog.nstruct; i++)
  {
    strct1 = &opt.catalog.strctProps[i];
    strct1->dummyd = 0;
    strct1->dummyi = 0;
  }

  for (i = 1; i <= opt.catalog.nstruct; i++)
  {
    strct1 = &opt.catalog.strctProps[i];
    strct2 = &opt.catalog.strctProps[strct1->HostID];

    isolated   [i] = 0;
    looselyint [i] = 0;
    highlyint  [i] = 0;
    central    [i] = 0;

    if (strct1->Type > 7)
    {
      if (strct1->HostID == -1)
      {
        strct1->dummyd = strct1->TotMass;
        strct1->dummyi = strct1->ID;
        central[i] = 1;
        if (strct1->NumSubs == 0)
          isolated[i] = 1;
        else
          highlyint[i] = 1;
        continue;
      }

      if (strct2->dummyi == 0)
      {
        strct2->dummyd = strct1->TotMass;
        strct2->dummyi = strct1->ID;
      }

      //
      // Determine if galaxy is TRULY isolated
      //
      if (strct2->NumSubs == 1)
         isolated[i] = 1;
      else
      {
        if ((strct1->Type == 10) && (strct1->NumSubs == 0))
          looselyint[i] = 1;
        else
          highlyint[i] = 1;
      }

      //
      // Save mass of central galaxy
      //
      if (strct1->TotMass > strct2->dummyd)
      {
        strct2->dummyd = strct1->TotMass;
        strct2->dummyi = strct1->ID;
      }
    }
  }

  //
  // Tag central galaxies
  //
  for (i = 1; i <= opt.catalog.nstruct; i++)
  {
    strct1 = &opt.catalog.strctProps[i];
    strct2 = &opt.catalog.strctProps[strct1->HostID];

    if (strct1->Type > 7)
      if (strct1->ID == strct2->dummyi)
        central[i] = 1;
  }

  //
  // Merge central with IHSM
  //
  Structure   tmpstrct;
  FILE      * fpart;
  int         totpart;
  int         n;
  char        buff [NAME_LENGTH];

  for (i = 1; i <= opt.catalog.nstruct; i++)
  {
    strct1 = &opt.catalog.strctProps[i];
    if ((central[i] == 1) && (strct1->HostID > 0) && (strct1->Type > 7))
    {
      strct1 = &opt.catalog.strctProps[i];
      strct2 = &opt.catalog.strctProps[strct1->HostID];

      /*
      if (strct1->TotMass > 11)
      {
        sprintf (buff, "galaxy_%07d", i);
        fpart = fopen(buff, "w");
        for (j = 0; j < strct1->NumPart; j++)
        {
          fprintf (fpart, "%e  ", strct1->Part[j].Pos[0]);
          fprintf (fpart, "%e  ", strct1->Part[j].Pos[1]);
          fprintf (fpart, "%e\n", strct1->Part[j].Pos[2]);
        }
        fclose (fpart);

        sprintf (buff, "diffuse_%07d", i);
        fpart = fopen(buff, "w");
        for (j = 0; j < strct2->NumPart; j++)
        {
          fprintf (fpart, "%e  ", strct2->Part[j].Pos[0]);
          fprintf (fpart, "%e  ", strct2->Part[j].Pos[1]);
          fprintf (fpart, "%e\n", strct2->Part[j].Pos[2]);
        }
        fclose (fpart);
      }
      */

      totpart = strct1->NumPart + strct2->NumPart;

      tmpstrct.Part = (Particle *) malloc (totpart * (sizeof(Particle)));
      for (j = 0, n = 0; j < strct1->NumPart; j++, n++)
        Particle_copy (&strct1->Part[j], &tmpstrct.Part[n]);
      for (j = 0; j < strct2->NumPart; j++, n++)
        Particle_copy (&strct2->Part[j], &tmpstrct.Part[n]);
      free (strct1->Part);

      strct1->NumPart = totpart;
      strct1->Part = (Particle *) malloc (totpart * (sizeof(Particle)));

      for (j = 0; j < totpart; j++)
        Particle_copy (&tmpstrct.Part[j], &strct1->Part[j]);

      /*
      if (strct1->TotMass > 11)
      {
        sprintf (buff, "both_%07d", i);
        fpart = fopen(buff, "w");
        for (j = 0; j < strct1->NumPart; j++)
        {
          fprintf (fpart, "%e  ", strct1->Part[j].Pos[0]);
          fprintf (fpart, "%e  ", strct1->Part[j].Pos[1]);
          fprintf (fpart, "%e\n", strct1->Part[j].Pos[2]);
        }
        fclose (fpart);
      }
      */
      free (tmpstrct.Part);
    }
  }


  double  cmx;
  double  cmy;
  double  cmz;

  double  totmass;
  double  mass100kpc3d;
  double  mass100kpc2d;

  double  halftotmass;
  double  halfmass100kpc3d;
  double  halfmass100kpc2d;

  double  r;
  double  r100kpc3d;
  double  r100kpc2d;

  double r2d;

  FILE * f;
  char buffer [NAME_LENGTH];
  sprintf (buffer, "sizemass_eagle_velociraptor.dat");
  f = fopen (buffer, "w");

  for (i = 1; i < opt.catalog.nstruct; i++)
  {
    strct = &opt.catalog.strctProps[i];

    totmass      = 0.0;
    mass100kpc3d = 0.0;
    mass100kpc2d = 0.0;

    halftotmass      = 0.0;
    halfmass100kpc3d = 0.0;
    halfmass100kpc2d = 0.0;

    r            = 0.0;
    r100kpc3d    = 0.0;
    r100kpc2d    = 0.0;

    if (strct->Type > 7)
    {
      strct->Pos[0] *= 1000;
      strct->Pos[1] *= 1000;
      strct->Pos[2] *= 1000;

      Structure_correct_periodicity       (strct, &opt.simulation);
      Structure_shift_to_centre_of_mass   (strct);
      Structure_get_particle_radius       (strct);

      //
      //  Calculate 3D Sizes
      //
      qsort (strct->Part, strct->NumPart, sizeof(Particle), Particle_rad_compare);

      for (k = 0; k < strct->NumPart; k++)
      {
        // Total Mass
        totmass += strct->Part[k].Mass;

        // Mass Inside 3D 100 kpc
        if (strct->Part[k].Radius < 100.0)
          mass100kpc3d += strct->Part[k].Mass;
      }

      halftotmass      = totmass      / 2.0;
      halfmass100kpc3d = mass100kpc3d / 2.0;
      halfmass100kpc2d = mass100kpc3d / 2.0;
      strct->dummyd    = 0;

      for (k = 0; k < strct->NumPart; k++)
      {
        strct->dummyd += strct->Part[k].Mass;

        if (strct->dummyd < halftotmass)      r         = strct->Part[k].Radius;
        if (strct->dummyd < halfmass100kpc3d) r100kpc3d = strct->Part[k].Radius;
      }

      //
      //  Calculate 2D Sizes
      //
      for (k = 0; k < strct->NumPart; k++)
        strct->Part[k].Pos[2] = 0.0;

      Structure_get_particle_radius (strct);
      qsort (strct->Part, strct->NumPart, sizeof(Particle), Particle_rad_compare);
      /*
      for (k = 0; k < strct->NumPart; k++)
      {
        if (strct->Part[k].Radius < 100.0)
          mass100kpc2d += strct->Part[k].Mass;
      }

      halfmass100kpc2d  = mass100kpc2d / 2.0;
      */
      strct->dummyd     = 0;

      for (k = 0; k < strct->NumPart; k++)
      {
        strct->dummyd += strct->Part[k].Mass;
        if (strct->dummyd < halfmass100kpc2d)
         r100kpc2d = strct->Part[k].Radius;
      }
      fprintf (f, "%d  ", strct->ID);
      fprintf (f, "%e  ", totmass);
      fprintf (f, "%e  ", r);
      fprintf (f, "%e  ", mass100kpc3d);
      fprintf (f, "%e  ", r100kpc3d);
      //fprintf (f, "%e  ", mass100kpc2d);
      fprintf (f, "%e  ", r100kpc2d);
      fprintf (f, "\n");
    }
  }
  fclose (f);

  free (isolated);
  free (looselyint);
  free (highlyint);
  free (central);

  //Catalog_free (&opt.catalog);

  return (0);
}


//
//  Parameters
//
void test_params (Options * opt)
{
  int   i;
  int   dummy;
  char  buffer [NAME_LENGTH];

  opt->param.file = fopen (opt->param.name, "r");
  if (opt->param.file == NULL)
  {
    printf ("Couldn't open file  %s\n", opt->param.name);
    printf ("Exiting...\n");
    exit (0);
  }

  // Catalog
  fscanf (opt->param.file, "%s", buffer);  Archive_name   (&opt->catalog.archive, buffer);
                                           Archive_prefix (&opt->catalog.archive, buffer);
  fscanf (opt->param.file, "%s", buffer);  Archive_format (&opt->catalog.archive, buffer);
  fscanf (opt->param.file, "%s", buffer);  Archive_path   (&opt->catalog.archive, buffer);
  fscanf (opt->param.file, "%d", &dummy);  Archive_nfiles (&opt->catalog.archive, dummy);

  // Simulation
  fscanf (opt->param.file, "%s", buffer);  Archive_name   (&opt->simulation.archive, buffer);
                                           Archive_prefix (&opt->simulation.archive, buffer);
  fscanf (opt->param.file, "%s", buffer);  Archive_format (&opt->simulation.archive, buffer);
  fscanf (opt->param.file, "%s", buffer);  Archive_path   (&opt->simulation.archive, buffer);
  fscanf (opt->param.file, "%d", &dummy);  Archive_nfiles (&opt->simulation.archive, dummy);

  // Output
  fscanf (opt->param.file, "%s", buffer);  Archive_name   (&opt->output, buffer);
                                           Archive_prefix (&opt->output, buffer);
  fscanf (opt->param.file, "%s", buffer);  Archive_format (&opt->output, buffer);
  fscanf (opt->param.file, "%s", buffer);  Archive_path   (&opt->output, buffer);
  fscanf (opt->param.file, "%d", &dummy);  Archive_nfiles (&opt->output, dummy);

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
    {"input",   0, NULL, 'i'},
    {0,         0, NULL, 0}
  };

  while ((myopt = getopt_long (argc, argv, "i:vh", lopts, &index)) != -1)
  {
    switch (myopt)
    {
      case 'i':
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
    printf ("  test.c                                                                 \n");
    printf ("                                                                         \n");
    printf ("  Author:            Rodrigo Can\\~as                                    \n");
    printf ("                                                                         \n");
    printf ("  Last edition:      16 - 03 - 2018                                      \n");
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
