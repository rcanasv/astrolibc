/*
 *
 *  \file    ctlgmatch.c
 *  \brief   This code cross matches structures (galaxies, dmhs, ...)
 *           in the same snapshot identified by two different codes.
 *
 *           Currently only VELOCIraptor - HaloMaker cross match is
 *           available.
 *
 */


#include "ctlgmatch.h"


int main (int argc, char ** argv)
{

  int i, j, k;
  Options opt;

  ctlgMatch_options (argc, argv, &opt);

printf ("Reading parameters\n");

  ctlgMatch_params  (&opt);

printf ("Parameters successfully read\n");

  //
  //  Load catalogs
  //
  for (i = 0; i < opt.numCatalogs; i++)
  {
printf ("Initializing catalogs\n");
    Catalog_init (&opt.catalog[i]);
printf ("Loading catalogs\n");
    Catalog_load (&opt.catalog[i]);
    Catalog_get_particle_properties (&opt.catalog[i], &opt.simulation[i]);
  }


  //
  //  Tag isolated galaxies in VELOCIraptor
  //
  int * isolated = (int *) malloc ((opt.catalog[0].nstruct+1)*sizeof(int));

  for (i = 1; i <= opt.catalog[0].nstruct; i++)
  {
    isolated[i] = 0;

    if ((opt.catalog[0].strctProps[i].Type == 10) && \
        (opt.catalog[0].strctProps[i].NumSubs == 0))
      isolated[i] = 1;
  }

  //
  //  Read TreeFrog
  //
  stf_read_treefrog (&opt.mtree, &opt.catalog[0]);

  //
  //  Print common properties
  //
  Structure * strct1;
  Structure * strct2;
  FILE * f = fopen (opt.output.name, "w");
  for (i = 1; i <= opt.catalog[0].nstruct; i++)
  {
    strct1 = &opt.catalog[0].strctProps[i];
    if ((strct1->Type > 7) && (strct1->NumMatch))
    {
      if (strct1->MatchMrrts[0] > 0.1)
      {
        strct2 = &opt.catalog[1].strctProps[strct1->MatchIDs[0]];
        fprintf (f, "%7d   ", strct1->ID);
        fprintf (f, "%7d   ", strct2->ID);
        fprintf (f, "%7d   ", isolated[strct1->ID]);
        fprintf (f, "%e    ", strct1->MatchMrrts[0]);
        fprintf (f, "%7d   ", strct1->Type);
        fprintf (f, "%e    ", strct1->TotMass);
        fprintf (f, "%e    ", strct2->TotMass);
        fprintf (f, "%e    ", strct1->Rsize);
        fprintf (f, "%e    ", strct2->Rsize*1000);
        fprintf (f, "%e    ", strct1->Vdisp);
        fprintf (f, "%e    ", strct2->Vdisp);
        fprintf (f, "%e    ", strct1->Pos[0]);
        fprintf (f, "%e    ", strct2->Pos[0]);
        fprintf (f, "%e    ", strct1->Pos[1]);
        fprintf (f, "%e    ", strct2->Pos[1]);
        fprintf (f, "%e    ", strct1->Pos[2]);
        fprintf (f, "%e    ", strct2->Pos[2]);
        fprintf (f, "\n");
      }
    }
  }
  fclose (f);
  free (isolated);

  //
  //  Free memory
  //
  for (i = 0; i < opt.numCatalogs; i++)
    Catalog_free (&opt.catalog[i]);

  return (0);
}


//
//  Parameters
//
void ctlgMatch_params (Options * opt)
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

  // Catalogs to load
  fscanf (opt->param.file, "%d", &opt->numCatalogs);
  opt->catalog    = (Catalog *)    malloc (opt->numCatalogs * sizeof(Catalog));
  opt->simulation = (Simulation *) malloc (opt->numCatalogs * sizeof(Simulation));
  for (i = 0; i < opt->numCatalogs; i++)
  {
    fscanf (opt->param.file, "%s", buffer);  Archive_name   (&opt->catalog[i].archive, buffer);
                                             Archive_prefix (&opt->catalog[i].archive, buffer);
    fscanf (opt->param.file, "%s", buffer);  Archive_format (&opt->catalog[i].archive, buffer);
    fscanf (opt->param.file, "%s", buffer);  Archive_path   (&opt->catalog[i].archive, buffer);
    fscanf (opt->param.file, "%d", &dummy);  Archive_nfiles (&opt->catalog[i].archive, dummy);


    fscanf (opt->param.file, "%s", buffer);  Archive_name   (&opt->simulation[i].archive, buffer);
                                             Archive_prefix (&opt->simulation[i].archive, buffer);
    fscanf (opt->param.file, "%s", buffer);  Archive_format (&opt->simulation[i].archive, buffer);
    fscanf (opt->param.file, "%s", buffer);  Archive_path   (&opt->simulation[i].archive, buffer);
    fscanf (opt->param.file, "%d", &dummy);  Archive_nfiles (&opt->simulation[i].archive, dummy);
  }

  // TreeFrog input
  fscanf (opt->param.file, "%s", buffer);  Archive_name   (&opt->mtree, buffer);
                                           Archive_prefix (&opt->mtree, buffer);
  fscanf (opt->param.file, "%s", buffer);  Archive_format (&opt->mtree, buffer);
  fscanf (opt->param.file, "%s", buffer);  Archive_path   (&opt->mtree, buffer);
  fscanf (opt->param.file, "%d", &dummy);  Archive_nfiles (&opt->mtree, dummy);

  // ASCII Output
  fscanf (opt->param.file, "%s", buffer);  Archive_name   (&opt->output, buffer);
                                           Archive_prefix (&opt->output, buffer);
  fscanf (opt->param.file, "%s", buffer);  Archive_format (&opt->output, buffer);
  fscanf (opt->param.file, "%s", buffer);  Archive_path   (&opt->output, buffer);
  fscanf (opt->param.file, "%d", &dummy);  Archive_nfiles (&opt->output, dummy);
  fclose (opt->param.file);
}


//
//  Options
//
int ctlgMatch_options (int argc, char ** argv, Options * opt)
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
      	ctlgMatch_usage (0, argv);
        break;

      default:
      	ctlgMatch_usage (1, argv);
    }
  }

  if (flag == 0)
    ctlgMatch_usage (1, argv);

}


//
//  Usage
//
void ctlgMatch_usage (int opt, char ** argv)
{
  if (opt == 0)
  {
    printf ("                                                                         \n");
    printf ("  ctlgmatch.c                                                            \n");
    printf ("                                                                         \n");
    printf ("  Author:            Rodrigo Can\\~as                                    \n");
    printf ("                                                                         \n");
    printf ("  Last edition:      29 - 11 - 2017                                      \n");
    printf ("                                                                         \n");
    printf ("  Description:       This code cross matches two catalogs of structures  \n");
    printf ("                     (e.g. dark matter haloes, galaxies, streams, etc)   \n");
    printf ("                     of the same snapshot, and produces a list of the    \n");
    printf ("                     matched structures' properties.                     \n");
    printf ("                                                                         \n");
    printf ("                     This is meant to be useful to compare a catalog     \n");
    printf ("                     produced by e.g. VELOCIraptor and HaloMaker.        \n");
    printf ("                                                                         \n");
    printf ("                     Currently the code supports the following formats:  \n");
    printf ("                                                                         \n");
    printf ("                         Finders:                                        \n");
    printf ("                                     VELOCIraptor  (Elahi+2011)          \n");
    printf ("                                     HaloMaker     (Tweed+2009)          \n");
    printf ("                                                                         \n");
    printf ("                         Simulations:                                    \n");
    printf ("                                                                         \n");
    printf ("                                     Gadget-2 ++   (Springel 2005)       \n");
    printf ("                                     RAMSES        (Teyssier 2002)       \n");
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
  for (i = 0; i < opt.numCatalogs; i++)
  {
    printf ("%s\n", opt.catalog[i].archive.name);
    printf ("%s\n", opt.catalog[i].archive.path);
    printf ("%s\n", opt.catalog[i].archive.format);
    printf ("%d\n", opt.catalog[i].nstruct);
    printf ("%d\n", opt.catalog[i].nprocs);
    printf ("%d\n", opt.catalog[i].iprops);
    for (j = 1; j <= 10; j++)
    {
      printf ("%5d   ", opt.catalog[i].strctProps[j].ID);
      printf ("%4d   ", opt.catalog[i].strctProps[j].DirectHostID);
      printf ("%4d   ", opt.catalog[i].strctProps[j].HostID);
      printf ("%4d   ", opt.catalog[i].strctProps[j].NumSubs);
      printf ("%4d   ", opt.catalog[i].strctProps[j].Type);
      printf ("%8d   ", opt.catalog[i].strctProps[j].NumPart);
      printf ("%e    ", opt.catalog[i].strctProps[j].TotMass);
      printf ("%e    ", opt.catalog[i].strctProps[j].Pos[0]);
      printf ("%e    ", opt.catalog[i].strctProps[j].Pos[1]);
      printf ("%e    ", opt.catalog[i].strctProps[j].Pos[2]);
      printf ("%e    ", opt.catalog[i].strctProps[j].Vel[0]);
      printf ("%e    ", opt.catalog[i].strctProps[j].Vel[1]);
      printf ("%e    ", opt.catalog[i].strctProps[j].Vel[2]);
      printf ("\n");
    }
  }
  */
