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


  options_ctlgMatch (argc, argv, &opt);


  //
  //  Read parameter file
  //
  opt.param.file = fopen (opt.param.name, "r");
  if (opt.param.file == NULL)
  {
    printf ("Couldn't open file  %s\n", opt.param.name);
    printf ("Exiting...\n");
    exit (0);
  }
  fscanf (opt.param.file, "%d", &opt.numCatalogs);
  opt.catalog = (Catalog *) malloc (opt.numCatalogs * sizeof(Catalog));
  for (i = 0; i < opt.numCatalogs; i++)
  {
    fscanf (opt.param.file, "%s", buffer);  Archive_name   (&opt.catalog[i].archive, buffer);
                                            Archive_prefix (&opt.catalog[i].archive, buffer);
    fscanf (opt.param.file, "%s", buffer);  Archive_format (&opt.catalog[i].archive, buffer);
    fscanf (opt.param.file, "%s", buffer);  Archive_path   (&opt.catalog[i].archive, buffer);
  }
  fclose (opt.param.file);


  //
  //  Load catalogs
  //
  for (i = 0; i < opt.numCatalogs; i++)
  {
    Catalog_init (&opt.catalog[i]);
    Catalog_load (&opt.catalog[i]);
  }


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
      printf ("%d   ", opt.catalog[i].strctProps[j].ID);
      printf ("%d   ", opt.catalog[i].strctProps[j].DirectHostID);
      printf ("%d   ", opt.catalog[i].strctProps[j].HostID);
      printf ("%d   ", opt.catalog[i].strctProps[j].NumSubs);
      printf ("%d   ", opt.catalog[i].strctProps[j].Type);
      printf ("%d   ", opt.catalog[i].strctProps[j].NumPart);
      printf ("%g   ", opt.catalog[i].strctProps[j].TotMass);
      printf ("\n");
    }
  }



  for (i = 0; i < opt.numCatalogs; i++)
    Catalog_free (&opt.catalog[i]);


  return (0);
}



//
//  Options
//
int options_ctlgMatch (int argc, char ** argv, Options * opt)
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
      	usage_ctlgMatch (0, argv);
        break;

      default:
      	usage_ctlgMatch (1, argv);
    }
  }

  if (flag == 0)
    usage_ctlgMatch (1, argv);

}


//
//  Usage
//
void usage_ctlgMatch (int opt, char ** argv)
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
