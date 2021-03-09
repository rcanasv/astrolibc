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

#include "astrolibc.h"
#include "allvars.h"
#include "stf.h"
#include "halomaker.h"
#include "ramses.h"


typedef struct Options
{
  int            verbose;
  Archive        param;
  Archive        output;
  Catalog        catalog;
} Options;


int   find_options (int argc, char ** argv, Options * opt);
void  find_usage   (int opt,  char ** argv);
void  find_params  (Options * opt);


int main (int argc, char ** argv)
{

  int          i, j, k, n;
  Options      opt;
  Structure  * strct;
  Structure  * strct1;
  Structure  * strct2;

  find_options (argc, argv, &opt);
  find_params  (&opt);

  //
  //  Load catalogs
  //
  Catalog_init (&opt.catalog);
  Catalog_load (&opt.catalog);


  FILE * f_csv;
  FILE * f_out;
  
  char   buffer [3000];
  char   data   [3000];
  char   numstr [3000];
  int    length;
  double nums [4];
  double dx, dy, dz, r;
  
  
  f_csv = fopen("illustris_oddballs_ids.txt", "r");
  f_out = fopen("illustris_oddballs_properties.txt", "w");
  fprintf (f_out, "     ID_ctlg");
  fprintf (f_out, "     ID_veloci");
  fprintf (f_out, "          Type");
  fprintf (f_out, "          Mass");
  fprintf (f_out, "       NumPart");
  fprintf (f_out, "         R_100");
  fprintf (f_out, "          D_CM");
  fprintf (f_out, "\n");

  while (fgets(buffer, sizeof(buffer), f_csv))
  {
    sscanf (buffer, "%s", data);
    length = strlen(data);
    for (i = 0, k = 0, j = 0; i < length; i++)
    {
      if (data[i] == ',' || i==length-1)
      {
        numstr[j] = '\0';
        j = 0;
        
        nums[k] = atof(numstr);
        k++;
      }
      else
      {
        numstr[j] = data[i];
        j++;
      }
    }
    for (i = 1; i < 4; i++)
    {
      nums[i] /= 0.704;
      printf ("%f ", nums[i]);
    }
    printf("\n");
    
    for (i = 1; i <= opt.catalog.nstruct; i++)
    {
      strct = &opt.catalog.strctProps[i];
      if (strct->Type > 7)
      {
        dx = strct->Pos[0] - nums[1];
        dy = strct->Pos[1] - nums[2];
        dz = strct->Pos[2] - nums[3];
        r = dx*dx + dy*dy + dz*dz;
        r = sqrt(r);
        
        if (r < 5) // in kpc
        {
          fprintf (f_out, "%12d  ", (int) nums[0]);
          fprintf (f_out, "%12d  ", strct->ID);
          fprintf (f_out, "%12d  ", strct->Type);
          fprintf (f_out, "%e  ", strct->TotMass*1e10);
          fprintf (f_out, "%12d  ", strct->NumPart);          
          fprintf (f_out, "%e  ", strct->Rsize*1000.0);
          fprintf (f_out, "%e  ", r);
          fprintf (f_out, "\n");
        }
      }
    }
  }
  fclose(f_csv);
  fclose(f_out);
 

  //
  //  Free memory
  //
  Catalog_free (&opt.catalog);

  return (0);
}


//
//  Parameters
//
void find_params (Options * opt)
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

  fscanf (opt->param.file, "%s  %s  %s  %s  %d", prfxbuff,namebuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->catalog.archive, namebuff);
  Archive_prefix (&opt->catalog.archive, prfxbuff);
  Archive_format (&opt->catalog.archive, frmtbuff);
  Archive_path   (&opt->catalog.archive, pathbuff);
  Archive_nfiles (&opt->catalog.archive, nflsbuff);

  // ASCII Output
  fscanf (opt->param.file, "%s  %s  %s  %s  %d", prfxbuff,namebuff, frmtbuff, pathbuff, &nflsbuff);
  Archive_name   (&opt->output, namebuff);
  Archive_prefix (&opt->output, prfxbuff);
  Archive_format (&opt->output, frmtbuff);
  Archive_path   (&opt->output, pathbuff);
  Archive_nfiles (&opt->output, nflsbuff);

  fclose (opt->param.file);
}


//
//  Options
//
int find_options (int argc, char ** argv, Options * opt)
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
    {"aperture",0, NULL, 'a'},
    {0,         0, NULL, 0}
  };

  while ((myopt = getopt_long (argc, argv, "i:a:vh", lopts, &index)) != -1)
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
        find_usage (0, argv);
        break;

      default:
        find_usage (1, argv);
    }
  }

  if (flag == 0)
    find_usage (1, argv);

}


//
//  Usage
//
void find_usage (int opt, char ** argv)
{
  if (opt == 0)
  {
    printf ("                                                                         \n");
    printf ("  find.c                                                                 \n");
    printf ("                                                                         \n");
    printf ("  Author:            Rodrigo Can\\~as                                    \n");
    printf ("                                                                         \n");
    printf ("  Last edition:      7 - 5 - 2019                                        \n");
    printf ("                                                                         \n");
    printf ("  Description:                                                           \n");
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