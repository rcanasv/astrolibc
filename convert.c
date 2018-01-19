/*
 *
 *  \file    convert.c
 *  \brief   This code converts from a given simulation or finder
 *           to another.
 *
 *           Currently supported conversions are:
 *
 *                    HaloMaker   --->   VELOCIraptor (ascii)
 *
 *
 */


#include "convert.h"
#include "stf.h"


int main (int argc, char ** argv)
{
  int i, j, k;
  Options opt;


  options_convert (argc, argv, &opt);

  //
  //  Read parameter file
  //
  Archive * input;
  Archive * output;

  opt.param.file = fopen (opt.param.name, "r");
  if (opt.param.file == NULL)
  {
   printf ("Couldn't open file  %s\n", opt.param.name);
   printf ("Exiting...\n");
   exit (0);
  }
  fscanf (opt.param.file, "%s", opt.type);

  if (strcmp(opt.type, "finder") == 0 || \
      strcmp(opt.type, "Finder") == 0)
  {
    input  = &opt.icatalog.archive;
    output = &opt.ocatalog.archive;
  }


  // Input file
  fscanf (opt.param.file, "%s", buffer);  Archive_name   (input, buffer);
                                          Archive_prefix (input, buffer);
  fscanf (opt.param.file, "%s", buffer);  Archive_format (input, buffer);
  fscanf (opt.param.file, "%s", buffer);  Archive_path   (input, buffer);

  // Output file
  fscanf (opt.param.file, "%s", buffer);  Archive_name   (output, buffer);
                                          Archive_prefix (output, buffer);
  fscanf (opt.param.file, "%s", buffer);  Archive_format (output, buffer);
  fscanf (opt.param.file, "%s", buffer);  Archive_path   (output, buffer);
  fclose (opt.param.file);


  if (strcmp(opt.type, "finder") ||
      strcmp(opt.type, "Finder"))
    convert_finder (&opt);




  return (0);
}



//
//
//
int convert_finder (Options * opt)
{
  int i, j, k;

  //
  //  Load catalog to convert
  //
  Catalog_init (&opt->icatalog);
  Catalog_load (&opt->icatalog);

  //
  //  Copy info to output catalog
  //
  opt->ocatalog.nstruct    = opt->icatalog.nstruct;
  opt->ocatalog.nprocs     = 1;
  opt->ocatalog.iprops     = 1;
  opt->ocatalog.iparts     = 1;
  opt->ocatalog.strctProps = opt->icatalog.strctProps;
  opt->ocatalog.strctParts = opt->icatalog.strctParts;


  stf_write_catalog_group  (&opt->ocatalog);
  stf_write_catalog_particles (&opt->ocatalog);

  //stf_write_group_catalogs  (&opt->ocatalog);
  //stf_write_group_particles (&opt->ocatalog);

/*
  Catalog * ctlg = &opt->ocatalog;
  for (i = 317; i <= 317; i++)
  {
    printf ("%d\n", ctlg->strctProps[i].ID);
    for (j = 0; j < ctlg->strctProps[i].NumPart; j++)
      printf ("     %u \n", -1*ctlg->strctProps[i].PIDs[j]);
  }
*/

  Catalog_free (&opt->icatalog);
}



//
//  Options
//
int options_convert (int argc, char ** argv, Options * opt)
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
     	 usage_convert (0, argv);
       break;

     default:
     	 usage_convert (1, argv);
   }
 }

 if (flag == 0)
   usage_convert (1, argv);

}


//
//  Usage
//
void usage_convert (int opt, char ** argv)
{
 if (opt == 0)
 {
   printf ("                                                                         \n");
   printf ("  convert.c                                                              \n");
   printf ("                                                                         \n");
   printf ("  Author:            Rodrigo Can\\~as                                    \n");
   printf ("                                                                         \n");
   printf ("  Last edition:      18 - 01 - 2018                                      \n");
   printf ("                                                                         \n");
   printf ("  Description:       This code converts file formats from Simulations    \n");
   printf ("                     or finders                                          \n");
   printf ("                                                                         \n");
   printf ("                     Currently the code supports the following formats:  \n");
   printf ("                                                                         \n");
   printf ("                         Finders:                                        \n");
   printf ("                                     VELOCIraptor  (Elahi+2011)          \n");
   printf ("                                     HaloMaker     (Tweed+2009)          \n");
   printf ("                                                                         \n");
/*   printf ("                         Simulations:                                    \n");
   printf ("                                                                         \n");
   printf ("                                     Gadget-2 ++   (Springel 2005)       \n");
   printf ("                                     RAMSES        (Teyssier 2002)       \n");
   printf ("                                                                         \n");*/
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
