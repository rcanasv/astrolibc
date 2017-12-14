/*
 *  \file stf.c
 *  \brief This file contains VELOCIraptor-stf functions.
 *
 */

#include "particle.h"
#include "structure.h"
#include "catalog.h"
#include "stf.h"


//---------- Read Propertes File ----------//

void stf_read_properties (Catalog * stf)
{
  int    i, j, k;

  int    dummyi;
  long   dummyl;
  float  dummyf;
  double dummyd;

  FILE * f;
  char   propts_fname [NAME_LENGTH];
  int    mystructs;

  //
  // Open properties file to read total number of structures
  // and number of processors if stf was run with MPI
  //
  sprintf (propts_fname, "%s/%s.properties", stf->archive.path, stf->archive.prefix);
  if ((f = fopen (propts_fname, "r")) == NULL)
  {
    sprintf (propts_fname, "%s/%s.properties.0", stf->archive.path, stf->archive.prefix);
    if ((f = fopen (propts_fname, "r")) == NULL)
    {
      printf ("ERROR: Cannot open file  %s\n", propts_fname);
      exit (0);
    }
  }
  fgets  (longbuffer, sizeof(longbuffer), f);  sscanf (longbuffer, "%d  %d", &dummyi, &(stf->nprocs));
  fgets  (longbuffer, sizeof(longbuffer), f);  sscanf (longbuffer, "%d  %d", &dummyi, &(stf->nstruct));
  fclose (f);

  //
  // Allocate memory for structure properties
  //
  if ( (stf->strctProps = (Structure *) malloc ((stf->nstruct+1) * sizeof(Structure))) == NULL)
  {
    printf ("Couldn't allocate memory for structure properties.\n");
    printf ("Exiting...\n");
    exit (0);
  }
  else
    stf->iprops = 1;

  for (i = 1; i <= stf->nstruct; i++)
  {
    stf->strctProps[i].SubIDs    = NULL;
    stf->strctProps[i].ProgIDs   = NULL;
    stf->strctProps[i].ProgMrrts = NULL;
  }

  //
  // Read file(s) and store desired structure properties
  //
  int offst = 1;

  for (i = 0; i < stf->nprocs; i++)
  {
    if (stf->nprocs == 1)
      sprintf (propts_fname, "%s/%s.properties", stf->archive.path, stf->archive.prefix);
    else
      sprintf (propts_fname, "%s/%s.properties.%d", stf->archive.path, stf->archive.prefix, i);

    //printf ("Openning file  %s  \n", propts_fname);

    if ((f = fopen (propts_fname, "r")) == NULL)
    {
      printf ("ERROR: Cannot open file  %s", propts_fname);
      exit (0);
    }

    fgets  (longbuffer, NAME_LENGTH, f);
    fgets  (longbuffer, NAME_LENGTH, f);
    sscanf (longbuffer, "%d  %d", &mystructs, &dummyi);
    fgets  (longbuffer, 3000, f);

    for (j = 0; j < mystructs; j++)
    {
      fgets (longbuffer, 3000, f);
      sscanf (longbuffer, "%d  %d  %d  %d  %d  %d  %d  %lf  %lf  %lf  %lf  %lf  %lf  %lf          \
              %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf      \
              %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf      \
              %lf  %lf",                                                                          \
              &(stf->strctProps[j+offst].ID), &dummyi, &(stf->strctProps[j+offst].DirectHostID),  \
              &(stf->strctProps[j+offst].HostID), &(stf->strctProps[j+offst].NumSubs),            \
              &(stf->strctProps[j+offst].Type), &(stf->strctProps[j+offst].NumPart), &dummyd,     \
              &(stf->strctProps[j+offst].Pos[0]), &(stf->strctProps[j+offst].Pos[1]),             \
              &(stf->strctProps[j+offst].Pos[2]), &dummyd, &dummyd, &dummyd,                      \
              &(stf->strctProps[j+offst].Vel[0]), &(stf->strctProps[j+offst].Vel[1]),             \
              &(stf->strctProps[j+offst].Vel[2]), &dummyd, &dummyd, &dummyd,                      \
              &(stf->strctProps[j+offst].TotMass), &dummyd, &dummyd, &dummyd, &dummyd,            \
              &(stf->strctProps[j+offst].Efrac), &dummyd, &(stf->strctProps[j+offst].Rsize),      \
              &dummyd, &dummyd, &dummyd, &(stf->strctProps[j+offst].RHalfMass),                   \
              &(stf->strctProps[j+offst].Rvmax), &(stf->strctProps[j+offst].Vmax),                \
              &(stf->strctProps[j+offst].Vdisp), &dummyd, &dummyd, &dummyd, &dummyd, &dummyd,     \
              &dummyd, &dummyd, &dummyd, &dummyd, &(stf->strctProps[j+offst].Lambda),             \
              &(stf->strctProps[j+offst].L[0]), &(stf->strctProps[j+offst].L[1]),                 \
              &(stf->strctProps[j+offst].L[2])                                                    \
            );
    }
    offst += mystructs;
    fclose(f);
  }
}


/*

int read_stf_filesofgroup (char * prefix, int strct_id, int ** files_of_strct)
{
  int i, j, k;
  FILE * f;

  char buffer [NAME_LENGTH];

  sprintf (buffer, "%s.filesofgroup", prefix);

  int tmpid;
  int nfiles;

  f = fopen (buffer, "r");
  do
  {
    fgets (buffer, NAME_LENGTH, f);
    sscanf (buffer, "%d  %d", &tmpid, &nfiles);
    fgets (buffer, NAME_LENGTH, f);
    get_n_num_from_string (buffer, nfiles, files_of_strct);
  }
  while (tmpid != strct_id);

  fclose (f);

  return nfiles;
}




int load_treefrog (char * tffile, int strct_id, int ** prog_ids, float ** prog_mrrts)
{
  int i, j, k;
  FILE * f;

  int tmpid;
  int nprogs;

  char buffer [LONG_LENGTH];

  f = fopen (tffile, "r");

  fgets (buffer, LONG_LENGTH, f);
  fgets (buffer, LONG_LENGTH, f);
  fgets (buffer, LONG_LENGTH, f);
  fgets (buffer, LONG_LENGTH, f);

  do
  {
    fgets (buffer, LONG_LENGTH, f);
    sscanf (buffer, "%d  %d", &tmpid, &nprogs);

    if (tmpid != strct_id)
      for (i = 0; i < nprogs; i++)
        fgets (buffer, LONG_LENGTH, f);
  }
  while (tmpid != strct_id);


  *prog_ids   = (int *)   malloc (nprogs * sizeof(int));
  *prog_mrrts = (float *) malloc (nprogs * sizeof(float));

  for (i = 0; i < nprogs; i++)
  {
    fgets (buffer, LONG_LENGTH, f);
    sscanf (buffer, "%d  %f", &prog_ids[0][i], &prog_mrrts[0][i]);
  }
  fclose (f);

  return nprogs;
}




int load_stf_extended_output (char * prefix, int filenum)
{
    FILE * f;
    char buffer[NAME_LENGTH];
    sprintf(buffer, "%s.extended.%d", prefix, filenum);
    int nparts = 0;
    int i;

    if ((f = fopen(buffer, "r")) == NULL)
      return 0;

    while (fgets(buffer, NAME_LENGTH, f) != NULL)
      nparts++;
    rewind(f);

//     printf("nparts %d\n", nparts);

    extended_oIndex   = (int *) malloc (nparts * sizeof(int));
    extended_IdStruct = (int *) malloc (nparts * sizeof(int));
    extended_IdHost   = (int *) malloc (nparts * sizeof(int));
    extended_IdIGM    = (int *) malloc (nparts * sizeof(int));

    for (i = 0; i < nparts; i++)
    {
      fgets(buffer, NAME_LENGTH, f);
      sscanf(buffer, "%d  %d  %d  %d  ", &extended_oIndex[i], &extended_IdStruct[i], &extended_IdHost[i], &extended_IdIGM[i]);
    }

    fclose (f);

    return nparts;
}

void free_extended_arrays (void)
{
  free (extended_oIndex);
  free (extended_IdStruct);
  free (extended_IdHost);
  free (extended_IdIGM);
}
*/