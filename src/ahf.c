/*
 *  \file ahf.c
 *  \brief This file contains AHF functions.
 *
 */

#include "ahf.h"

//---------- Read Propertes File ----------//

void ahf_read_properties (Catalog * ahf)
{
  int    i, j, k;

  int       dummyi;
  long      dummyl;
  float     dummyf;
  double    dummyd;

  char    * buffc;
  int     * buffi;
  long    * buffl;
  float   * bufff;
  double  * buffd;

  FILE * f;
  char   parts_fname  [LONG_LENGTH];
  char   halos_fname  [LONG_LENGTH];
  char   longbuffer   [LONG_LENGTH];
  int    mystructs;

  int    oTask;

  //
  // Open properties file to read total number of structures
  // and number of processors if stf was run with MPI
  //
  ahf->nstruct = 0;
  sprintf (parts_fname, "%s/%s.AHF_particles", ahf->archive.path, ahf->archive.prefix);
  if ((f = fopen (parts_fname, "r")) != NULL)
  {
    fgets  (longbuffer, sizeof(longbuffer), f);
    sscanf (longbuffer, "%d", &(ahf->nstruct));
    fclose (f);
  }

  // Try bz2 compressed file
  sprintf (parts_fname, "%s/%s.AHF_particles.bz2", ahf->archive.path, ahf->archive.prefix);
  if ((f = fopen (parts_fname, "r")) != NULL)
  {
    //read_bz2_file (parts_fname, &buff, 0);
    fclose (f);
    read_bz2_file (parts_fname, &buffc, 20);
    sscanf (buffc, "%d", &(ahf->nstruct));
  }

  // Number of particles should have already been read, otherwise error
  if (ahf->nstruct == 0)
  {
    printf ("ERROR: Cannot open particles file  %s\n", parts_fname);
    exit (0);
  }
  printf ("nstruct  %d\n", ahf->nstruct);

  // Allocate memory for structure properties
  if ( (ahf->strctProps = (Structure *) malloc ((ahf->nstruct+1) * sizeof(Structure))) == NULL)
  {
    printf ("Couldn't allocate memory for structure properties.\n");
    printf ("Exiting...\n");
    exit (0);
  }
  else
    ahf->iprops = 1;

  for (i = 1; i <= ahf->nstruct; i++)
  {
    // Pointers
    ahf->strctProps[i].SubIDs     = NULL;
    ahf->strctProps[i].MatchIDs   = NULL;
    ahf->strctProps[i].MatchMrrts = NULL;
    ahf->strctProps[i].Part       = NULL;
    // Flags
    ahf->strctProps[i].flg_PartRadius           = 0;
    ahf->strctProps[i].flg_SortedByRadius       = 0;
    ahf->strctProps[i].flg_CorrectedPeriodicity = 0;
    ahf->strctProps[i].flg_ShiftedCM            = 0;
  }

  //
  // Read file(s) and store desired structure properties
  //
  int offst = 1;

  ahf->nprocs = 1;
  mystructs = ahf->nstruct;
  oTask = 0;
  sprintf (parts_fname, "%s/%s.AHF_halos", ahf->archive.path, ahf->archive.prefix);
  if ((f = fopen (parts_fname, "r")) == NULL)
  {
    sprintf (parts_fname, "%s/%s.AHF_halos.0", ahf->archive.path, ahf->archive.prefix);
    if ((f = fopen (parts_fname, "r")) == NULL)
    {
      printf ("ERROR: Cannot open file  %s\n", parts_fname);
      exit (0);
    }
  }


  for (i = 0; i < ahf->nprocs; i++)
  {
    fgets  (longbuffer, LONG_LENGTH, f); // AHF Column names
    for (j = 0; j < mystructs; j++)
    {
      fgets (longbuffer, LONG_LENGTH, f);
      ahf->strctProps[j+offst].oTask = oTask;
      sscanf (longbuffer, "%ld  %ld  %d  %lf  %d  %lf  %lf  %lf  %lf  %lf  %lf   \
              %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  \
              %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  \
              %lf  %lf  %lf  %lf  %lf  %lf",                                                                          \
              &(ahf->strctProps[j+offst].ID),            \
              &(ahf->strctProps[j+offst].HostID),        \
              &(ahf->strctProps[j+offst].NumSubs),       \
              &(ahf->strctProps[j+offst].Mvir),          \
              &(ahf->strctProps[j+offst].NumPart),       \
              &(ahf->strctProps[j+offst].Pos[0]),        \
              &(ahf->strctProps[j+offst].Pos[1]),        \
              &(ahf->strctProps[j+offst].Pos[2]),        \
              &(ahf->strctProps[j+offst].Vel[0]),        \
              &(ahf->strctProps[j+offst].Vel[1]),        \
              &(ahf->strctProps[j+offst].Vel[2]),        \
              &(ahf->strctProps[j+offst].Rvir),          \
              &(ahf->strctProps[j+offst].Rmax),          \
              &dummyd,                                   \
              &(ahf->strctProps[j+offst].mbpOffset),     \
              &(ahf->strctProps[j+offst].comOffset),     \
              &(ahf->strctProps[j+offst].Vmax),          \
              &dummyd,                                   \
              &(ahf->strctProps[j+offst].Vdisp),         \
              &dummyd,                                   \
              &dummyd,                                   \
              &dummyd,                                   \
              &dummyd,                                   \
              &dummyd,                                   \
              &dummyd,                                   \
              &dummyd,                                   \
              &dummyd,                                   \
              &dummyd,                                   \
              &dummyd,                                   \
              &dummyd,                                   \
              &dummyd,                                   \
              &dummyd,                                   \
              &dummyd,                                   \
              &dummyd,                                   \
              &dummyd,                                   \
              &dummyd,                                   \
              &dummyd,                                   \
              &dummyd,                                   \
              &(ahf->strctProps[j+offst].Ekin),          \
              &(ahf->strctProps[j+offst].Epot),          \
              &dummyd,                                   \
              &dummyd,                                   \
              &(ahf->strctProps[j+offst].cNFW)           \
            );
    }
    offst += mystructs;
    fclose(f);
  }
}


void ahf_catalog_get_particle_list (Catalog * ahf)
{
  int     i, j, k;
  int     check;
  int     npart;
  char  * fbuff;

  int       dummyi;
  long      dummyl;
  float     dummyf;
  double    dummyd;

  FILE * f = NULL;
  char   parts_fname  [LONG_LENGTH];
  char   halos_fname  [LONG_LENGTH];
  char   longbuffer   [LONG_LENGTH];
  int    mystructs;

  int    oTask;

  Structure * strct;

  // Allocate memory for PIDS
  for (i = 1; i <= ahf->nstruct; i++)
  {
    strct = &ahf->strctProps[i];
    if ((strct->PIDs = (long *) malloc (strct->NumPart * sizeof(long)))==NULL)
    {
      printf ("Couldn't allocate memory for Struct Particle IDs\n");
      exit (0);
    }
    strct->iIDs = 1;
  }

  // Open particles file
  sprintf (parts_fname, "%s/%s.AHF_particles", ahf->archive.path, ahf->archive.prefix);
  if ((f = fopen (parts_fname, "r")) != NULL)
  {
    // Read
    fgets  (longbuffer, NAME_LENGTH, f); // Number of structures

    // Loop over structures to get Particle IDs
    for (k = 1; k <= ahf->nstruct; k++)
    {
      fgets (longbuffer, NAME_LENGTH, f);
      check = sscanf (longbuffer, "%d  %ld", &npart, &dummyl);
      strct = &ahf->strctProps[k];
      for (i = 0; i < npart; i++)
      {
        fgets (longbuffer, NAME_LENGTH, f);
        check = sscanf (longbuffer, "%ld  %d", &strct->PIDs[i], &dummyi);
      }
    }
    fclose (f);
  }

  // Try bz2 compressed file
  sprintf (parts_fname, "%s/%s.AHF_particles.bz2", ahf->archive.path, ahf->archive.prefix);
  if ((f = fopen (parts_fname, "r")) != NULL)
  {
    fclose (f);
    read_bz2_file (parts_fname, &fbuff, 0);
    void * init = fbuff;

    // Read
    mysgets (longbuffer, NAME_LENGTH, &fbuff); // Number of structures

    // Loop over structures to get Particle IDs
    for (k = 1; k <= ahf->nstruct; k++)
    {
      mysgets (longbuffer, NAME_LENGTH, &fbuff);
      check = sscanf (longbuffer, "%d  %ld", &npart, &dummyl);
      strct = &ahf->strctProps[k];
      for (i = 0; i < npart; i++)
      {
        mysgets (longbuffer, NAME_LENGTH, &fbuff);
        check = sscanf (longbuffer, "%ld  %d", &strct->PIDs[i], &dummyi);
      }
    }

    fbuff = init;
    free (fbuff);
  }

  // Number of particles should have already been read, otherwise error
  if (f == NULL)
  {
    printf ("ERROR: Cannot open particles file  %s\n", parts_fname);
    exit (0);
  }

}
