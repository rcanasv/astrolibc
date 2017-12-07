/*
 *  \file   halomakerio.c
 *  \brief  This file contains Gadget input-output routines
 *
 *
 */

#include "base.h"
#include "halomaker.h"
#include "halomakerio.h"


//---------- Read Propertes File ----------//

void halomaker_read_properties (Catalog * hmkr)
{
  int    i, j, k;

  int    dummyi;
  long   dummyl;
  float  dummyf;
  double dummyd;

  FILE * f;
  char   treebricks_fname [NAME_LENGTH];
  int    mystructs;

  //
  // Open properties file to read total number of structures
  // and number of processors if hmkr was run with MPI
  //
  sprintf (treebricks_fname, "%s/%s.properties", hmkr->archive.path, hmkr->archive.prefix);
  if ((f = fopen (treebricks_fname, "r")) == NULL)
  {
    sprintf (propts_fname, "%s/%s.properties.0", hmkr->archive.path, hmkr->archive.prefix);
    if ((f = fopen (propts_fname, "r")) == NULL)
    {
      printf ("ERROR: Cannot open file  %s\n", propts_fname);
      exit (0);
    }
  }

  SKIP  fread (&hmkr->nparts,   sizeof(int),   1, f);   SKIP
  SKIP  fread (&hmkr->massres,  sizeof(float), 1, f);   SKIP
  SKIP  fread (&hmkr->aexp,     sizeof(float), 1, f);   SKIP
  SKIP  fread (&hmkr->OmegaM,   sizeof(float), 1, f);   SKIP
  SKIP  fread (&hmkr->AgeUniv,  sizeof(float), 1, f);   SKIP
  SKIP  fread (&hmkr->nstructs, sizeof(int),   1, f);
        fread (&hmkr->nsubs,    sizeof(int),   1, f);   SKIP

  hmkr->ntotal = hmkr->nstructs + hmkr->nsubs;



  fgets  (longbuffer, sizeof(longbuffer), f);  sscanf (longbuffer, "%d  %d", &dummyi, &(hmkr->nprocs));
  fgets  (longbuffer, sizeof(longbuffer), f);  sscanf (longbuffer, "%d  %d", &dummyi, &(hmkr->nstruct));
  fclose (f);

  //
  // Allocate memory for structure properties
  //
  if ( (hmkr->strctProps = (Structure *) malloc ((hmkr->nstruct+1) * sizeof(Structure))) == NULL)
  {
    printf ("Couldn't allocate memory for structure properties.\n");
    printf ("Exiting...\n");
    exit (0);
  }
  else
    hmkr->iprops = 1;

  for (i = 1; i <= hmkr->nstruct; i++)
  {
    hmkr->strctProps[i].SubIDs    = NULL;
    hmkr->strctProps[i].ProgIDs   = NULL;
    hmkr->strctProps[i].ProgMrrts = NULL;
  }

  //
  // Read file(s) and store desired structure properties
  //
  int offst = 1;

  for (i = 0; i < hmkr->nprocs; i++)
  {
    if (hmkr->nprocs == 1)
      sprintf (propts_fname, "%s/%s.properties", hmkr->archive.path, hmkr->archive.prefix);
    else
      sprintf (propts_fname, "%s/%s.properties.%d", hmkr->archive.path, hmkr->archive.prefix, i);

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
              &(hmkr->strctProps[j+offst].ID), &dummyi, &(hmkr->strctProps[j+offst].DirectHostID),  \
              &(hmkr->strctProps[j+offst].HostID), &(hmkr->strctProps[j+offst].NumSubs),            \
              &(hmkr->strctProps[j+offst].Type), &(hmkr->strctProps[j+offst].NumPart), &dummyd,     \
              &(hmkr->strctProps[j+offst].Pos[0]), &(hmkr->strctProps[j+offst].Pos[1]),             \
              &(hmkr->strctProps[j+offst].Pos[2]), &dummyd, &dummyd, &dummyd,                      \
              &(hmkr->strctProps[j+offst].Vel[0]), &(hmkr->strctProps[j+offst].Vel[1]),             \
              &(hmkr->strctProps[j+offst].Vel[2]), &dummyd, &dummyd, &dummyd,                      \
              &(hmkr->strctProps[j+offst].TotMass), &dummyd, &dummyd, &dummyd, &dummyd,            \
              &(hmkr->strctProps[j+offst].Efrac), &dummyd, &(hmkr->strctProps[j+offst].Rsize),      \
              &dummyd, &dummyd, &dummyd, &(hmkr->strctProps[j+offst].RHalfMass),                   \
              &(hmkr->strctProps[j+offst].Rvmax), &(hmkr->strctProps[j+offst].Vmax),                \
              &(hmkr->strctProps[j+offst].Vdisp), &dummyd, &dummyd, &dummyd, &dummyd, &dummyd,     \
              &dummyd, &dummyd, &dummyd, &dummyd, &(hmkr->strctProps[j+offst].Lambda),             \
              &(hmkr->strctProps[j+offst].L[0]), &(hmkr->strctProps[j+offst].L[1]),                 \
              &(hmkr->strctProps[j+offst].L[2])                                                    \
            );
    }
    offst += mystructs;
    fclose(f);
  }
}




/*
//
//
// ---  Read Galfile --- //
//
//
void read_gal_file (char * filename)
{
  int     dummy;
  int     i, j;
  FILE *  f;

  if ((f = fopen (filename, "r")) == NULL)
  {
    printf ("Can't open file named   %s \n", filename);
    exit(0);
  }

  HMKR_SKIP  fread (&gal_number, sizeof(int),    1, f);  HMKR_SKIP
  HMKR_SKIP  fread (&gal_level,  sizeof(int),    1, f);  HMKR_SKIP
  HMKR_SKIP  fread (&gal_mass,   sizeof(double), 1, f);  HMKR_SKIP
  HMKR_SKIP  fread (&gal_pos,    sizeof(double), 3, f);  HMKR_SKIP
  HMKR_SKIP  fread (&gal_vel,    sizeof(double), 3, f);  HMKR_SKIP
  HMKR_SKIP  fread (&gal_ang,    sizeof(double), 3, f);  HMKR_SKIP
  HMKR_SKIP  fread (&nlist,      sizeof(int),    1, f);  HMKR_SKIP

/*
  printf ("\n");
  printf ("my_number   %d \n", gal_number);
  printf ("Level       %d \n", gal_level);
  printf ("Mass        %8.5f \n", gal_mass);
  printf ("Pos         %8.5f \t %8.5f \t %8.5f \n", gal_pos[0], gal_pos[1], gal_pos[2]);
  printf ("Vel         %8.5f \t %8.5f \t %8.5f \n", gal_vel[0], gal_vel[1], gal_vel[2]);
  printf ("Ang Mom     %8.5f \t %8.5f \t %8.5f \n", gal_ang[0], gal_ang[1], gal_ang[2]);
  printf ("Nlist       %d \n", nlist);
  printf ("\n");
//

  //
  // Data is stored in double
  //
  if (!(p = malloc (nlist * sizeof(struct pdata_d))))
  {
    printf ("Cannot allocate memory for particle information\n");
    printf ("Exiting\n");
    exit (0);
  }

  // Positions
  for (j = 0; j < 3; j++)
  {
    HMKR_SKIP
    if (dummy == 8 * nlist)
      for (i = 0; i < nlist; i++)
      {
        fread (&p[i].Pos[j], sizeof(double), 1, f);
        p[i].Pos[j] = p[i].Pos[j] * 1000;
      }
    HMKR_SKIP
  }

  // Velocities
  for (j = 0; j < 3; j++)
  {
    HMKR_SKIP
    if (dummy == 8 * nlist)
      for (i = 0; i < nlist; i++)
      	fread (&p[i].Vel[j], sizeof(double), 1, f);
    HMKR_SKIP
  }

  // Masses
  HMKR_SKIP
  if (dummy == 8 * nlist)
    for (i = 0; i < nlist; i++)
    {
      fread (&p[i].Mass, sizeof(double), 1, f);

      // converts to M/10**10 Msun
      p[i].Mass *= 10.0;
    }
  HMKR_SKIP

  // Ids
  HMKR_SKIP
  if (dummy == 4 * nlist)
    for (i = 0; i < nlist; i++)
      fread (&p[i].Id, sizeof(int), 1, f);
  HMKR_SKIP

  // Age
  HMKR_SKIP
  if (dummy == 8 * nlist)
    for (i = 0; i < nlist; i++)
      fread (&p[i].Age, sizeof(double), 1, f);
  HMKR_SKIP

  // Metallicity
  HMKR_SKIP
  if (dummy == 8 * nlist)
    for (i = 0; i < nlist; i++)
      fread (&p[i].Metal, sizeof(double), 1, f);
  HMKR_SKIP

  // Chemical Elements
  HMKR_SKIP
  if (dummy == 8 * nlist)
    for (i = 0; i < nlist; i++)
      fread (&p[i].Chem, sizeof(double), 1, f);
  HMKR_SKIP

  HMKR_SKIP HMKR_PSKIP
  fclose(f);
}
*/
