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
    stf->strctProps[i].SubIDs     = NULL;
    stf->strctProps[i].MatchIDs   = NULL;
    stf->strctProps[i].MatchMrrts = NULL;
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



void stf_write_catalog_group (Catalog * stf)
{
  int    i, j, k;
  FILE * f;
  char   fname[NAME_LENGTH];

  Catalog * ctlg = stf;

  sprintf (fname, "%s/%s.catalog_groups.0", stf->archive.path, stf->archive.prefix);
  f = fopen (fname, "w");

  // ThisTask   NProcs
  fprintf (f, "%d %d\n", 0, 1);

  // NThisTask  NTot
  fprintf (f, "%d %d\n", stf->nstruct, stf->nstruct);

  // Structure number of particles
  for (i = 1; i <= stf->nstruct; i++)
    fprintf (f, "%d\n", stf->strctProps[i].NumPart);

  // Offset bound
  stf->strctProps[i].dummyi = 0;
  for (i = 2; i <= stf->nstruct; i++)
    stf->strctProps[i].dummyi = stf->strctProps[i-1].NumPart + stf->strctProps[i-1].dummyi;

  for (i = 1; i <= stf->nstruct; i++)
    fprintf (f, "%d\n", stf->strctProps[i].dummyi);

  // Offset unbound
  for (i = 1; i <= stf->nstruct; i++)
    fprintf (f, "%d\n", 0);

  fclose (f);
}


void stf_write_catalog_particles (Catalog * stf)
{
  int    i, j, k;
  FILE * f;
  char   fname[NAME_LENGTH];

  Catalog * ctlg = stf;

  int    npartTot = 0;

  //
  // Write catalog_particles
  //
  sprintf (fname, "%s/%s.catalog_particles.0", stf->archive.path, stf->archive.prefix);
  f = fopen (fname, "w");

  for (i = 1; i <= stf->nstruct; i++)
    npartTot += stf->strctProps[i].NumPart;

  // ThisTask   NProcs
  fprintf (f, "%d %d\n", 0, 1);

  // NThisTask  NTot
  fprintf (f, "%d %d\n", npartTot, npartTot);

  // Structure number of particles
  for (i = 1; i <= stf->nstruct; i++)
    for (j = 0; j < stf->strctProps[i].NumPart; j++)
      fprintf (f, "%u\n", -1*stf->strctProps[i].PIDs[j]);

  fclose (f);


  //
  // Write catalog_particles.unbound
  //
  sprintf (fname, "%s/%s.catalog_particles.unbound.0", stf->archive.path, stf->archive.prefix);
  f = fopen (fname, "w");

  // ThisTask   NProcs
  fprintf (f, "%d %d\n", 0, 1);

  // NThisTask  NTot
  fprintf (f, "%d %d\n", 0, 0);

  fclose (f);
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
*/



void stf_read_treefrog (Archive * tfrog, Catalog * stf)
{
  int i, j, k;
  FILE * f;

  int tmpid;
  int nmatch;

  char buffer [LONG_LENGTH];

  f = fopen (tfrog->name, "r");

  fgets (buffer, LONG_LENGTH, f);
  fgets (buffer, LONG_LENGTH, f);
  fgets (buffer, LONG_LENGTH, f);
  fgets (buffer, LONG_LENGTH, f);

  int    * MatchIDs;
  double * MatchMrrts;

  for (i = 1; i <= stf->nstruct; i++)
  {
    fgets (buffer, LONG_LENGTH, f);
    sscanf (buffer, "%d  %d", &tmpid, &nmatch);

    stf->strctProps[i].NumMatch = nmatch;

    if (nmatch > 0)
    {
      stf->strctProps[i].MatchIDs   = (int *)   malloc (nmatch * sizeof(int));
      stf->strctProps[i].MatchMrrts = (float *) malloc (nmatch * sizeof(float));
      for (j = 0; j < nmatch; j++)
      {
        fgets  (buffer, LONG_LENGTH, f);
        sscanf (buffer, "%d  %f", &stf->strctProps[i].MatchIDs[j],    \
                                  &stf->strctProps[i].MatchMrrts[j]);
      }
    }
  }

  fclose (f);
}



void  stf_catalog_get_particle_properties (Catalog * stf, Simulation * sim)
{

  int    i, j, k;
  FILE * f;
  char   fname [NAME_LENGTH];


  stfExtendedOutput * xtndd;
  int                 ninextended;
  int                 id;
  int                 indx;


  Particle  * part;
  Structure * strct;


  for (i = 1; i <= stf->nstruct; i++)
  {
    strct = &stf->strctProps[i];

    strct->dummyi = 0;
    strct->iPart  = 0;
    strct->Part   = NULL;

    if (strct->Part = (Particle *) malloc (stf->nstruct * sizeof(Particle)) == NULL)
    {
      printf ("Error not enough memory to allocate Particle to Structure\n");
      exit(0);
    }
    else
      strct->iPart = 1;
  }


  sprintf (fname, "%s/%s.filesofgroup", stf->archive.path, stf->archive.name);
  if (f = fopen(fname,"r") != NULL)
  {
    fclose (f);

    //
    // Go through every file of simulation/extendedOutput
    //
    for (i = 0; i < sim->archive.nfiles; i++)
    {
      ninextended = stf_catalog_load_extended_output (stf, i, &xtndd);

      if (ninextended)
      {
        Simulation_load_particles (sim, i, &part);

        for (j = 0; j < ninextended; j++)
        {
          id    = xtndd.IdStruct[j];
          indx  = xtndd.oIndex[j];
          strct = &stf->strctProps[id];

          strct->Part[strct->dummyi++] = part[j];
        }
        free (xtndd);
        free (part);
      }
    }
  }
  else
  {
    // 1.  Read .catalog_* files

    // 2.  Open Simulation file and assign particles
    ;
  }
}



void  stf_structure_get_particle_properties (Structure * strct, Archive * arx)
{

  FILE * f;
  char fname [NAME_LENGTH];

  sprintf (fname, "%s/%s.filesofgroup", stf->archive.path, stf->archive.name);
  if (f = fopen(fname,"r") != NULL)
  {
    fclose (f);
    // Extended Output files should (in principle) exist
    stf_structure_load_extended_output (strct, arx);
  }
  else
  {
    // 1.  Read .catalog_* files

    // 2.  Open Simulation file and assign particles
  }
  return;
}





int stf_load_extended_output (Catalog * stf,  int filenum, stfExtendedOutput ** xtndd)
{

  int    i, j, k;
  FILE * f;
  char   fname  [NAME_LENGTH];
  char   buffer [NAME_LENGTH];
  int    nparts;

  stfExtendedOutput * extended;

  nparts = 0;

  sprintf(fname, "%s/%s.extended.%d", stf->archive.path, stf->archive.prefix, filenum);
  if ((f = fopen(fname, "r")) != NULL)
  {
    while (fgets(buffer, NAME_LENGTH, f) != NULL)
      nparts++;
    rewind(f);

    if (nparts > 0)
    {
      extended = (stfExtendedOutput *) malloc (nparts * sizeof(stfExtendedOutput));
      for (i = 0; i < nparts; i++)
      {
        fgets(buffer, NAME_LENGTH, f);
        sscanf(buffer, "%d  %d  %d  %d  ",                   \
                       &extended[i].oIndex, &extended[i].IdStruct, \
                       &extended[i].IdHost, &extended[i].IdIGM);
      }
      *(xtndd) = extended;
    }
    else
      extended = NULL;

    fclose (f);
  }
  else
  {
    printf ("Couldn't open file  %s\n", fname);
    exit (0);
  }

  return nparts;
}
