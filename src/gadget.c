/*
 *  \file gadget.c
 *  \brief This file contains Gadget- format actions
 *
 */

 #include "gadget.h"

void gadget_init (Simulation * gdt)
{
  FILE  * fd;
  int     n, k, pc, pc_new;
  int     dummy;
  int     ntot_withmasses;
  char    hname [4];
  char    fname [NAME_LENGTH];

  gheader header1;

  pc = 0;

  sprintf (fname, "%s/%s", gdt->archive.path, gdt->archive.prefix);
  if ((fd = fopen (fname, "r")) == NULL)
  {
    sprintf (fname, "%s/%s.0", gdt->archive.path, gdt->archive.prefix);
    if ((fd = fopen (fname, "r")) == NULL)
    {
      printf ("Couldn't open file %s ... Exiting\n", fname);
      exit (0);
    }
  }


  if (gdt->format == GADGET_HEAD)
  {
    fread(&dummy, sizeof(dummy), 1, fd);     //printf ("%d\n", dummy);
    fread(&hname[0], sizeof(hname), 1, fd);  //printf ("%s\n", hname);
    fread(&dummy, sizeof(dummy), 1, fd);     //printf ("%d\n", dummy);
    fread(&dummy, sizeof(dummy), 1, fd);     //printf ("%d\n", dummy);
  }

  fread(&dummy,   sizeof(dummy),   1, fd);  //printf ("%d\n", dummy);
  fread(&header1, sizeof(header1), 1, fd);
  fread(&dummy,   sizeof(dummy),   1, fd);  //printf ("%d\n", dummy);


  gdt->npart = 0;
  for(k = 0; k < 6; k++)
    gdt->npart += header1.npartTotal[k];

  gdt->Lbox = header1.BoxSize;
  for (k = 0; k < 6; k++)
  {
    gdt->MassTable[k]     = header1.mass[k];
    gdt->NpartThisFile[k] = header1.npart[k];
    gdt->NpartTot[k]      = header1.npartTotal[k];
  }
  gdt->NfilesPerSnapshot     = header1.num_files;
  gdt->Cooling               = header1.flag_cooling;
  gdt->Feedback              = header1.flag_feedback;
  gdt->SFR                   = header1.flag_sfr;
  gdt->h                     = header1.HubbleParam;
  gdt->HubbleParam           = gdt->h * 100.0;
  gdt->cosmology.HubbleParam = gdt->h * 100.0;
  gdt->cosmology.OmegaM      = header1.Omega0;
  gdt->cosmology.OmegaL      = header1.OmegaLambda;
  gdt->a                     = header1.time;
  gdt->z                     = header1.redshift;

  /*
  fread(&dummy, sizeof(dummy), 1, fd);     printf ("%d \n", dummy);
  fread(&hname, sizeof(dummy), 1, fd);     printf ("%s \n", hname); sprintf(hname,"    ");
  fread(&dummy, sizeof(dummy), 1, fd);     printf ("%d \n", dummy);
  fread(&dummy, sizeof(dummy), 1, fd);     printf ("%d \n", dummy);

  fread(&dummy, sizeof(dummy), 1, fd);     printf ("%d \n", dummy);
  for (k=0; k < 29; k++)
  {
    fread(&hname, sizeof(dummy), 1, fd);     printf ("%s \n", hname);sprintf(hname,"    ");
    fread(&hname, sizeof(double), 1, fd);    printf ("%s \n", hname);sprintf(hname,"    ");
    fread(&dummy, sizeof(dummy), 1, fd);     printf ("%d \n", dummy);
    fread(&dummy, sizeof(dummy), 1, fd);     printf ("%d \n", dummy);
    fread(&dummy, sizeof(dummy), 1, fd);     printf ("%d \n", dummy);
    fread(&dummy, sizeof(dummy), 1, fd);     printf ("%d \n", dummy);
    fread(&dummy, sizeof(dummy), 1, fd);     printf ("%d \n", dummy);
    fread(&dummy, sizeof(dummy), 1, fd);     printf ("%d \n", dummy);
    fread(&dummy, sizeof(dummy), 1, fd);     printf ("%d \n", dummy);
  }
  fread(&dummy, sizeof(dummy), 1, fd);     printf ("%d \n", dummy);

  fread(&dummy, sizeof(dummy), 1, fd);     printf ("%d \n", dummy);
  fread(&hname, sizeof(dummy), 1, fd);     printf ("%s \n", hname);sprintf(hname,"    ");
  fread(&dummy, sizeof(dummy), 1, fd);     printf ("%d \n", dummy);
  fread(&dummy, sizeof(dummy), 1, fd);     printf ("%d \n", dummy);
  */
  fclose (fd);
  return;
}

int gadget_get_npart_ThisFile (Simulation * gdt, int filenum)
{
  FILE   * fd;

  int      k;
  int      NumPart;
  int      dummy;

  gheader  header1;

  char     hname[4];
  char     fname[NAME_LENGTH];

  // Open File
  sprintf (fname, "%s/%s", gdt->archive.path, gdt->archive.prefix);
  if ((fd = fopen (fname, "r")) == NULL)
  {
    sprintf (fname, "%s/%s.%d", gdt->archive.path, gdt->archive.prefix, filenum);
    if ((fd = fopen (fname, "r")) == NULL)
    {
      printf ("Couldn't open file %s ... Exiting\n", fname);
      exit (0);
    }
  }

  // Read Header
  if (gdt->format == GADGET_HEAD)
  {
    fread(&dummy, sizeof(dummy), 1, fd);
    fread(&hname[0], sizeof(hname), 1, fd);
    fread(&dummy, sizeof(dummy), 1, fd);
    fread(&dummy, sizeof(dummy), 1, fd);
  }
  fread(&dummy,   sizeof(dummy),   1, fd);
  fread(&header1, sizeof(header1), 1, fd);
  fread(&dummy,   sizeof(dummy),   1, fd);

  fclose (fd);

  // Get total number of particles in THIS file
  for(k = 0, NumPart = 0; k < 6; k++)
    NumPart += header1.npart[k];

  return NumPart;
}

void gadget_load_particles (Simulation * gdt, int filenum, Particle ** part)
{
  FILE  * fd;
  int     n, k, pc, pc_new;
  int     NumPart;
  int     dummy;
  int     ntot_withmasses;
  char    fname [NAME_LENGTH];
  char    hname [4];

  int     i;
  float   ftmp3 [3];
  float   ftmp;

  Particle * P;
  gheader    header1;

  pc = 0;

  //
  // Open File
  //
  sprintf (fname, "%s/%s", gdt->archive.path, gdt->archive.prefix);
  if ((fd = fopen (fname, "r")) == NULL)
  {
    sprintf (fname, "%s/%s.%d", gdt->archive.path, gdt->archive.prefix, filenum);
    if ((fd = fopen (fname, "r")) == NULL)
    {
      printf ("Couldn't open file %s ... Exiting\n", fname);
      exit (0);
    }
  }

  //
  // Read Header
  //
  if (gdt->format == GADGET_HEAD)
  {
    fread(&dummy, sizeof(dummy), 1, fd);     //printf ("%d\n", dummy);
    fread(&hname[0], sizeof(hname), 1, fd);  //printf ("%s\n", hname);
    fread(&dummy, sizeof(dummy), 1, fd);     //printf ("%d\n", dummy);
    fread(&dummy, sizeof(dummy), 1, fd);     //printf ("%d\n", dummy);
  }

  fread(&dummy,   sizeof(dummy),   1, fd);  //printf ("%d\n", dummy);
  fread(&header1, sizeof(header1), 1, fd);
  fread(&dummy,   sizeof(dummy),   1, fd);  //printf ("%d\n", dummy);


  //
  // Allocate Memory
  //
  // Get total number of particles in THIS file
  for(k = 0, NumPart = 0; k < 6; k++)
    NumPart += header1.npart[k];

  // Get total number of particles with mass array in THIS file
  for(k = 0, ntot_withmasses = 0; k < 6; k++)
    if(header1.mass[k] == 0)
      ntot_withmasses += header1.npart[k];

  printf("Allocating memory...");
  if((*(part) = (Particle *) malloc (NumPart * sizeof(Particle))) == NULL)
  {
    fprintf(stderr, "failed to allocate memory.\n");
    exit(0);
  }
  P = *(part);
  printf("done\n");


  //
  // Read Positions
  //
  if (gdt->format == GADGET_HEAD)
  {
    fread(&dummy, sizeof(dummy), 1, fd);  //printf ("%d\n", dummy);
    fread(&hname, sizeof(hname), 1, fd);  //printf ("%s\n", hname);
    fread(&dummy, sizeof(dummy), 1, fd);  //printf ("%d\n", dummy);
    fread(&dummy, sizeof(dummy), 1, fd);  //printf ("%d\n", dummy);
  }

  fread(&dummy,   sizeof(dummy),   1, fd);  //printf ("%d\n", dummy);
  for(k = 0, pc_new = pc; k < 6; k++)
    for(n = 0; n < header1.npart[k]; n++)
    {
      //fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
      fread(&ftmp3[0], sizeof(float), 3, fd);
      for (i = 0; i < 3; i++)
         P[pc_new].Pos[i] = ftmp3[i] / header1.HubbleParam;
      pc_new++;
    }
  fread(&dummy,   sizeof(dummy),   1, fd);  //printf ("%d\n", dummy);


  //
  // Read Velocities
  //
  if (gdt->format == GADGET_HEAD)
  {
    fread(&dummy, sizeof(dummy), 1, fd);  //printf ("%d\n", dummy);
    fread(&hname, sizeof(hname), 1, fd);  //printf ("%s\n", hname);
    fread(&dummy, sizeof(dummy), 1, fd);  //printf ("%d\n", dummy);
    fread(&dummy, sizeof(dummy), 1, fd);  //printf ("%d\n", dummy);
  }
  fread(&dummy,   sizeof(dummy),   1, fd);  //printf ("%d\n", dummy);
  for(k = 0, pc_new = pc; k < 6; k++)
    for(n = 0; n < header1.npart[k]; n++)
    {
      //fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
      fread(&ftmp3[0], sizeof(float), 3, fd);
      for (i = 0; i < 3; i++)
         P[pc_new].Vel[i] = ftmp3[i];
      pc_new++;
    }
  fread(&dummy,   sizeof(dummy),   1, fd);  //printf ("%d\n", dummy);


  //
  // Read Ids
  //
  if (gdt->format == GADGET_HEAD)
  {
    fread(&dummy, sizeof(dummy), 1, fd);  //printf ("%d\n", dummy);
    fread(&hname, sizeof(hname), 1, fd);  //printf ("%s\n", hname);
    fread(&dummy, sizeof(dummy), 1, fd);  //printf ("%d\n", dummy);
    fread(&dummy, sizeof(dummy), 1, fd);  //printf ("%d\n", dummy);
  }
  fread(&dummy,   sizeof(dummy),   1, fd);  //printf ("%d\n", dummy);
  for(k = 0, pc_new = pc; k < 6; k++)
    for(n = 0; n < header1.npart[k]; n++)
    {
      fread(&P[pc_new].Id, sizeof(int), 1, fd);
      pc_new++;
    }
  fread(&dummy,   sizeof(dummy),   1, fd);  //printf ("%d\n", dummy);


  //
  // Read Masses
  //
  if (gdt->format == GADGET_HEAD)
  {
    fread(&dummy, sizeof(dummy), 1, fd);  //printf ("%d\n", dummy);
    fread(&hname, sizeof(hname), 1, fd);  //printf ("%s\n", hname);
    fread(&dummy, sizeof(dummy), 1, fd);  //printf ("%d\n", dummy);
    fread(&dummy, sizeof(dummy), 1, fd);  //printf ("%d\n", dummy);
  }
  if(ntot_withmasses>0)
    fread(&dummy,   sizeof(dummy),   1, fd);  //printf ("%d\n", dummy);
  for(k = 0, pc_new = pc; k < 6; k++)
  {
    for(n = 0; n < header1.npart[k]; n++)
    {
      P[pc_new].Type = k;
      if(header1.mass[k] == 0)
      {
        //fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
        fread(&ftmp, sizeof(float), 1, fd);
        P[pc_new].Mass = ftmp / header1.HubbleParam;
      }
      else
        P[pc_new].Mass= header1.mass[k] / header1.HubbleParam;
      pc_new++;
    }
  }
  if(ntot_withmasses>0)
    fread(&dummy,   sizeof(dummy),   1, fd);  //printf ("%d\n", dummy);


  //
  // Need to add here gas and star properties
  //

  // Close File
  printf("Successful snapshot reading\n");
  fclose(fd);

  return;
}


void gadget_write_snapshot (Particle * P, int NPartTot, gheader * header, Archive * output)
{
  FILE * snap_file;

  int dummy;
  int k, pc_new, n;
  int pc = 0;
  int ntot_withmasses;

  double      ref_mass [6];
  int         id_ref   [6];
  int     **  ids;

  int i, bob, offset;

  if((snap_file = fopen(output->name,"w")) == NULL)
  {
    printf("Couldn't open file %s\n", output->name);
    exit(0);
  }

  //! Initialize
  for(k = 0; k < 6; k++)
  {
    header->npart[k]      = 0;
    header->npartTotal[k] = 0;
    header->mass[k]       = 0;
    ref_mass[k]           = 0;
    id_ref[k]             = 0;
  }
  header->num_files       = 1;


  //! Count Particles by type
  for(k = 0; k < NPartTot; k++)
  {
    header->npart[P[k].Type]++;
    header->npartTotal[P[k].Type]++;
  }

  ids = (int **) malloc (6 * sizeof(int *));
  for(k = 0; k < 6; k++)
   ids[k] = (int *) malloc (header->npart[k] * sizeof(int));

  //! CHECK MASSES
  for(n = 0; n < NPartTot; n++)
  {
    header->mass[P[n].Type] += P[n].Mass;
    ref_mass[P[n].Type] = P[n].Mass;

    bob = P[n].Type;

    ids[bob][id_ref[bob]] = n;
    id_ref[bob]++;
  }

  for(k = 0; k < 6; k++)
  {
    if((header->mass[k] - ref_mass[k] * header->npart[k] == 0) && header->npart[k] > 0)
      header->mass[k] = ref_mass[k];
    else
      header->mass[k] = 0;
  }


  //! START WRITING
  fflush(stdout);
  header->num_files = 1;

  dummy = 256;
  fwrite(&dummy,  sizeof(dummy),   1, snap_file);
  fwrite(header, sizeof(gheader), 1, snap_file);
  fwrite(&dummy,  sizeof(dummy),   1, snap_file);


   /*
   for(k = 0; k < 6; k++)
     printf("\tType %i Particles \t%i\n", k, header->npart[k]);
   for(k = 0; k < 6; k++)
     printf("\tType %i Mass      \t%g\n", k, header->mass[k]);
   printf("\tTime                \t%g\n", header->time);
   printf("\tRedshift            \t%g\n", header->redshift);
   printf("\tSFR                 \t%i\n", header->flag_sfr);
   printf("\tFeedback            \t%i\n", header->flag_feedback);
   for(k = 0; k < 6; k++)
     printf("\tType %i npart Total \t%i\n", k, header->npartTotal[k]);
   printf("\tCooling             \t%i\n", header->flag_cooling);
   printf("\tnum files           \t%i\n", header->num_files);
   printf("\tBox Size            \t%g\n", header->BoxSize);
   printf("\tOmega0              \t%g\n", header->Omega0);
   printf("\tOmegaL              \t%g\n", header->OmegaLambda);
   printf("\tHubbleParam         \t%g\n", header->HubbleParam);
   printf("\tFill bytes          \t%s\n", header->fill);
   */

  //!------ Pos
  dummy = 3 * sizeof(P[0].Pos[0]) * NPartTot;
  fwrite(&dummy, sizeof(int), 1, snap_file);
  for (k = 0; k < 6; k++)
    for(n = 0; n < id_ref[k]; n++)
      fwrite(&P[ids[k][n]].Pos[0], sizeof(P[0].Pos[0]), 3, snap_file);
  fwrite(&dummy, sizeof(int), 1, snap_file);

  //!------ Vel
  fwrite(&dummy, sizeof(int), 1, snap_file);
  for(k = 0; k < 6; k++)
    for (n = 0; n < id_ref[k]; n++)
      fwrite(&P[ids[k][n]].Vel[0], sizeof(P[0].Vel[0]), 3, snap_file);
  fwrite(&dummy, sizeof(int), 1, snap_file);

  //!------- ID
  dummy = 4 * NPartTot;
  fwrite(&dummy, sizeof(int), 1, snap_file);
  for(k = 0; k < 6; k++)
    for (n = 0; n < id_ref[k]; n++)
      fwrite(&P[ids[k][n]].Id, sizeof(int), 1, snap_file);
  fwrite(&dummy, sizeof(int), 1, snap_file);

   //!------- Mass
  dummy = 0;
  for(k = 0; k < 6; k++)
    if((header->mass[k] == 0) && (header->npart[k] > 0))
      dummy += header->npart[k];
  dummy *= sizeof(P[0].Mass);

  offset = 0;
  if(dummy != 0)
  {
    //printf("writing mass block\n");
    fwrite(&dummy, sizeof(dummy), 1, snap_file);
    for(k = 0; k < 6; k++)
    {
      if((header->npart[k] > 0) && (header->mass[k] == 0))
      {
        for(n = 0; n < id_ref[k]; n++)
          fwrite(&P[ids[k][n]].Mass, sizeof(P[0].Mass), 1, snap_file);
      }
    }
    fwrite(&dummy, sizeof(dummy), 1, snap_file);
  }

   // printf("done writing\n");

  //! free memory
  for(k = 0; k < 6; k++)
    if(header->npart[k] > 0)
      free(ids[k]);
  free(ids);

  //! Close file
  fclose(snap_file);
}
