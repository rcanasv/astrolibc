/*
 *  \file gadget.c
 *  \brief This file contains Gadget- format actions
 *
 */

 #include "gadget.h"

// Initalize
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

// Get number of particles in this file
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

// Here we assume that f is properly opened
void gadget_get_header (Simulation * gdt, FILE * f, gheader * header)
{
  int   dummy;
  char  hname [4];

  rewind(f);
  if (gdt->format == GADGET_HEAD)
  {
    fread(&dummy,    sizeof(dummy), 1, f);
    fread(&hname[0], sizeof(hname), 1, f);
    fread(&dummy,    sizeof(dummy), 1, f);
    fread(&dummy,    sizeof(dummy), 1, f);
  }
  fread(&dummy,   sizeof(dummy),   1, f);
  fread(header, sizeof(*header), 1, f);
  fread(&dummy,   sizeof(dummy),   1, f);
}

// Read Property
long gadget_get_property (Simulation * gdt, FILE * f, Particle * P, char * prop)
{
  int   k, pc_new;
  int   dummy;
  char  hname [4];

  if (gdt->format == GADGET_HEAD)
  {
    rewind(f);
    do
    {
      fread(&dummy, sizeof(dummy), 1, f);
      fread(&hname, sizeof(hname), 1, f);
      fread(&dummy, sizeof(dummy), 1, f);
      fread(&dummy, sizeof(dummy), 1, f);

      if (strncmp (hname, prop, 4))
      {
        fread(&dummy, sizeof(dummy), 1, f);
        fseek(f, dummy, SEEK_CUR);
        fread(&dummy, sizeof(dummy), 1, f);
      }
    } while (strncmp (hname, prop, 4));
  }
  return ftell(f);
}

// Get Info
int gadget_get_info (Simulation * gdt, FILE * f, ginfo ** info)
{
  int    k, n, pc_new;
  int    dummy;
  char   hname [4];
  ginfo  * data;

  rewind(f);
  do
  {
    fread(&dummy, sizeof(dummy), 1, f);
    fread(&hname, sizeof(hname), 1, f);
    fread(&dummy, sizeof(dummy), 1, f);
    fread(&dummy, sizeof(dummy), 1, f);

    if (strncmp (hname, "INFO", 4))
    {
      fread(&dummy, sizeof(dummy), 1, f);
      fseek(f, dummy, SEEK_CUR);
      fread(&dummy, sizeof(dummy), 1, f);
    }
  } while (strncmp (hname, "INFO", 4));

  fread(&dummy, sizeof(dummy), 1, f);
  n = dummy / 40;
  *(info) = (ginfo *) malloc (n * sizeof(ginfo));
  data = *(info);
  for (k = 0; k < n; k++)
  {
    fread(&data[k].name, sizeof(int), 1, f); data[k].name[4] = '\0';
    fread(&data[k].type, sizeof(int), 2, f); data[k].type[8] = '\0';
    fread(&data[k].ndim, sizeof(int), 1, f);
    fread(&data[k].flag, sizeof(int), 6, f);
    data[k].read = 0;
  }
  fread(&dummy, sizeof(dummy), 1, f);

  return n;
}

// Load particles
void gadget_load_particles (Simulation * gdt, int filenum, Particle ** part)
{
  FILE  * fd;
  int     n, k, pc, pc_new;
  int     NumPart;
  int     nblocks;
  int     dummy;
  int     ntot_withmasses;
  char    fname [NAME_LENGTH];
  char    hname [4];

  int     i;
  float   ftmp3 [3];
  float   ftmp;

  Particle * P;
  gheader    header1;
  ginfo    * info;

  pc = 0;

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

  // Read INFO
  if (gdt->format == GADGET_HEAD)
    nblocks = gadget_get_info (gdt, fd, &info);

  // Read Header
  rewind(fd);
  gadget_get_header (gdt, fd, &header1);


  // Get npartot in THIS file
  for(k = 0, NumPart = 0; k < 6; k++)
    NumPart += header1.npart[k];

  // Get npart with mass array in THIS file
  for(k = 0, ntot_withmasses = 0; k < 6; k++)
    if(header1.mass[k] == 0)
      ntot_withmasses += header1.npart[k];

  // Allocate memory
  if((*(part) = (Particle *) malloc (NumPart * sizeof(Particle))) == NULL)
  {
    fprintf(stderr, "failed to allocate memory.\n");
    exit(0);
  }
  printf ("Memory allocated\n");
  P = *(part);


  //
  // Set Type of Particles
  //
  memset(P, 0, NumPart*sizeof(Particle));
  for (k = 0, pc_new = pc; k < 6; k++)
    for (n = 0; n < header1.npart[k]; n++)
      P[pc_new++].Type = k;

  //
  // Position
  //
  fseek(fd, gadget_get_property (gdt, fd, P, "POS "), SEEK_SET);
  fread(&dummy,   sizeof(dummy),   1, fd);
  for(k = 0, pc_new = 0; k < 6; k++)
    for(n = 0; n < header1.npart[k]; n++)
    {
      fread(&ftmp3[0], sizeof(float), 3, fd);
      for (i = 0; i < 3; i++)
         P[pc_new].Pos[i] = ftmp3[i] / header1.HubbleParam;
      pc_new++;
    }
  fread(&dummy,   sizeof(dummy),   1, fd);
  for (k = 0; k < nblocks; k++) if (!strncmp(info[k].name, "POS ", 4)) info[k].read = 1;

  //
  // Velocities
  //
  fseek(fd, gadget_get_property (gdt, fd, P, "VEL "), SEEK_SET);
  fread(&dummy,   sizeof(dummy),   1, fd);
  for(k = 0, pc_new = pc; k < 6; k++)
    for(n = 0; n < header1.npart[k]; n++)
    {
      fread(&ftmp3[0], sizeof(float), 3, fd);
      for (i = 0; i < 3; i++)
         P[pc_new].Vel[i] = ftmp3[i];
      pc_new++;
    }
  fread(&dummy,   sizeof(dummy),   1, fd);
  for (k = 0; k < nblocks; k++) if (!strncmp(info[k].name, "VEL ", 4)) info[k].read = 1;

  //
  // IDs
  //
  fseek(fd, gadget_get_property (gdt, fd, P, "ID  "), SEEK_SET);
  fread(&dummy,   sizeof(dummy),   1, fd);  //printf ("%d\n", dummy);
  for(k = 0, pc_new = pc; k < 6; k++)
    for(n = 0; n < header1.npart[k]; n++)
      fread(&P[pc_new++].Id, sizeof(int), 1, fd);
  fread(&dummy,   sizeof(dummy),   1, fd);  //printf ("%d\n", dummy);
  for (k = 0; k < nblocks; k++) if (!strncmp(info[k].name, "ID  ", 4)) info[k].read = 1;

  //
  // Mass
  //
  if(ntot_withmasses>0)
  {
    fseek(fd, gadget_get_property (gdt, fd, P, "MASS"), SEEK_SET);
    fread(&dummy, sizeof(dummy), 1, fd);
    for(k = 0, pc_new = pc; k < 6; k++)
      for(n = 0; n < header1.npart[k]; n++)
      {
        if(header1.mass[k] == 0)
        {
          fread(&ftmp, sizeof(float), 1, fd);
          P[pc_new].Mass = ftmp / header1.HubbleParam;
        }
        else
          P[pc_new].Mass = header1.mass[k] / header1.HubbleParam;
        pc_new++;
      }
    fread(&dummy, sizeof(dummy), 1, fd);
    for (k = 0; k < nblocks; k++) if (!strncmp(info[k].name, "MASS", 4)) info[k].read = 1;
  }



  //
  // U internal energy
  //
  register int iblock;
  for (iblock = 0; iblock < nblocks; iblock++)
  {
    // Internal Energy U
    if (!strncmp(info[iblock].name, "U   ", 4))
    {
      fseek(fd, gadget_get_property (gdt, fd, P, "U   "), SEEK_SET);
      fread(&dummy, sizeof(dummy), 1, fd);
      for(k = 0, pc_new = pc; k < 6; k++)
        if (info[iblock].flag[k])
        {
          for(n = 0; n < header1.npart[k]; n++)
            fread(&P[pc_new++].U, sizeof(float), 1, fd);
        }
        else
          pc_new += header1.npart[k];
      fread(&dummy, sizeof(dummy), 1, fd);
      info[iblock].read = 1;
    }

    // Density Rho
    if (!strncmp(info[iblock].name, "RHO ", 4))
    {
      fseek(fd, gadget_get_property (gdt, fd, P, "RHO "), SEEK_SET);
      fread(&dummy, sizeof(dummy), 1, fd);
      for(k = 0, pc_new = pc; k < 6; k++)
        if (info[iblock].flag[k])
        {
          for(n = 0; n < header1.npart[k]; n++)
            fread(&P[pc_new++].Rho, sizeof(float), 1, fd);
        }
        else
          pc_new += header1.npart[k];
      fread(&dummy, sizeof(dummy), 1, fd);
      info[iblock].read = 1;
    }

    // NE
    if (!strncmp(info[iblock].name, "NE  ", 4))
    {
      fseek(fd, gadget_get_property (gdt, fd, P, "NE  "), SEEK_SET);
      fread(&dummy, sizeof(dummy), 1, fd);
      for(k = 0, pc_new = pc; k < 6; k++)
        if (info[iblock].flag[k])
        {
          for(n = 0; n < header1.npart[k]; n++)
            fread(&P[pc_new++].Rho, sizeof(float), 1, fd);
        }
        else
          pc_new += header1.npart[k];
      fread(&dummy, sizeof(dummy), 1, fd);
      info[iblock].read = 1;
    }

    // NH
    if (!strncmp(info[iblock].name, "NH  ", 4))
    {
      fseek(fd, gadget_get_property (gdt, fd, P, "NH  "), SEEK_SET);
      fread(&dummy, sizeof(dummy), 1, fd);
      for(k = 0, pc_new = pc; k < 6; k++)
        if (info[iblock].flag[k])
        {
          for(n = 0; n < header1.npart[k]; n++)
            fread(&P[pc_new++].Rho, sizeof(float), 1, fd);
        }
        else
          pc_new += header1.npart[k];
      fread(&dummy, sizeof(dummy), 1, fd);
      info[iblock].read = 1;
    }

    // HSML
    if (!strncmp(info[iblock].name, "HSML", 4))
    {
      fseek(fd, gadget_get_property (gdt, fd, P, "HSML"), SEEK_SET);
      fread(&dummy, sizeof(dummy), 1, fd);
      for(k = 0, pc_new = pc; k < 6; k++)
        if (info[iblock].flag[k])
        {
          for(n = 0; n < header1.npart[k]; n++)
            fread(&P[pc_new++].HSML, sizeof(float), 1, fd);
        }
        else
          pc_new += header1.npart[k];
      fread(&dummy, sizeof(dummy), 1, fd);
      info[iblock].read = 1;
    }

    // SFR
    if (!strncmp(info[iblock].name, "SFR ", 4))
    {
      fseek(fd, gadget_get_property (gdt, fd, P, "SFR "), SEEK_SET);
      fread(&dummy, sizeof(dummy), 1, fd);
      for(k = 0, pc_new = pc; k < 6; k++)
        if (info[iblock].flag[k])
        {
          for(n = 0; n < header1.npart[k]; n++)
            fread(&P[pc_new++].SFR, sizeof(float), 1, fd);
        }
        else
          pc_new += header1.npart[k];
      fread(&dummy, sizeof(dummy), 1, fd);
      info[iblock].read = 1;
    }

    // Age
    if (!strncmp(info[iblock].name, "AGE ", 4))
    {
      fseek(fd, gadget_get_property (gdt, fd, P, "AGE "), SEEK_SET);
      fread(&dummy, sizeof(dummy), 1, fd);
      for(k = 0, pc_new = pc; k < 6; k++)
        if (info[iblock].flag[k])
        {
          for(n = 0; n < header1.npart[k]; n++)
            fread(&P[pc_new++].Age, sizeof(float), 1, fd);
        }
        else
          pc_new += header1.npart[k];
      fread(&dummy, sizeof(dummy), 1, fd);
      info[iblock].read = 1;
    }

    // Metallicity Z
    if (!strncmp(info[iblock].name, "Z   ", 4))
    {
      fseek(fd, gadget_get_property (gdt, fd, P, "Z   "), SEEK_SET);
      fread(&dummy, sizeof(dummy), 1, fd);
      for(k = 0, pc_new = pc; k < 6; k++)
        if (info[iblock].flag[k])
        {
          for(n = 0; n < header1.npart[k]; n++)
            fread(&P[pc_new++].Z, sizeof(float), 1, fd);
        }
        else
          pc_new += header1.npart[k];
      fread(&dummy, sizeof(dummy), 1, fd);
      info[iblock].read = 1;
    }

    // Density NH
    if (!strncmp(info[iblock].name, "TEMP", 4))
    {
      fseek(fd, gadget_get_property (gdt, fd, P, "TEMP"), SEEK_SET);
      fread(&dummy, sizeof(dummy), 1, fd);
      for(k = 0, pc_new = pc; k < 6; k++)
        if (info[iblock].flag[k])
        {
          for(n = 0; n < header1.npart[k]; n++)
            fread(&P[pc_new++].T, sizeof(float), 1, fd);
        }
        else
          pc_new += header1.npart[k];
      fread(&dummy, sizeof(dummy), 1, fd);
      info[iblock].read = 1;
    }
  }

  for (k = 0; k < nblocks; k++)
    printf ("%d  %s   read  %d\n", k, info[k].name, info[k].read);

  // Close File
  printf("Successful snapshot reading\n");
  fclose(fd);

  return;
}

// Write gadget file
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
  header->flag_cooling    = 0;
  header->time            = 1;
  header->redshift        = 0;
  header->flag_sfr        = 0;
  header->flag_feedback   = 0;
  header->BoxSize         = 1000000;
  header->Omega0          = 0.28;
  header->OmegaLambda     = 0.72;
  header->HubbleParam     = 0.677;


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

   for(k = 0; k < 6; k++)
     printf("\ttype  %i idref \t%i\n", k, id_ref[k]);

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
  dummy = sizeof(P[0].Id) * NPartTot;
  fwrite(&dummy, sizeof(int), 1, snap_file);
  for(k = 0; k < 6; k++)
    for (n = 0; n < id_ref[k]; n++)
      fwrite(&P[ids[k][n]].Id, sizeof(P[0].Id), 1, snap_file);
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

  //!------- U
  if (id_ref[0])
  {
    dummy = sizeof(P[0].U) * id_ref[0];
    fwrite(&dummy, sizeof(int), 1, snap_file);
    for (n = 0; n < id_ref[0]; n++)
      fwrite(&P[ids[0][n]].Id, sizeof(P[0].U), 1, snap_file);
    fwrite(&dummy, sizeof(int), 1, snap_file);

    //!------- RHO
    dummy = sizeof(P[0].Rho) * id_ref[0];
    fwrite(&dummy, sizeof(int), 1, snap_file);
    for (n = 0; n < id_ref[0]; n++)
      fwrite(&P[ids[0][n]].Rho, sizeof(P[0].Rho), 1, snap_file);
    fwrite(&dummy, sizeof(int), 1, snap_file);

    //!------- T
    dummy = sizeof(P[0].T) * id_ref[0];
    fwrite(&dummy, sizeof(int), 1, snap_file);
    for (n = 0; n < id_ref[0]; n++)
      fwrite(&P[ids[0][n]].T, sizeof(P[0].T), 1, snap_file);
    fwrite(&dummy, sizeof(int), 1, snap_file);
  }
  printf("done writing\n");

  //! free memory
  for(k = 0; k < 6; k++)
    if(header->npart[k] > 0)
      free(ids[k]);
  free(ids);

  //! Close file
  fclose(snap_file);
}


/*
for (k = 0; k < nblocks; k++)
{
  printf ("%s\n", info[k].name);
  printf ("%s\n", info[k].type);
  printf ("ndim  %d\n", info[k].ndim);
  for (i = 0; i < 6; i++)
    printf ("  ptype  %d\n", info[k].flag[i]);
}

for(k = 0; k < 6; k++)
  printf("\tType %i Particles \t%i\n", k, header1.npart[k]);
for(k = 0; k < 6; k++)
  printf("\tType %i Mass      \t%g\n", k, header1.mass[k]);
printf("\tTime                \t%g\n", header1.time);
printf("\tRedshift            \t%g\n", header1.redshift);
printf("\tSFR                 \t%i\n", header1.flag_sfr);
printf("\tFeedback            \t%i\n", header1.flag_feedback);
for(k = 0; k < 6; k++)
  printf("\tType %i npart Total \t%i\n", k, header1.npartTotal[k]);
printf("\tCooling             \t%i\n", header1.flag_cooling);
printf("\tnum files           \t%i\n", header1.num_files);
printf("\tBox Size            \t%g\n", header1.BoxSize);
printf("\tOmega0              \t%g\n", header1.Omega0);
printf("\tOmegaL              \t%g\n", header1.OmegaLambda);
printf("\tHubbleParam         \t%g\n", header1.HubbleParam);
printf("\tFill bytes          \t%s\n", header1.fill);

*/


/*
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
 */
