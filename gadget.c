/*
 *  \file gadget.c
 *  \brief This file contains Gadget- format actions
 *
 */

 #include "gadget.h"


void gadget_load_particles (Simulation * gdt, int filenum, Particle ** part)
{
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
    ref_mass[k]          = 0;
    id_ref[k]            = 0;
  }
  header->time           = 1;
  header->redshift       = 0;
  header->flag_sfr       = 0;
  header->flag_feedback  = 0;
  header->flag_cooling   = 0;
  header->num_files      = 1;
  header->BoxSize        = 100000;
  header->Omega0         = 1;
  header->OmegaLambda    = 1;
  header->HubbleParam    = 1;


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
//     printf("%i\n", k);
//     printf("      Npart  %i\n", header->npart[k]);
//     printf("      Mass   %g\n", header->mass[k]);
  }

//   offset = 0;
//   for(k = 0; k < 6; k++)
//   {
//     if(header->npart[k] > 0)
//       for(i = 0; i < header->npart[k]; i++)
// 	ID[offset + i] = ids[k][i];
//
//     offset += header->npart[k];
//   }

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
  dummy = 3 * sizeof(float) * NPartTot;
  fwrite(&dummy, sizeof(int), 1, snap_file);
  for(k = 0; k < NPartTot; k++)
    fwrite(&P[k].Pos[0], sizeof(float), 3, snap_file);
  fwrite(&dummy, sizeof(int), 1, snap_file);

  //!------ Vel
  fwrite(&dummy, sizeof(int), 1, snap_file);
  for(k = 0; k < NPartTot; k++)
    fwrite(&P[k].Vel[0], sizeof(float), 3, snap_file);
  fwrite(&dummy, sizeof(int), 1, snap_file);

  //!------- ID
  dummy = 4 * NPartTot;
  fwrite(&dummy, sizeof(int), 1, snap_file);
  for(k = 0; k < NPartTot; k++)
    fwrite(&P[k].Id, sizeof(int), 1, snap_file);
  fwrite(&dummy, sizeof(int), 1, snap_file);

   //!------- Mass
  dummy = 0;
  for(k = 0; k < 6; k++)
    if((header->mass[k] == 0) && (header->npart[k] > 0))
      dummy += header->npart[k];
  dummy *= sizeof(float);

  offset = 0;
  if(dummy != 0)
  {
//     printf("writing mass block\n");
    fwrite(&dummy, sizeof(dummy), 1, snap_file);
    for(k = 0; k < 6; k++)
    {
      if((header->npart[k] > 0) && (header->mass[k] == 0))
      {
        for(i = 0; i < header->npart[k]; i++)
          fwrite(&P[i].Mass, sizeof(float), 1, snap_file);
      }
    }
    fwrite(&dummy, sizeof(dummy), 1, snap_file);
  }

//   printf("done writing\n");

  //! free memory
  for(k = 0; k < 6; k++)
    if(header->npart[k] > 0)
      free(ids[k]);
  free(ids);

  //! Close file
  fclose(snap_file);
}
