/*
 *  \file   gadgetio.c
 *  \brief  This file contains Gadget input-output routines
 *
 *
 */

#include "base.h"
#include "gadget.h"
#include "gadgetio.h"


//
//
//---------- Read Gadget Snapshot ---------//
//
//
void read_gadget_snapshot(char * snapshot)
{
 FILE *fd;
 int k, dummy, ntot_withmasses;
 int n, pc, pc_new;

 pc = 0;

 if((fd = fopen(snapshot, "r")) == NULL)
 {
   printf("Couldn't open file\n");
   exit(0);
 }

 fflush(stdout);

 fread(&dummy, sizeof(dummy), 1, fd);
 fread(&header1, sizeof(header1), 1, fd);
 fread(&dummy, sizeof(dummy), 1, fd);

 for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
   NumPart += header1.npart[k];

 for(k = 0, ntot_withmasses = 0; k < 6; k++)
   if(header1.mass[k] == 0)
     ntot_withmasses += header1.npart[k];

 printf("Allocating memory...");

 if(!(p = malloc(NumPart * sizeof(pdata_s))))
   {
     fprintf(stderr, "failed to allocate memory.\n");
     exit(0);
   }

 if(!(ID = malloc(NumPart * sizeof(int))))
   {
     fprintf(stderr, "failed to allocate memory.\n");
     exit(0);
   }

 printf("done\n");

 SKKIP;

 for(k = 0, pc_new = pc; k < 6; k++)
   for(n = 0; n < header1.npart[k]; n++)
   {
     fread(&p[pc_new].Pos[0], sizeof(float), 3, fd);
     pc_new++;
   }

 SKKIP;
 SKKIP;

 for(k = 0, pc_new = pc; k < 6; k++)
   for(n = 0; n < header1.npart[k]; n++)
   {
     fread(&p[pc_new].Vel[0], sizeof(float), 3, fd);
     pc_new++;
   }

 SKKIP;
 SKKIP;

 for(k = 0, pc_new = pc; k < 6; k++)
   for(n = 0; n < header1.npart[k]; n++)
   {
     fread(&p[pc_new].Id, sizeof(int), 1, fd);
     ID[pc_new] = p[pc_new].Id;
     pc_new++;
   }
 SKKIP;

 if(ntot_withmasses>0)
   SKKIP;
 for(k = 0, pc_new = pc; k < 6; k++)
 {
   for(n=0; n < header1.npart[k]; n++)
   {
     p[pc_new].Type = k;
     if(header1.mass[k] == 0)
	fread(&p[pc_new].Mass, sizeof(float), 1, fd);
     else
	p[pc_new].Mass= header1.mass[k];
     pc_new++;
   }
 }
 if(ntot_withmasses>0)
   SKKIP;


 int check;
 for(k = 0, pc_new = pc; k < 6; k++)
 {
   check = 0;
   for(n = 0; n < header1.npart[k]-1; n++)
   {
     if(p[pc_new + 1].Mass == p[pc_new].Mass)
	check++;
     pc_new++;
   }
   if(header1.npart[k] > 0) pc_new++;
   if((check == header1.npart[k]-1) && (header1.mass[k] == 0))
     header1.mass[k] = p[pc_new-1].Mass;
 }
//   printf("NumPart =   %i\n", NumPart);
//   printf("PCNEW =     %i\n", pc_new);
 printf("Successful snapshot reading\n");

 fclose(fd);
}



 //
 //
 //---------- Write Snapshot File ----------//
 //
 //
 void write_snapshot(struct pdata_s * pg, int particles, gheader header, char * output_name)
 {
   FILE * snap_file;

   int dummy;
   int k, pc_new, n;
   int pc = 0;
   int ntot_withmasses;

   double ref_mass[6];
   int id_ref[6];
   int ** ids;

   int i, bob, offset;

   if((snap_file = fopen(output_name,"w")) == NULL)
   {
     printf("Couldn't open file\n");
     exit(0);
   }

   //! Initialize
   for(k = 0; k < 6; k++)
   {
     header.npart[k]      = 0;
     header.npartTotal[k] = 0;
     header.mass[k]       = 0;
     ref_mass[k]          = 0;
     id_ref[k]            = 0;
   }
   header.time           = 0.9823335728;
   header.redshift       = 0.0179841427;
   header.flag_sfr       = 0;
   header.flag_feedback  = 0;
   header.flag_cooling   = 0;
   header.num_files      = 1;
   header.BoxSize        = 100000;
   header.Omega0         = 0.2720000148;
   header.OmegaLambda    = 0.7279999852;
   header.HubbleParam    = 0.7040000153;


   //! Count Particles by type
   for(k = 0; k < particles; k++)
   {
     header.npart[pg[k].Type]++;
     header.npartTotal[pg[k].Type]++;
   }

   ids = (int **)malloc(6 * sizeof(int *));
   for(k = 0; k < 6; k++)
    ids[k] = (int *)malloc((header.npart[k]) * sizeof(int));

   //! CHECK MASSES
   for(n = 0; n < particles; n++)
   {
     header.mass[pg[n].Type] += pg[n].Mass;
     ref_mass[pg[n].Type] = pg[n].Mass;

     bob = pg[n].Type;

     ids[bob][id_ref[bob]] = n;
     id_ref[bob]++;
   }

   for(k = 0; k < 6; k++)
   {
     if((header.mass[k] - ref_mass[k] * header.npart[k] == 0) && header.npart[k] > 0)
       header.mass[k] = ref_mass[k];
     else
       header.mass[k] = 0;
 //     printf("%i\n", k);
 //     printf("      Npart  %i\n", header.npart[k]);
 //     printf("      Mass   %g\n", header.mass[k]);
   }

 //   offset = 0;
 //   for(k = 0; k < 6; k++)
 //   {
 //     if(header.npart[k] > 0)
 //       for(i = 0; i < header.npart[k]; i++)
 // 	ID[offset + i] = ids[k][i];
 //
 //     offset += header.npart[k];
 //   }

   //! START WRITING
   fflush(stdout);
   header.num_files = 1;

   dummy = 256;
   fwrite(&dummy, sizeof(dummy), 1, snap_file);
   fwrite(&header, sizeof(header), 1, snap_file);
   fwrite(&dummy, sizeof(dummy), 1, snap_file);

 //   for(k = 0; k < 6; k++)
 //     printf("\tType %i Particles \t%i\n", k, header.npart[k]);
 //   for(k = 0; k < 6; k++)
 //     printf("\tType %i Mass      \t%g\n", k, header.mass[k]);
 //   printf("\tTime                \t%g\n", header.time);
 //   printf("\tRedshift            \t%g\n", header.redshift);
 //   printf("\tSFR                 \t%i\n", header.flag_sfr);
 //   printf("\tFeedback            \t%i\n", header.flag_feedback);
 //   for(k = 0; k < 6; k++)
 //     printf("\tType %i npart Total \t%i\n", k, header.npartTotal[k]);
 //   printf("\tCooling             \t%i\n", header.flag_cooling);
 //   printf("\tnum files           \t%i\n", header.num_files);
 //   printf("\tBox Size            \t%g\n", header.BoxSize);
 //   printf("\tOmega0              \t%g\n", header.Omega0);
 //   printf("\tOmegaL              \t%g\n", header.OmegaLambda);
 //   printf("\tHubbleParam         \t%g\n", header.HubbleParam);
 //   printf("\tFill bytes          \t%s\n", header.fill);
 //

   //!------ Pos
   dummy = 3 * 4 * particles;
   fwrite(&dummy, sizeof(dummy), 1, snap_file);
   for(k = 0; k < particles; k++)
     fwrite(&pg[k].Pos[0], sizeof(float), 3, snap_file);
   fwrite(&dummy, sizeof(dummy), 1, snap_file);

   //!------ Vel
   fwrite(&dummy, sizeof(dummy), 1, snap_file);
   for(k = 0; k < particles; k++)
     fwrite(&pg[k].Vel[0], sizeof(float), 3, snap_file);
   fwrite(&dummy, sizeof(dummy), 1, snap_file);

   //!------- ID
   dummy = 4 * particles;
   fwrite(&dummy, sizeof(dummy), 1, snap_file);
   for(k = 0; k < particles; k++)
     fwrite(&pg[k].Id, sizeof(int), 1, snap_file);
   fwrite(&dummy, sizeof(dummy), 1, snap_file);

    //!------- Mass
   dummy = 0;
   for(k = 0; k < 6; k++)
     if((header.mass[k] == 0) && (header.npart[k] > 0))
       dummy += header.npart[k];
   dummy *= sizeof(float);

   offset = 0;
   if(dummy != 0)
   {
 //     printf("writing mass block\n");
     fwrite(&dummy, sizeof(dummy), 1, snap_file);
     for(k = 0; k < 6; k++)
     {
       if((header.npart[k] > 0) && (header.mass[k] == 0))
       {
 	for(i = 0; i < header.npart[k]; i++)
 	  fwrite(&pg[k].Mass, sizeof(float), 1, snap_file);
       }
     }
     fwrite(&dummy, sizeof(dummy), 1, snap_file);
   }

 //   printf("done writing\n");

   //! free memory
   for(k = 0; k < 6; k++)
     if(header.npart[k] > 0)
       free(ids[k]);
   free(ids);

   //! Close file
   fclose(snap_file);
 }
