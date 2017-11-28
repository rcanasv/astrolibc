/*
 *  \file   halomakerio.c
 *  \brief  This file contains Gadget input-output routines
 *
 *
 */

#include "base.h"
#include "halomaker.h"
#include "halomakerio.h"


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
*/

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
