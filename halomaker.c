/*
 *  \file   halomaker.c
 *  \brief  This file contains Gadget input-output routines
 *
 *
 */

#include "base.h"
#include "particle.h"
#include "structure.h"
#include "catalog.h"
#include "halomaker.h"


//---------- Read Propertes File ----------//

void halomaker_read_properties (Catalog * hmkr)
{
  int     i, j, k;

  int     dummy;
  int     dummyi;
  long    dummyl;
  float   dummyf;
  double  dummyd;

  int     nparts;
  float   massres;
  float   aexp;
  float   omegam;
  float   ageuniv;
  int     nstruct;
  int     nsubs;

  FILE  * f;
  char    treebricks_fname [NAME_LENGTH];

  int     iReadParticles = 0;
  int     iGetProfile    = 0;

  //
  // Open treebricks file to read properties
  //
  sprintf (treebricks_fname, "%s/%s", hmkr->archive.path, hmkr->archive.prefix);
  if ((f = fopen (treebricks_fname, "r")) == NULL)
  {
    printf ("ERROR: Cannot open file  %s\n", treebricks_fname);
    exit (0);
  }

  HMKR_SKIP    fread (&nparts,   sizeof(int),   1, f);   HMKR_SKIP
  HMKR_SKIP    fread (&massres,  sizeof(float), 1, f);   HMKR_SKIP
  HMKR_SKIP    fread (&aexp,     sizeof(float), 1, f);   HMKR_SKIP
  HMKR_SKIP    fread (&omegam,   sizeof(float), 1, f);   HMKR_SKIP
  HMKR_SKIP    fread (&ageuniv,  sizeof(float), 1, f);   HMKR_SKIP
  HMKR_SKIP    fread (&nstruct,  sizeof(int),   1, f);
               fread (&nsubs,    sizeof(int),   1, f);   HMKR_SKIP


  hmkr->cosmology.Aexp    = aexp;
  hmkr->cosmology.OmegaM  = omegam;
  hmkr->cosmology.AgeUniv = ageuniv;

  hmkr->nprocs  = 1;
  hmkr->nstruct = nstruct + nsubs;


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
  // Loop over structures
  //
  for (i = 1; i <= hmkr->nstruct; i++)
  {
    // Number of particles
    HMKR_SKIP    fread (&dummyi, sizeof(int), 1, f);    HMKR_SKIP
    hmkr->strctProps[i].NumPart = dummyi;

    // List of particles in galaxy i
    HMKR_SKIP
    if (iReadParticles)
    {
      // hmo->strctProps[i].partIDs = (int *) malloc (hmo->strctProps[i].nparts * sizeof(int));
      // fread (hmo->strctProps[i].partIDs,  sizeof(int), hmo->strctProps[i].nparts, f);
    }
    else
      fseek (f, dummy, SEEK_CUR);
    HMKR_SKIP

    // Structure ID
    HMKR_SKIP    fread (&dummyi, sizeof(int), 1, f);    HMKR_SKIP
    hmkr->strctProps[i].ID = dummyi;

    // TimeStep Number
    HMKR_SKIP    fread (&dummyi, sizeof(int), 1, f);    HMKR_SKIP
    hmkr->strctProps[i].Timestep = dummyi;

    // Level, Parent ID, ID first sub, Nsubstructures, ID next sub
    HMKR_SKIP    fread (&dummyi, sizeof(int), 1, f);
                 hmkr->strctProps[i].Level = dummyi;
                 fread (&dummyi, sizeof(int), 1, f);
                 hmkr->strctProps[i].DirectHostID = dummyi;
                 fread (&dummyi, sizeof(int), 1, f);
                 hmkr->strctProps[i].IDfirstSub = dummyi;
                 fread (&dummyi, sizeof(int), 1, f);
                 hmkr->strctProps[i].NumSubs = dummyi;
                 fread (&dummyi, sizeof(int), 1, f);
                 hmkr->strctProps[i].IDnextSub = dummyi;
                                                        HMKR_SKIP

    if (hmkr->strctProps[i].DirectHostID == hmkr->strctProps[i].ID)
    {
      hmkr->strctProps[i].Type         = 10;
      hmkr->strctProps[i].HostID       = -1;
      hmkr->strctProps[i].DirectHostID = -1;
    }
    else
    {
      hmkr->strctProps[i].Type   = 20;
      hmkr->strctProps[i].HostID = hmkr->strctProps[i].DirectHostID;
    }

    // Mass
    HMKR_SKIP    fread (&dummyf, sizeof(float), 1, f);  HMKR_SKIP
    hmkr->strctProps[i].TotMass = dummyf * 1e+11;       // Mass is now in Solar Masses

    // Galaxy Position
    HMKR_SKIP    fread (&dummyf, sizeof(float), 1, f);
                 hmkr->strctProps[i].Pos[0] = dummyf * 1000;
                 fread (&dummyf, sizeof(float), 1, f);
                 hmkr->strctProps[i].Pos[1] = dummyf * 1000;
                 fread (&dummyf, sizeof(float), 1, f);
                 hmkr->strctProps[i].Pos[2] = dummyf * 1000;
                                                        HMKR_SKIP
    // Position is now in -Lbox/2 to Lbox/2 in kpc
    double Lbox = 100000 * hmkr->cosmology.Aexp / 0.704;

    // Galaxy Velocity
    HMKR_SKIP    fread (&dummyf, sizeof(float), 1, f);
                 hmkr->strctProps[i].Vel[0] = dummyf;
                 fread (&dummyf, sizeof(float), 1, f);
                 hmkr->strctProps[i].Vel[1] = dummyf;
                 fread (&dummyf, sizeof(float), 1, f);
                 hmkr->strctProps[i].Vel[2] = dummyf;   HMKR_SKIP

    // Galaxy Angular Momentum
    HMKR_SKIP    fread (&dummyf, sizeof(float), 1, f);
                 hmkr->strctProps[i].L[0] = dummyf;
                 fread (&dummyf, sizeof(float), 1, f);
                 hmkr->strctProps[i].L[1] = dummyf;
                 fread (&dummyf, sizeof(float), 1, f);
                 hmkr->strctProps[i].L[2] = dummyf;     HMKR_SKIP

    // Distance to most distant particle
    HMKR_SKIP    fread (&dummyf, sizeof(float), 1, f);
                 hmkr->strctProps[i].Rsize = dummyf;

    //  Inertia tensor eigvals
                 fread (&dummyf, sizeof(float), 1, f);
                 hmkr->strctProps[i].Ix[0] = dummyf;
                 fread (&dummyf, sizeof(float), 1, f);
                 hmkr->strctProps[i].Ix[1] = dummyf;
                 fread (&dummyf, sizeof(float), 1, f);
                 hmkr->strctProps[i].Ix[2] = dummyf;    HMKR_SKIP

    // Kinetic Potential and Total Energy
    HMKR_SKIP    fread (&dummyf, sizeof(float), 1, f);
                 hmkr->strctProps[i].Ekin = dummyf;
                 fread (&dummyf, sizeof(float), 1, f);
                 hmkr->strctProps[i].Epot = dummyf;
                 fread (&dummyf, sizeof(float), 1, f);
                 hmkr->strctProps[i].Etot = dummyf;     HMKR_SKIP

    // Spin Parameter
    HMKR_SKIP    fread (&dummyf, sizeof(float), 1, f);  HMKR_SKIP
    hmkr->strctProps[i].Lambda = dummyf;

    // Velocity dispersion, bulge velocity disp, bulge mass
    HMKR_SKIP    fread (&dummyf, sizeof(float), 1, f);
                 hmkr->strctProps[i].Vdisp = dummyf;
                 fread (&dummyf, sizeof(float), 1, f);
                 fread (&dummyf, sizeof(float), 1, f);  HMKR_SKIP

   // Virial Radius, Virial Mass, Virial Temp, Sound Speed
   HMKR_SKIP    fread (&dummyf, sizeof(float), 1, f);
                hmkr->strctProps[i].Rvir = dummyf;
                fread (&dummyf, sizeof(float), 1, f);
                hmkr->strctProps[i].Mvir = dummyf;
                fread (&dummyf, sizeof(float), 1, f);
                hmkr->strctProps[i].Tvir = dummyf;
                fread (&dummyf, sizeof(float), 1, f);
                hmkr->strctProps[i].Csvir = dummyf;     HMKR_SKIP

   // Central Density (NFW), Characteristic Radius (NFW)
   HMKR_SKIP    fread (&dummyf, sizeof(float), 1, f);
                hmkr->strctProps[i].RhoNFW = dummyf;
                fread (&dummyf, sizeof(float), 1, f);
                hmkr->strctProps[i].ReNFW = dummyf;     HMKR_SKIP

   // Stellar surface density profiles
   HMKR_SKIP    fread (&dummyi, sizeof(int), 1, f);     HMKR_SKIP

   HMKR_SKIP
   if (iGetProfile)
   {
     // hmo->strctProps[i].rbin = (float *) malloc (hmo->strctProps[i].nbins * sizeof(float));
     // fread (hmo->strctProps[i].rbin, sizeof(float), hmo->strctProps[i].nbins, f);
   }
   else
     fseek (f, dummy, SEEK_CUR);
   HMKR_SKIP

   HMKR_SKIP
   if (iGetProfile)
   {
     // hmo->strctProps[i].sbin = (float *) malloc (hmo->strctProps[i].nbins * sizeof(float));
     // fread (hmo->strctProps[i].sbin, sizeof(float), hmo->strctProps[i].nbins, f);
   }
   else
     fseek (f, dummy, SEEK_CUR);
   HMKR_SKIP
 }

  fclose(f);
}




// ---  Read Galfile --- //

/*
void halomaker_read_galfile (Structure * strct)
{
  int     dummy;
  int     i, j;
  FILE *  f;

  int      gal_number;
  int      gal_level;
  double   gal_mass;
  double   gal_pos[3];
  double   gal_vel[3];
  double   gal_ang[3];
  int      nlist;

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
