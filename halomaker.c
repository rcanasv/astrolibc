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


  hmkr->a                 = aexp;
  hmkr->cosmology.OmegaM  = omegam;
  hmkr->AgeUniv           = ageuniv;

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
    hmkr->strctProps[i].SubIDs     = NULL;
    hmkr->strctProps[i].MatchIDs   = NULL;
    hmkr->strctProps[i].MatchMrrts = NULL;
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

    hmkr->strctProps[i].dummyi = 0;

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
    double Lbox   = 100000 * hmkr->a / 0.704;
    double Lbox_2 = Lbox / 2.0;

    HMKR_SKIP    fread (&dummyf, sizeof(float), 1, f);
                 hmkr->strctProps[i].Pos[0] = dummyf * 1000 + Lbox_2;
                 fread (&dummyf, sizeof(float), 1, f);
                 hmkr->strctProps[i].Pos[1] = dummyf * 1000 + Lbox_2;
                 fread (&dummyf, sizeof(float), 1, f);
                 hmkr->strctProps[i].Pos[2] = dummyf * 1000 + Lbox_2;
                                                        HMKR_SKIP

    /*
    while (hmkr->strctProps[i].Pos[0] < 0) hmkr->strctProps[i].Pos[0] + Lbox;
    while (hmkr->strctProps[i].Pos[1] < 0) hmkr->strctProps[i].Pos[1] + Lbox;
    while (hmkr->strctProps[i].Pos[2] < 0) hmkr->strctProps[i].Pos[2] + Lbox;

    while (hmkr->strctProps[i].Pos[0] > Lbox) hmkr->strctProps[i].Pos[0] - Lbox;
    while (hmkr->strctProps[i].Pos[1] > Lbox) hmkr->strctProps[i].Pos[1] - Lbox;
    while (hmkr->strctProps[i].Pos[2] > Lbox) hmkr->strctProps[i].Pos[2] - Lbox;
    */

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


void halomaker_read_particles (Catalog * hmkr)
{
  int     i, j, k;

  int     dummy;
  int     dummyi;
  long    dummyl;
  float   dummyf;
  double  dummyd;

  FILE  * f;
  char    treebricks_fname [NAME_LENGTH];

  //
  // Open treebricks file to read properties
  //
  sprintf (treebricks_fname, "%s/%s", hmkr->archive.path, hmkr->archive.prefix);
  if ((f = fopen (treebricks_fname, "r")) == NULL)
  {
    printf ("ERROR: Cannot open file  %s\n", treebricks_fname);
    exit (0);
  }

  HMKR_SKIP    fseek (f, dummy, SEEK_CUR);    HMKR_SKIP
  HMKR_SKIP    fseek (f, dummy, SEEK_CUR);    HMKR_SKIP
  HMKR_SKIP    fseek (f, dummy, SEEK_CUR);    HMKR_SKIP
  HMKR_SKIP    fseek (f, dummy, SEEK_CUR);    HMKR_SKIP
  HMKR_SKIP    fseek (f, dummy, SEEK_CUR);    HMKR_SKIP
  HMKR_SKIP    fseek (f, dummy, SEEK_CUR);    HMKR_SKIP

  //
  // Loop over structures
  //
  for (i = 1; i <= hmkr->nstruct; i++)
  {
    // Number of particles
    HMKR_SKIP    fseek (f, dummy, SEEK_CUR);    HMKR_SKIP

    // List of particles in galaxy i
    HMKR_SKIP
      hmkr->strctProps[i].PIDs = (int *) malloc (hmkr->strctProps[i].NumPart * sizeof(int));
      fread (hmkr->strctProps[i].PIDs,  sizeof(int)*hmkr->strctProps[i].NumPart, 1, f);
    HMKR_SKIP

    // Structure ID
    HMKR_SKIP    fseek (f, dummy, SEEK_CUR);    HMKR_SKIP

    // TimeStep Number
    HMKR_SKIP    fseek (f, dummy, SEEK_CUR);    HMKR_SKIP

    // Level, Parent ID, ID first sub, Nsubstructures, ID next sub
    HMKR_SKIP    fseek (f, dummy, SEEK_CUR);    HMKR_SKIP

    // Mass
    HMKR_SKIP    fseek (f, dummy, SEEK_CUR);    HMKR_SKIP

    // Galaxy Position
    HMKR_SKIP    fseek (f, dummy, SEEK_CUR);    HMKR_SKIP

    // Galaxy Velocity
    HMKR_SKIP    fseek (f, dummy, SEEK_CUR);    HMKR_SKIP

    // Galaxy Angular Momentum
    HMKR_SKIP    fseek (f, dummy, SEEK_CUR);    HMKR_SKIP

    // Distance to most distant particle
    HMKR_SKIP    fseek (f, dummy, SEEK_CUR);    HMKR_SKIP

    // Kinetic Potential and Total Energy
    HMKR_SKIP    fseek (f, dummy, SEEK_CUR);    HMKR_SKIP

    // Spin Parameter
    HMKR_SKIP    fseek (f, dummy, SEEK_CUR);    HMKR_SKIP

    // Velocity dispersion, bulge velocity disp, bulge mass
    HMKR_SKIP    fseek (f, dummy, SEEK_CUR);    HMKR_SKIP

    // Virial Radius, Virial Mass, Virial Temp, Sound Speed
    HMKR_SKIP    fseek (f, dummy, SEEK_CUR);    HMKR_SKIP

    // Central Density (NFW), Characteristic Radius (NFW)
    HMKR_SKIP    fseek (f, dummy, SEEK_CUR);    HMKR_SKIP

    // Stellar surface density profiles
    HMKR_SKIP    fseek (f, dummy, SEEK_CUR);    HMKR_SKIP  // nbins
    HMKR_SKIP    fseek (f, dummy, SEEK_CUR);    HMKR_SKIP  // radius
    HMKR_SKIP    fseek (f, dummy, SEEK_CUR);    HMKR_SKIP  // Stellar surface density
  }

  fclose(f);
}




// ---  Read Galfile --- //
void halomaker_read_galfile (Simulation * sim, Structure * strct)
{
  int        dummy;
  double     dummyd;
  int        i, j;

  FILE     * f;
  Particle * P;
  char       fname[NAME_LENGTH];

  int        gal_number;
  int        gal_level;
  double     gal_mass;
  double     gal_pos[3];
  double     gal_vel[3];
  double     gal_ang[3];
  int        nlist;


  sprintf (fname, "%s/gal_stars_%07d", sim->archive.name, strct->ID);
  if ((f = fopen (fname, "r")) == NULL)
  {
    printf ("Can't open file named   %s \n", fname);
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
  if (strct->NumPart != nlist)
  {
    printf ("Error number of particles doesn't match\n");
    printf ("Exiting\n");
    exit (0);
  }

  if (!(strct->Part = malloc (nlist * sizeof(Particle))))
  {
    printf ("Cannot allocate memory for particle information\n");
    printf ("Exiting\n");
    exit (0);
  }
  else
    strct->iPart = 1;

  P = strct->Part;

  // Positions
  for (j = 0; j < 3; j++)
  {
    HMKR_SKIP
    if (dummy == 8 * nlist)
    {
      for (i = 0; i < nlist; i++)
      {
        fread (&dummyd, sizeof(double), 1, f);
        P[i].Pos[j] = dummyd;
      }
    }
    else
    {
      printf ("Number of Particles and Positions does not coincide\n");
      return;
    }
    HMKR_SKIP
  }

  // Velocities
  for (j = 0; j < 3; j++)
  {
    HMKR_SKIP
    if (dummy == 8 * nlist)
    {
      for (i = 0; i < nlist; i++)
      {
        fread (&dummyd, sizeof(double), 1, f);
        P[i].Vel[j] = dummyd;
      }
    }
    else
    {
      printf ("Number of Particles and Velocities does not coincide\n");
      return;
    }
    HMKR_SKIP
  }

  // Masses
  HMKR_SKIP
  if (dummy == 8 * nlist)
  {
    for (i = 0; i < nlist; i++)
    {
      fread (&dummyd, sizeof(double), 1, f);
      P[i].Mass = dummyd;
    }
  }
  else
  {
    printf ("Number of Particles and Masses does not coincide\n");
    return;
  }
  HMKR_SKIP

  // Ids
  HMKR_SKIP
  if (dummy == 4 * nlist)
  {
    for (i = 0; i < nlist; i++)
      fread (&P[i].Id, sizeof(int), 1, f);
  }
  else
  {
    printf ("Number of Particles and Ids does not coincide\n");
    return;
  }
  HMKR_SKIP

  // Age
  HMKR_SKIP
  if (dummy == 8 * nlist)
  {
    for (i = 0; i < nlist; i++)
    {
      fread (&dummyd, sizeof(double), 1, f);
      P[i].Age = dummyd;
    }
  }
  else
  {
    printf ("Number of Particles and Birth Epoch does not coincide\n");
    return;
  }
  HMKR_SKIP

  // Metallicity
  HMKR_SKIP
  if (dummy == 8 * nlist)
  {
    for (i = 0; i < nlist; i++)
    {
      fread (&dummyd, sizeof(double), 1, f);
      P[i].Metal = dummyd;
    }
  }
  else
  {
    printf ("Number of Particles and Metallicity does not coincide\n");
    return;
  }
  HMKR_SKIP

  // Chemical Elements
  HMKR_SKIP
  fseek (f, dummy, SEEK_CUR);
  /*
  if (dummy == 8 * nlist)
  {
    for (i = 0; i < nlist; i++)
      fread (&gal->Part[i].Chem, sizeof(double), 1, f);
  }
  else
  {
    printf ("Number of Particles and Chemical Elements does not coincide\n");
    return;
  }
  */
  HMKR_SKIP

  fclose(f);

  //
  // Unit conversion is done here
  //
  for (i = 0; i < nlist; i++)
  {
    // Shift to centre of mass position
    P[i].Pos[0] -= gal_pos[0];
    P[i].Pos[1] -= gal_pos[1];
    P[i].Pos[2] -= gal_pos[2];

    // Shift to centre of mass velocity
    P[i].Vel[0] -= gal_vel[0];
    P[i].Vel[1] -= gal_vel[1];
    P[i].Vel[2] -= gal_vel[2];

    // convert from Mpc to kpc
    P[i].Pos[0] *= 1000;
    P[i].Pos[1] *= 1000;
    P[i].Pos[2] *= 1000;
    
    // Move galaxy to box  0 - Lbox
    P[i].Pos[0] += strct->Pos[0];
    P[i].Pos[1] += strct->Pos[1];
    P[i].Pos[2] += strct->Pos[2];
    
    // converts to Msun
    P[i].Mass   *= 1e+11;
  }
}



void halomaker_catalog_get_particle_properties (Catalog * hmkr, Simulation * sim)
{
  int         i, j;
  char        fname[NAME_LENGTH];
  Structure * strct;

  if (sim->format == GALFILE)
  {
    for (i = 1; i <= hmkr->nstruct; i++)
    {
      strct = &hmkr->strctProps[i];
      halomaker_read_galfile (sim, strct);
    }
  }
  else
  {
    // Would need to load all files and then distribute particles
    // according to their id to the corresponding structure
    ;
  }
}


void halomaker_structure_get_particle_properties (Catalog * hmkr, Simulation * sim, int * strct_to_get)
{
  int i;
  Structure * strct;
  if (sim->format == GALFILE)
  {
    for (i = 0; i < hmkr->nstruct; i++)
      if (strct_to_get[i])
      {
        strct = &hmkr->strctProps[i];
        halomaker_read_galfile (sim, strct);
      }
  }
  /*
  Else would need to open ramses file and look for IDs....
  */
}
