/*
 *  \file   hdf5sim.c
 *  \brief  This file contains functions for HDF5.
 *
 */

#include "base.h"
#include "hdf5sim.h"


void hdf5_sim_init_groups (Simulation * sim, HDF5_SimGroup * group)
{

  // Particle types names
  if ((sim->format % 10) < 5)
  {
    // Convention for EAGLE, Illustris,
    // TNG, SIMBA, and MUFASA
    strcpy (group->Header,     "Header");
    strcpy (group->GasPart,    "PartType0");
    strcpy (group->DarkPart,   "PartType1");
    strcpy (group->ExtraPart,  "PartType2");
    strcpy (group->TracerPart, "PartType3");
    strcpy (group->StarPart,   "PartType4");
    strcpy (group->BHPart,     "PartType5");
  }
}


void hdf5_sim_init_header (Simulation * sim, HDF5_SimHeader * header)
{
  switch (sim->format)
  {
    case EAGLE:
      strcpy (header->Lbox,              "BoxSize");
      strcpy (header->Ez,                "E(z)");
      strcpy (header->a,                 "ExpansionFactor");
      strcpy (header->Cooling,           "Flag_Cooling");
      strcpy (header->Double,            "Flag_DoublePrecision");
      strcpy (header->Feedback,          "Flag_Feedback");
      strcpy (header->IcInfo,            "Flag_IC_Info");
      strcpy (header->Metals,            "Flag_Metals");
      strcpy (header->SFR,               "Flag_Sfr");
      strcpy (header->Age,               "Flag_StellarAge");
      strcpy (header->Hz,                "H(z)");
      strcpy (header->HubbleParam,       "HubbleParam");
      strcpy (header->MassTable,         "MassTable");
      strcpy (header->NfilesPerSnapshot, "NumFilesPerSnapshot");
      strcpy (header->NpartThisFile,     "NumPart_ThisFile");
      strcpy (header->NpartTot,          "NumPart_Total");
      strcpy (header->NpartTotHW,        "NumPart_Total_HighWord");
      strcpy (header->OmegaM,            "Omega0");
      strcpy (header->OmegaB,            "OmegaBaryon");
      strcpy (header->OmegaL,            "OmegaLambda");
      strcpy (header->z,                 "Redshift");
      strcpy (header->RunLabel,          "RunLabel");
      strcpy (header->Time,              "Time");
      sim->to_kpc = 1000.0;
      break;

    case ILLUSTRIS:
      break;

    case GIZMO_SIMBA:
      strcpy (header->Lbox,              "BoxSize");
      strcpy (header->MassTable,         "MassTable");
      strcpy (header->NpartThisFile,     "NumPart_ThisFile");
      strcpy (header->NpartTot,          "NumPart_Total");
      strcpy (header->NpartTotHW,        "NumPart_Total_HighWord");
      strcpy (header->OmegaM,            "Omega0");
      strcpy (header->OmegaL,            "OmegaLambda");
      strcpy (header->z,                 "Redshift");
      strcpy (header->Time,              "Time");
      strcpy (header->HubbleParam,       "HubbleParam");
      strcpy (header->NfilesPerSnapshot, "NumFilesPerSnapshot");
      sim->to_kpc = 1.0;
      break;
  }
}


void hdf5_sim_init_dataset (Simulation * sim, HDF5_PartDset * dataset)
{
  strcpy (dataset->Position, "Coordinates");
  strcpy (dataset->Velocity, "Velocity");
  strcpy (dataset->ID,       "ParticleIDs");
  strcpy (dataset->Mass,     "Mass");

  if (sim->format == ILLUSTRIS    ||
      sim->format == TNG          ||
      sim->format == GIZMO_SIMBA)
  {
    strcpy (dataset->Velocity, "Velocities");
    strcpy (dataset->Mass,     "Masses");
  }

  strcpy (dataset->Density, "Density");
  strcpy (dataset->U, "InternalEnergy");
}


void hdf5_sim_init (Simulation * snapshot)
{
  int     i;
  FILE  * f;
  char    fname  [LONG_LENGTH];
  char    buffer [LONG_LENGTH];

  hid_t     id_file;
  hid_t     id_group;

  herr_t    status;

  HDF5_SimGroup    group;
  HDF5_SimHeader   header;

  //
  //  Initialize Depending on  Format
  //
  hdf5_sim_init_groups (snapshot, &group);
  hdf5_sim_init_header (snapshot, &header);

  //
  // Read Header
  //
  sprintf (fname, "%s/%s.hdf5", snapshot->archive.path, snapshot->archive.prefix);
  if ((id_file = H5Fopen (fname, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
  {
    sprintf (fname, "%s/%s.0.hdf5", snapshot->archive.path, snapshot->archive.prefix);
    if ((id_file = H5Fopen (fname, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
    {
      printf ("Couldn't open file %s\n", fname);
      exit (0);
    }
  }
  printf("Opening file  %s\n",fname);


  id_group = H5Gopen (id_file, group.Header, H5P_DEFAULT);
  hdf5_get_attribute (id_group, header.Lbox,          &snapshot->Lbox,                  sizeof(snapshot->Lbox));
  hdf5_get_attribute (id_group, header.HubbleParam,   &snapshot->h,                     sizeof(snapshot->cosmology.HubbleParam));
  hdf5_get_attribute (id_group, header.OmegaM,        &snapshot->cosmology.OmegaM,      sizeof(snapshot->cosmology.OmegaM));
  if (snapshot->format == EAGLE)
    hdf5_get_attribute (id_group, header.OmegaB,        &snapshot->cosmology.OmegaB,      sizeof(snapshot->cosmology.OmegaB));
  hdf5_get_attribute (id_group, header.OmegaL,        &snapshot->cosmology.OmegaL,      sizeof(snapshot->cosmology.OmegaL));
  hdf5_get_attribute (id_group, header.Time,          &snapshot->Time,                  sizeof(snapshot->Time));
  hdf5_get_attribute (id_group, header.z,             &snapshot->z,                     sizeof(snapshot->z));
  hdf5_get_attribute (id_group, header.NpartThisFile, &snapshot->NpartThisFile,         sizeof(snapshot->NpartThisFile[0]));
  hdf5_get_attribute (id_group, header.NpartTot,      &snapshot->NpartTot,              sizeof(snapshot->NpartTot[0]));
  hdf5_get_attribute (id_group, header.MassTable,     &snapshot->MassTable,             sizeof(snapshot->MassTable[0]));


  snapshot->cosmology.HubbleParam = snapshot->h * 100.0;
  snapshot->Lbox = snapshot->Lbox * snapshot->to_kpc / (1.0 + snapshot->z) / snapshot->h;
  snapshot->a = snapshot->Time;

  status = H5Gclose (id_group);
  status = H5Fclose (id_file);


  //
  // Adjust Units
  //

  //
  // Display header values
  //
  printf ("BoxSize          %g\n", snapshot->Lbox);
  printf ("HubbleParam      %g\n", snapshot->cosmology.HubbleParam);
  printf ("Om0              %g\n", snapshot->cosmology.OmegaM);
  printf ("OmB              %g\n", snapshot->cosmology.OmegaB);
  printf ("OmL              %g\n", snapshot->cosmology.OmegaL);
  printf ("time             %g\n", snapshot->Time);
  printf ("z                %g\n", snapshot->z);
  for (i = 0; i < 6; i++)
    printf ("NumPartThisFile  %u\n", snapshot->NpartThisFile[i]);
  for (i = 0; i < 6; i++)
    printf ("NumPartTot       %u\n", snapshot->NpartTot[i]);
  for (i = 0; i < 6; i++)
    printf ("Mass             %g\n", snapshot->MassTable[i]);
}



int hdf5_sim_get_npart_ThisFile (Simulation * sim, int filenum)
{
  int     i;
  FILE  * f;
  char    fname  [LONG_LENGTH];
  char    buffer [LONG_LENGTH];

  hid_t     id_file;
  hid_t     id_group;

  herr_t    status;

  HDF5_SimGroup    group;
  HDF5_SimHeader   header;

  int   nPartThisFile;

  //
  //  Initialize Depending on  Format
  //
  hdf5_sim_init_groups (sim, &group);
  hdf5_sim_init_header (sim, &header);

  //
  // Read Header
  //
  sprintf (fname, "%s/%s.hdf5", sim->archive.path, sim->archive.prefix);
  if ((id_file = H5Fopen (fname, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
  {
    sprintf (fname, "%s/%s.0.hdf5", sim->archive.path, sim->archive.prefix);
    if ((id_file = H5Fopen (fname, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
    {
      printf ("Couldn't open file %s\n", fname);
      exit (0);
    }
  }
  printf("Opening file  %s\n",fname);

  id_group = H5Gopen (id_file, group.Header, H5P_DEFAULT);
  hdf5_get_attribute (id_group, header.NpartThisFile, &sim->NpartThisFile, sizeof(sim->NpartThisFile[0]));

  nPartThisFile = 0;
  for (i = 0; i < 6; i++)
    nPartThisFile += sim->NpartThisFile[i];

  status = H5Gclose (id_group);
  status = H5Fclose (id_file);

  return nPartThisFile;
}


void hdf5_sim_load_particles (Simulation * sim, int filenum, Particle ** part)
{
  int      i, j, k, n;
  char     fname  [LONG_LENGTH];
  char     buffer [LONG_LENGTH];

  int      dummy;
  int      dummyi;
  float    dummyf;
  double   dummyd;
  char     dummys [NAME_LENGTH];

  hid_t    id_file;
  hid_t    id_group;
  hid_t    id_dataset;
  hid_t    id_dataspace;
  hid_t    id_attribute;

  herr_t   status;

  Particle * P;
  double   * posbuff;
  double   * velbuff;
  double   * massbuff;
  long     * idbuff;
  double   * rhobuff;
  double   * ubuff;

  HDF5_SimGroup    group;
  HDF5_SimHeader   header;
  HDF5_PartDset    dataset;

  int  nPartThisFile;
  int  nOffset;

  //
  //  Initialize Depending on  Format
  //
  hdf5_sim_init_groups  (sim, &group);
  hdf5_sim_init_header  (sim, &header);
  hdf5_sim_init_dataset (sim, &dataset);

  //
  // Read Header
  //
  sprintf (fname, "%s/%s.hdf5", sim->archive.path, sim->archive.prefix);
  if ((id_file = H5Fopen (fname, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
  {
    sprintf (fname, "%s/%s.%d.hdf5", sim->archive.path, sim->archive.prefix, filenum);
    if ((id_file = H5Fopen (fname, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
    {
      printf ("Couldn't open file %s\n", fname);
      exit (0);
    }
  }

  //
  // Update local number of particles
  //

  id_group = H5Gopen (id_file, group.Header, H5P_DEFAULT);
  hdf5_get_attribute (id_group, header.NpartThisFile, &sim->NpartThisFile, sizeof(sim->NpartThisFile[0]));
  hdf5_get_attribute (id_group, header.MassTable,     &sim->MassTable,     sizeof(sim->MassTable[0]));
  status = H5Gclose  (id_group);

  nPartThisFile = 0;
  for (k = 0; k < 6; k++)
    nPartThisFile += sim->NpartThisFile[k];

  //
  // From buffer to Particle
  //
  if (((*(part)  = (Particle *) malloc (  nPartThisFile * sizeof(Particle)))  == NULL) || \
      ((posbuff  = (double   *) malloc (3*nPartThisFile * sizeof(double))  )  == NULL) || \
      ((velbuff  = (double   *) malloc (3*nPartThisFile * sizeof(double))  )  == NULL) || \
      ((idbuff   = (long     *) malloc (  nPartThisFile * sizeof(long))    )  == NULL) || \
      ((massbuff = (double   *) malloc (  nPartThisFile * sizeof(double))  )  == NULL)
     )
  {
    printf ("Couldn't allocate memory for Particle array\n");
    exit(0);
  }

  if (sim->NpartThisFile[0])
  {
      rhobuff = (double *) malloc (sim->NpartThisFile[0] * sizeof(double));
      ubuff   = (double *) malloc (sim->NpartThisFile[0] * sizeof(double));
  }

  char gname [6][NAME_LENGTH];
  strcpy (gname[0], group.GasPart);
  strcpy (gname[1], group.DarkPart);
  strcpy (gname[2], group.ExtraPart);
  strcpy (gname[3], group.TracerPart);
  strcpy (gname[4], group.StarPart);
  strcpy (gname[5], group.BHPart);

  nOffset = 0;
  for (k = 0; k < 6; k++)
  {
    if (sim->NpartThisFile[k])
    {
      id_group = H5Gopen (id_file, gname[k], H5P_DEFAULT);
      hdf5_get_data (id_group, dataset.Position,  &posbuff[3*nOffset],  sizeof(posbuff[0]));
      hdf5_get_data (id_group, dataset.Velocity,  &velbuff[3*nOffset],  sizeof(velbuff[0]));
      hdf5_get_data (id_group, dataset.Mass,      &massbuff[nOffset],   sizeof(massbuff[0]));
      hdf5_get_data (id_group, dataset.ID,        &idbuff[nOffset],     sizeof(idbuff[0]));

      if (k == 0 && sim->NpartThisFile[0])
      {
        hdf5_get_data (id_group, dataset.Density, &rhobuff[0],  sizeof(rhobuff[0]));
        hdf5_get_data (id_group, dataset.U,       &ubuff[0],    sizeof(ubuff[0]));
      }

      status = H5Gclose (id_group);
      nOffset += sim->NpartThisFile[k];
    }
  }
  status = H5Fclose (id_file);


  // Now assign values to respective particle properties
  P = *(part);
  for (k = 0, n = 0; k < 6; k++)
  {
    for (i = 0; i < sim->NpartThisFile[k]; i++, n++)
    {
      j = n*3;
      P[n].Pos[0] = posbuff[j]   * sim->to_kpc;
      P[n].Pos[1] = posbuff[j+1] * sim->to_kpc;
      P[n].Pos[2] = posbuff[j+2] * sim->to_kpc;

      P[n].Vel[0] = velbuff[j];
      P[n].Vel[1] = velbuff[j+1];
      P[n].Vel[2] = velbuff[j+2];

      P[n].Mass   = massbuff[n];
      P[n].Id     = idbuff[n];
      P[n].Type   = k;

      if (k == 0)
      {
        P[n].Rho = rhobuff[i];
        P[n].U   = ubuff[i];
      }

      if (i < 10)
        printf ("%e  %e  %e  %ld  %d\n", P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], P[n].Id, P[n].Type);
    }
  }

  // Free memory
  free (posbuff);
  free (velbuff);
  free (idbuff);
  free (massbuff);

  if (sim->NpartThisFile[0])
  {
    free (rhobuff);
    free (ubuff);
  }


  // Convert to human readable units
  double a = sim->a;
  double h = sim->h;

  for (k = 0, n = 0; k < 6; k++)
  {
    for (i = 0; i < sim->NpartThisFile[k]; i++, n++)
    {
     P[n].Pos[0] *= a / h;
     P[n].Pos[1] *= a / h;
     P[n].Pos[2] *= a / h;
     P[n].Mass   *= 1.0 / h;
    }
  }

  //end of function
}
