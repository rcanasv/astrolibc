/*
 *  \file   hdf5sim.c
 *  \brief  This file contains functions for HDF5.
 *
 */

#include "base.h"
#include "hdf5sim.h"


void hdf5_sim_init_groups (Simulation * sim, HDF5_SimGroup * group)
{
  switch (sim->format)
  {
    case EAGLE:
      strcpy (group->Header,     "Header");
      strcpy (group->GasPart,    "PartType0");
      strcpy (group->DarkPart,   "PartType1");
      strcpy (group->ExtraPart,  "PartType2");
      strcpy (group->TracerPart, "PartType3");
      strcpy (group->StarPart,   "PartType4");
      strcpy (group->BHPart,     "PartType5");
      break;

    case ILLUSTRIS:
      strcpy (group->Header,     "Header");
      strcpy (group->GasPart,    "PartType0");
      strcpy (group->DarkPart,   "PartType1");
      strcpy (group->ExtraPart,  "PartType2");
      strcpy (group->TracerPart, "PartType3");
      strcpy (group->StarPart,   "PartType4");
      strcpy (group->BHPart,     "PartType5");
      break;

    case GIZMO:
      strcpy (group->Header,     "Header");
      strcpy (group->GasPart,    "PartType0");
      strcpy (group->DarkPart,   "PartType1");
      strcpy (group->ExtraPart,  "PartType2");
      strcpy (group->TracerPart, "PartType3");
      strcpy (group->StarPart,   "PartType4");
      strcpy (group->BHPart,     "PartType5");
      break;
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
      break;

    case ILLUSTRIS:
      break;

    case GIZMO:
      break;
  }
}


void hdf5_sim_init (Simulation * snapshot)
{
  int     i;
  FILE  * f;
  char    fname  [NAME_LENGTH];
  char    buffer [NAME_LENGTH];

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
  sprintf (fname, "%s/%s", snapshot->archive.path, snapshot->archive.prefix);
  if ((id_file = H5Fopen (fname, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
  {
    printf ("Couldn't open file %s\n", fname);
    exit (0);
  }

  id_group = H5Gopen (id_file, group.Header, H5P_DEFAULT);
  hdf5_get_attribute (id_group, header.Lbox,          &snapshot->Lbox,                  sizeof(snapshot->Lbox));
  hdf5_get_attribute (id_group, header.HubbleParam,   &snapshot->cosmology.HubbleParam, sizeof(snapshot->cosmology.HubbleParam));
  hdf5_get_attribute (id_group, header.OmegaM,        &snapshot->cosmology.OmegaM,      sizeof(snapshot->cosmology.OmegaM));
  hdf5_get_attribute (id_group, header.OmegaB,        &snapshot->cosmology.OmegaB,      sizeof(snapshot->cosmology.OmegaB));
  hdf5_get_attribute (id_group, header.OmegaL,        &snapshot->cosmology.OmegaL,      sizeof(snapshot->cosmology.OmegaL));
  hdf5_get_attribute (id_group, header.Time,          &snapshot->Time,                  sizeof(snapshot->Time));
  hdf5_get_attribute (id_group, header.z,             &snapshot->z,                     sizeof(snapshot->z));
  hdf5_get_attribute (id_group, header.NpartThisFile, &snapshot->NpartThisFile,         sizeof(snapshot->NpartThisFile[0]));
  hdf5_get_attribute (id_group, header.NpartTot,      &snapshot->NpartTot,              sizeof(snapshot->NpartTot[0]));
  hdf5_get_attribute (id_group, header.MassTable,     &snapshot->MassTable,             sizeof(snapshot->MassTable[0]));
  status = H5Gclose (id_group);
  status = H5Fclose (id_file);

  //
  // Adjust Units
  //

  //
  // Display header values
  //
  /*
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
  */
}


void hdf5_sim_load_particles (Simulation * sim, int filenum, Particle ** part)
{

  int      i, j;
  char     fname  [NAME_LENGTH];
  char     buffer [NAME_LENGTH];

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

  HDF5_SimGroup    group;
  HDF5_SimHeader   header;
  HDF5_PartDset    dataset;

  //
  //  Initialize Depending on  Format
  //
  hdf5_sim_init_groups (snapshot, &group);
  hdf5_sim_init_header (snapshot, &header);

  double * posbuff  = (double *) malloc (3 * sim->NpartThisFile[4] * sizeof(double));
  double * velbuff  = (double *) malloc (3 * sim->NpartThisFile[4] * sizeof(double));
  int    * idbuff   = (int    *) malloc (    sim->NpartThisFile[4] * sizeof(int));
  double * massbuff = (double *) malloc (    sim->NpartThisFile[4] * sizeof(double));

  //
  // Read Header
  //
  sprintf (fname, "%s/%s", snapshot->archive.path, snapshot->archive.prefix);
  if ((id_file = H5Fopen (fname, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
  {
    printf ("Couldn't open file %s\n", fname);
    exit (0);
  }

  id_group = H5Gopen (id_file, group.StarPart, H5P_DEFAULT);
  hdf5_get_data (id_group, dataset.Position,  posbuff,  sizeof(posbuff[0]));
  hdf5_get_data (id_group, dataset.Velocity,  velbuff,  sizeof(velbuff[0]));
  hdf5_get_data (id_group, dataset.Mass,      massbuff, sizeof(massbuff[0]));
  hdf5_get_data (id_group, dataset.ID,        idbuff,   sizeof(idbuff[0]));
  status = H5Gclose (id_group);
  status = H5Fclose (id_file);

  //
  // From buffer to Particle
  //
  if ((*(part) = (Particle *) malloc (sim->NpartThisFile[4] * sizeof(Particle))) == NULL)
  {
    printf ("Couldn't allocate memory for Particle array\n");
    exit(0);
  }

  P = *(part);

  for (i = 0, j = 0; i < sim->NpartThisFile[4]; i++, j=j+3)
  {
    P[i].Pos[0] = posbuff  [j];
    P[i].Pos[1] = posbuff  [j+1];
    P[i].Pos[2] = posbuff  [j+2];
    P[i].Vel[0] = velbuff  [j];
    P[i].Vel[1] = velbuff  [j+1];
    P[i].Vel[2] = velbuff  [j+2];
    P[i].Mass   = massbuff [i];
    P[i].Id     = idbuff   [i];
  }

  free (posbuff);
  free (velbuff);
  free (idbuff);
  free (massbuff);

  for (i = 0; i < sim->NpartThisFile[4]; i++)
  {
    printf ("%f\n", P[i].Pos[0]);
    printf ("%f\n", P[i].Pos[1]);
    printf ("%f\n", P[i].Pos[2]);
  }



  //
  // Convert to human readable units
  //
  /*
  for (i = 0; i < ramses->npart; i++)
  {
    P[i].Pos[0] *= ramses->unit_l;
    P[i].Pos[1] *= ramses->unit_l;
    P[i].Pos[2] *= ramses->unit_l;

    P[i].Vel[0] *= ramses->unit_v;
    P[i].Vel[1] *= ramses->unit_v;
    P[i].Vel[2] *= ramses->unit_v;

    P[i].Mass   *= ramses->unit_m;
  }
  */
}
