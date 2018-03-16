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
      strcpy (groups->Header,     "Header");
      strcpy (groups->GasPart,    "PartType0");
      strcpy (groups->DarkPart,   "PartType1");
      strcpy (groups->ExtraPart,  "PartType2");
      strcpy (groups->TracerPart, "PartType3");
      strcpy (groups->StarPart,   "PartType4");
      strcpy (groups->BHPart,     "PartType5");
      break;

    case ILLUSRIS:
      strcpy (groups->Header,     "Header");
      strcpy (groups->GasPart,    "PartType0");
      strcpy (groups->DarkPart,   "PartType1");
      strcpy (groups->ExtraPart,  "PartType2");
      strcpy (groups->TracerPart, "PartType3");
      strcpy (groups->StarPart,   "PartType4");
      strcpy (groups->BHPart,     "PartType5");
      break;

    case GIZMO:
      strcpy (groups->Header,     "Header");
      strcpy (groups->GasPart,    "PartType0");
      strcpy (groups->DarkPart,   "PartType1");
      strcpy (groups->ExtraPart,  "PartType2");
      strcpy (groups->TracerPart, "PartType3");
      strcpy (groups->StarPart,   "PartType4");
      strcpy (groups->BHPart,     "PartType5");
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
      strcpy (header->NpartTot           "NumPart_Total");
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
  herr_t    id_status;

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
  sprintf (fname, "%s/info_%s.txt", snapshot->archive.path, snapshot->archive.prefix);
  if ((id_file = H5Fopen (argv[1], H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
  {
    printf ("Couldn't open file %s\n", fname);
    exit (0);
  }

  id_group = H5Gopen (id_file, group.Header, H5P_DEFAULT);
  get_attribute (id_group, header->Lbox,          &snapshot->Lbox,                  sizeof(snapshot->Lbox));
  get_attribute (id_group, header->HubbleParam,   &snapshot->cosmology.HubbleParam, sizeof(snapshot->cosmology.HubbleParam));
  get_attribute (id_group, header->OmegaM,        &snapshot->cosmology.OmegaM,      sizeof(snapshot->cosmology.OmegaM));
  get_attribute (id_group, header->OmegaB,        &snapshot->cosmology.OmegaB,      sizeof(snapshot->cosmology.OmegaB));
  get_attribute (id_group, header->OmegaL,        &snapshot->cosmology.OmegaL,      sizeof(snapshot->cosmology.OmegaL));
  get_attribute (id_group, header->Time,          &snapshot->Time,                  sizeof(snapshot->Time));
  get_attribute (id_group, header->z,             &snapshot->z,                     sizeof(snapshot->z));
  get_attribute (id_group, header->NpartThisFile, &snapshot->NpartThisFile,         sizeof(snapshot->NpartThisFile[0]));
  get_attribute (id_group, header->NpartTot,      &snapshot->NpartTot,              sizeof(snapshot->NpartTot[0]));
  get_attribute (id_group, header->MassTable,     &snapshot->MassTable,             sizeof(snapshot->MassTable[0]));
  status = H5Gclose (id_group);
  status = H5Fclose (id_file);

  //
  // Adjust Units
  //

  //
  // Display header values
  //
  printf ("BoxSize          %g\n", snapshot->Lbox);
  printf ("HubbleParam      %g\n", snapshot->HubbleParam);
  printf ("Om0              %g\n", snapshot->OmegaM);
  printf ("OmB              %g\n", snapshot->OmegaB);
  printf ("OmL              %g\n", snapshot->OmegaL);
  printf ("time             %g\n", snapshot->Time);
  printf ("z                %g\n", snapshot->z);
  for (i = 0; i < 6; i++)
    printf ("NumPartThisFile  %u\n", snapshot->NpartThisFile[i]);
  for (i = 0; i < 6; i++)
    printf ("NumPartTot       %u\n", snapshot->NpartTot[i]);
  for (i = 0; i < 6; i++)
    printf ("Mass             %g\n", snapshot->MassTable[i]);
}
