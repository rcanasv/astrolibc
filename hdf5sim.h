/*
 *  \file   hdf5.h
 *  \brief  Header file for HDF5
 *
 */


#ifndef HDF5SIM_H
#define HDF5SIM_H


#include "base.h"
#include "typedef.h"
#include "hdf5routines.h"
#include "format.h"


typedef struct HDF5_SimHeader
{
  char Lbox               [NAME_LENGTH];
  char Ez                 [NAME_LENGTH];
  char a                  [NAME_LENGTH];
  char Cooling            [NAME_LENGTH];
  char Double             [NAME_LENGTH];
  char Feedback           [NAME_LENGTH];
  char IcInfo             [NAME_LENGTH];
  char Metals             [NAME_LENGTH];
  char SFR                [NAME_LENGTH];
  char Age                [NAME_LENGTH];
  char Hz                 [NAME_LENGTH];
  char HubbleParam        [NAME_LENGTH];
  char MassTable          [NAME_LENGTH];
  char NfilesPerSnapshot  [NAME_LENGTH];
  char NpartThisFile      [NAME_LENGTH];
  char NpartTot           [NAME_LENGTH];
  char NpartTotHW         [NAME_LENGTH];
  char OmegaM             [NAME_LENGTH];
  char OmegaB             [NAME_LENGTH];
  char OmegaL             [NAME_LENGTH];
  char z                  [NAME_LENGTH];
  char RunLabel           [NAME_LENGTH];
  char Time               [NAME_LENGTH];
} HDF5_SimHeader;


typedef struct HDF5_SimGroup
{
  char Header      [NAME_LENGTH];
  char GasPart     [NAME_LENGTH];
  char DarkPart    [NAME_LENGTH];
  char ExtraPart   [NAME_LENGTH];
  char TracerPart  [NAME_LENGTH];
  char StarPart    [NAME_LENGTH];
  char BHPart      [NAME_LENGTH];
} HDF5_SimGroup;


typedef struct HDF5_PartDset
{
  char Position    [NAME_LENGTH];
  char Velocity    [NAME_LENGTH];
  char Mass        [NAME_LENGTH];
  char ID          [NAME_LENGTH];
} HDF5_PartDset;


void  hdf5_sim_init           (Simulation * sim);
void  hdf5_sim_init_groups    (Simulation * sim, HDF5_SimGroup  * group);
void  hdf5_sim_init_header    (Simulation * sim, HDF5_SimHeader * header);
void  hdf5_sim_load_particles (Simulation * sim, int filenum, Particle ** part);

//void  hdf5_sim_load_particles                (Simulation * sim, int filenum, Particle ** part);
//void  hdf5_sim_structure_calculate_star_age  (Simulation * sim, Structure * strct);
//void  hdf5_sim_catalog_calculate_star_age    (Simulation * sim, Catalog * ctlg);


#endif    /*  HDF5_H  */
