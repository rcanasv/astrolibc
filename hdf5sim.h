/*
 *  \file   hdf5.h
 *  \brief  Header file for HDF5
 *
 */


#ifndef HDF5_H
#define HDF5_H


#define HDF5GAS     0
#define HDF5DM      1
#define HDF5EXTRA   2
#define HDF5TRACER  3
#define HDF5STAR    4
#define HDF5BH      5


typedef HDF5_SimHeader
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
} HDF5_SimGroup;


typedef HDF5_SimGroup
{
  char Header     [NAME_LENGTH];
  char GasPart    [NAME_LENGTH];
  char DarkPart   [NAME_LENGTH];
  char ExtraPart  [NAME_LENGTH];
  char StarPart   [NAME_LENGTH];
  char BHPart     [NAME_LENGTH];
} HDF5_SimGroup;


void  hdf5_sim_init          (Simulation * sim);
void  hdf5_sim_init_groups   (Simulation * sim, HDF5_SimGroup * group);
void  hdf5_sim_init_header   (Simulation * sim, HDF5_SimHeader * header);


//void  hdf5_sim_load_particles                (Simulation * sim, int filenum, Particle ** part);
//void  hdf5_sim_structure_calculate_star_age  (Simulation * sim, Structure * strct);
//void  hdf5_sim_catalog_calculate_star_age    (Simulation * sim, Catalog * ctlg);


#endif    /*  HDF5_H  */
