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


herr_t finfo              (hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata);
void   hdf5_get_attribute (hid_t obj_id, char * attr_name, void * buffer, size_t vsize);
void   hdf5_get_data      (hid_t obj_id, char * dset_name, void * buffer, size_t vsize);

void  hdf5_sim_init                          (Simulation * sim);
void  hdf5_sim_load_particles                (Simulation * sim, int filenum, Particle ** part);
//void  hdf5_sim_structure_calculate_star_age  (Simulation * sim, Structure * strct);
//void  hdf5_sim_catalog_calculate_star_age    (Simulation * sim, Catalog * ctlg);

#endif    /*  HDF5_H  */
