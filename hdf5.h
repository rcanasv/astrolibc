/*
 *  \file   hdf5.h
 *  \brief  Header file for HDF5
 *
 */


#ifndef HDF5_H
#define HDF5_H


herr_t finfo              (hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata);
void   hdf5_get_attribute (hid_t obj_id, char * attr_name, void * buffer, size_t vsize);
void   hdf5_get_data      (hid_t obj_id, char * dset_name, void * buffer, size_t vsize);


#endif    /*  HDF5_H  */
