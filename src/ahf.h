/*
 *  \file stf.h
 *  \brief  This file contains all VELOCIraptor-stf related stuff
 *
 */

#ifndef AHF_H
#define AHF_H


#include "allvars.h"
#include "base.h"
#include "particle.h"
#include "structure.h"
#include "archive.h"
#include "catalog.h"
#include "misc.h"
#include "bzip.h"


void  ahf_read_properties           (Catalog * ahf);
void  ahf_catalog_get_particle_list (Catalog * ahf);

//int load_stf_extended_output (char * prefix, int filenum);
//void free_extended_arrays (void);
//int read_stf_filesofgroup (char * prefix, int strct_id, int ** files_of_strct);


#endif    /*  AHF_H  */
