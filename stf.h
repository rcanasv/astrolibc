/*
 *  \file stf.h
 *  \brief  This file contains all VELOCIraptor-stf related stuff
 *
 */

#ifndef STF_H
#define STF_H


#include "allvars.h"
#include "base.h"
#include "particle.h"
#include "structure.h"
#include "archive.h"
#include "catalog.h"
#include "misc.h"


typedef struct stfExtendedOutput
{
  int  oIndex;
  int  IdStruct;
  int  IdHost;
  int  IdIGM;
} stfExtendedOutput;


void  stf_read_properties                   (Catalog * stf);
void  stf_read_treefrog                     (Archive * tfrog, Catalog * stf);

void  stf_write_catalog_group               (Catalog * stf);
void  stf_write_catalog_particles           (Catalog * stf);

void  stf_get_particle_properties           (Catalog * stf, Simulation * sim);
void  stf_catalog_get_particle_properties   (Catalog * stf, Simulation * sim);
void  stf_structure_get_particle_properties (Catalog * stf, Simulation * sim, int * strcts_to_get);

int   stf_load_extended_output              (Catalog * stf, int filenum, stfExtendedOutput ** xtndd);
int   stf_get_files_to_read                 (Catalog * stf, int * strcts_to_get, int * files_to_read);

void  stf_catalog_fill_isolated             (Catalog * stf);

//int load_stf_extended_output (char * prefix, int filenum);
//void free_extended_arrays (void);
//int read_stf_filesofgroup (char * prefix, int strct_id, int ** files_of_strct);


#endif    /*  STF_H  */
