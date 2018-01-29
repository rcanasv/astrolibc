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

/*
typedef struct stfOutput
{
  Archive        stf;
  int            nstruct;
  int            nprocs;
  int            iprops;
  int            iparts;
  objProps     * strctProps;
  pdata_s     ** strctParts;
} stfOutput;
*/

typedef struct stfExtendedOutput
{
  int  oIndex;
  int  IdStruct;
  int  IdHost;
  int  IdIGM;
} stfExtendedOutput;


void  stf_read_properties                   (Catalog * stf);
void  stf_write_catalog_group               (Catalog * stf);
void  stf_write_catalog_particles           (Catalog * stf);
void  stf_read_treefrog                     (Archive * tfrog, Catalog * stf);
void  stf_get_particle_properties           (Catalog * stf, Simulation * sim);
void  stf_catalog_get_particle_properties   (Catalog * stf, Simulation * sim);
void  stf_structure_get_particle_properties (Catalog * stf, int id, Simulation * sim);
int   stf_load_extended_output              (Catalog * stf, int filenum, stfExtendedOutput ** xtndd);

//int load_stf_extended_output (char * prefix, int filenum);
//void free_extended_arrays (void);
//int read_stf_filesofgroup (char * prefix, int strct_id, int ** files_of_strct);


#endif    /*  STF_H  */
