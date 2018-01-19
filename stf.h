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


int * extended_oIndex;
int * extended_IdStruct;
int * extended_IdHost;
int * extended_IdIGM;

void stf_read_properties (Catalog * stf);

void stf_write_catalog_group (Catalog * stf);

void stf_write_catalog_particles (Catalog * stf);

//int load_stf_extended_output (char * prefix, int filenum);
//void free_extended_arrays (void);
//int read_stf_filesofgroup (char * prefix, int strct_id, int ** files_of_strct);
//int load_treefrog (char * tffile, int strct_id, int ** prog_ids, float ** prog_mrrts);


#endif    /*  STF_H  */
