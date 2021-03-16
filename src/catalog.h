/*
 *  \file   catalog.h
 *  \brief  Header file for Catalog structure type
 *
 */

#ifndef CATALOG_H
#define CATALOG_H


#include "base.h"
#include "typedef.h"
#include "archive.h"
#include "cosmology.h"
#include "structure.h"
#include "simulation.h"
#include "particle.h"
#include "format.h"

#include "simulation.h"
#include "stf.h"
#include "halomaker.h"
#include "ahf.h"


void Catalog_init (Catalog * ctlg);
void Catalog_free (Catalog * ctlg);

void Catalog_load            (Catalog * ctlg);
void Catalog_load_properties (Catalog * ctlg);
void Catalog_load_particles  (Catalog * ctlg);

void Catalog_fill_isolated   (Catalog * ctlg);
void Catalog_fill_SubIDS     (Catalog * ctlg);

void Catalog_get_particle_properties (Catalog * ctlg, Simulation * sim);
void Catalog_get_files_of_groups     (Catalog * ctlg);
void Catalog_get_particle_list       (Catalog * ctlg);

//void Catalog_fill_ProgIDs (Catalog * catalog, char * tffile);



#endif    /*  CATALOG_H  */
