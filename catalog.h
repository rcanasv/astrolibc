/*
 *  \file   catalog.h
 *  \brief  Header file for Catalog structure type
 *
 */

#ifndef CATALOG_H
#define CATALOG_H


#include "base.h"
#include "archive.h"
#include "particle.h"
#include "structure.h"
#include "cosmology.h"
#include "format.h"
#include "simulation.h"



typedef struct Catalog
{
  Archive        archive;
  Cosmology      cosmology;
  int            format;
  int            nstruct;
  int            nprocs;
  int            iprops;
  int            iparts;
  Structure    * strctProps;
  Particle    ** strctParts;
} Catalog;


void Catalog_init                    (Catalog * ctlg);
void Catalog_load                    (Catalog * ctlg);
void Catalog_load_properties         (Catalog * ctlg);
void Catalog_load_particles          (Catalog * ctlg);
void Catalog_free                    (Catalog * ctlg);
void Catalog_fill_SubIDS             (Catalog * ctlg);
void Catalog_get_particle_properties (Catalog * ctlg, Simulation * sim);

//void Catalog_fill_ProgIDs (Catalog * catalog, char * tffile);



#endif    /*  CATALOG_H  */
