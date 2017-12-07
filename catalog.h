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

typedef struct Catalog
{
  Archive        archive;
  Cosmology      cosmology;
  int            nstruct;
  int            nprocs;
  int            iprops;
  int            iparts;
  Structure    * strctProps;
  Particle    ** strctParts;
} Catalog;


void Catalog_init (Catalog * catalog);

void Catalog_load (Catalog * catalog);

void Catalog_free (Catalog * catalog);

void Catalog_fill_SubIDS (Catalog * catalog);

//void Catalog_fill_ProgIDs (Catalog * catalog, char * tffile);



#endif    /*  CATALOG_H  */
