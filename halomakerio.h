/*
 *  \file   halomakerio.h
 *  \brief  This header file contains HaloMaker input-output stuff
 *              -  Treebricks
 *              -  gal_XXXX_XXXX files
 *
 *
 *
 */

#ifndef HALOMAKERIO_H
#define HALOMAKERIO_H

#include "base.h"
#include "particle.h"
#include "halomaker.h"


#define HMKR_SKIP  fread  (&dummy, sizeof(int), 1, f);
#define HMKR_PSKIP printf ("%d \n", dummy);


//
//  Prototypes
//

// Read Gadget snapshot
void read_gadget_snapshot(char * snapshot);

// Write Gadget snapshot
void write_snapshot(struct pdata_s * pg, int particles, gheader header, char * output_name);




#endif    /*  GADGETIO_H */
