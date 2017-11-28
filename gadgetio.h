/*
 *  \file   gadgetio.h
 *  \brief  This header file contains Gadget input-output stuff
 *
 *
 *
 */

#ifndef GADGETIO_H
#define GADGETIO_H

#include "base.h"
#include "particle.h"
#include "gadget.h"


#define GSKIP fread(&dummy, sizeof(dummy), 1, fd);


//
//  Prototypes
//

// Read Gadget snapshot
void read_gadget_snapshot(char * snapshot);

// Write Gadget snapshot
void write_snapshot(pdata_s * pg, int particles, gheader header, char * output_name);



#endif    /*  GADGETIO_H */
