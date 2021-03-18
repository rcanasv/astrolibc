/*
 *  \file misc.h
 *  \brief
 *
 *
 */

 #ifndef MISC_H
 #define MISC_H


#include "allvars.h"
#include "base.h"
#include "particle.h"
#include "structure.h"
#include "archive.h"
#include "catalog.h"


int get_n_num_from_string (char * strng, int n_num, int ** nums);

int long_compare (const void * a, const void * b);

#endif    /*  MISC_H  */
