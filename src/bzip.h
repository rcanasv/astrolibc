/*
 *  \file   bzip.h
 *  \brief  Bzip file handling
 *
 */

#ifndef BZIP_H
#define BZIP_H

#include <bzlib.h>

#define BZ_BUF_SIZE  1000

void read_bz2_file (char * fname, char ** fbuff, int nbytes);

#endif /* BZIP_H_ */
