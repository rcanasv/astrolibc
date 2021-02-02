/*
 *  \file archive.h
 *  \brief This file contains definition of archive type.
 *         This is meant to be useful for snapshots and
 *         files. Decided to name Archive to avoid confusion
 *         with FILE type.
 *
 */


#ifndef ARCHIVE_H
#define ARCHIVE_H


#include "base.h"
#include "typedef.h"


void Archive_fill   (Archive * a, char * name, char * format, char * path, char * prefix, int nfiles);
void Archive_name   (Archive * a, char * name);
void Archive_format (Archive * a, char * format);
void Archive_path   (Archive * a, char * path);
void Archive_prefix (Archive * a, char * prefix);
void Archive_nfiles (Archive * a, int    nfiles);



#endif    /*  ARCHIVE_H  */
