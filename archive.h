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


typedef struct Archive
{
  FILE * file;
  char   name     [NAME_LENGTH];
  char   format   [NAME_LENGTH];
  char   path     [NAME_LENGTH];
  char   prefix   [NAME_LENGTH];
  int    nFiles;
} Archive;


void Archive_fill   (Archive * a, char * name, char * format, char * path, char * prefix, int nFiles);

void Archive_name   (Archive * a, char * name);

void Archive_format (Archive * a, char * format);

void Archive_path   (Archive * a, char * path);

void Archive_prefix (Archive * a, char * prefix);

void Archive_nFiles (Archive * a, int    nFiles);



#endif    /*  ARCHIVE_H  */
