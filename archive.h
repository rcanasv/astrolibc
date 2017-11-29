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
  char   format   [NAME_LENGTH];
  char   path     [NAME_LENGTH];
  char   prefix   [NAME_LENGTH];
  int    nFiles;
} Archive;


#endif    /*  ARCHIVE_H  */
