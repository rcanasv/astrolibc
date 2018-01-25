/*
 *  \file archive.c
 *  \brief This file contains definition of archive type.
 *         This is meant to be useful for snapshots and
 *         files. Decided to name Archive to avoid confusion
 *         with FILE type.
 *
 */



#include "base.h"
#include "archive.h"


void Archive_fill (Archive * a, char * name, char * format, char * path, char * prefix, int nfiles)
{
  strcpy (a->name,   name);
  strcpy (a->format, format);
  strcpy (a->path,   path);
  strcpy (a->prefix, prefix);
  a->nfiles = nfiles;
}


void Archive_name (Archive * a, char * name)
{
  strcpy (a->name, name);
}


void Archive_format (Archive * a, char * format)
{
  strcpy (a->format, format);
}


void Archive_path (Archive * a, char * path)
{
  strcpy (a->path, path);
}


void Archive_prefix (Archive * a, char * prefix)
{
  strcpy (a->prefix, prefix);
}


void Archive_nfiles (Archive * a, int nfiles)
{
  a->nfiles = nfiles;
}
