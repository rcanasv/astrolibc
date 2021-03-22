/*
 *  \file   bzip.c
 *  \brief  This file contains functions for bzip2
 *
 */

#include "base.h"
#include "bzip.h"


void read_bz2_file (char * fname, char ** fbuff, int nbytes)
{
  int      fsize;
  int      i;
  char   * ptr;

  // bzip2  variables
  FILE   * f;
  BZFILE * b;
  int      nBuf;
  char     buf [BZ_BUF_SIZE];
  int      bzerror;
  void   * unused = NULL;
  int      nUnused = 0;

  // Open file
  f = fopen (fname, "r");
  if (!f)
  {
    printf ("Couldn't open %s  bz2\n", fname);
    exit (0);
  }

  if (nbytes)
  {
    bzerror = BZ_OK;
    b = BZ2_bzReadOpen (&bzerror, f, 0, 0, unused, nUnused);
    if (bzerror != BZ_OK)
    {
      BZ2_bzReadClose (&bzerror, b);
      printf ("Couldn't read bz2 file\n");
      exit (0);
    }

    *(fbuff) = (char *) malloc (nbytes);
    ptr = *(fbuff);
    nBuf = BZ2_bzRead (&bzerror, b, ptr, nbytes);
    BZ2_bzReadClose (&bzerror, b);
  }
  else
  {
    // Prepar bzip Buffer
    bzerror = BZ_OK;
    fsize = 0;
    do
    {
      if (bzerror == BZ_STREAM_END)
        BZ2_bzReadGetUnused (&bzerror, b, &unused, &nUnused);

      b = BZ2_bzReadOpen (&bzerror, f, 0, 0, unused, nUnused);
      if (bzerror != BZ_OK)
      {
        BZ2_bzReadClose (&bzerror, b);
        printf ("Couldn't read bz2 file\n");
        exit (0);
      }

      while (bzerror == BZ_OK)
      {
        memset (buf, 0, BZ_BUF_SIZE);
        nBuf = BZ2_bzRead (&bzerror, b, buf, BZ_BUF_SIZE);
        fsize += nBuf;
      }
    } while (!feof(f));
    BZ2_bzReadClose (&bzerror, b);


    // Allocate memory buffer and info into it
    *(fbuff) = (char *) malloc (fsize);
    ptr = *(fbuff);

    // Read file again
    rewind (f);
    unused  = NULL;
    nUnused = 0;
    bzerror = BZ_OK;
    do
    {
      if (bzerror == BZ_STREAM_END)
        BZ2_bzReadGetUnused (&bzerror, b, &unused, &nUnused);

      b = BZ2_bzReadOpen (&bzerror, f, 0, 0, unused, nUnused);
      if (bzerror != BZ_OK)
      {
        BZ2_bzReadClose (&bzerror, b);
        printf ("Couldn't read bz2 file\n");
        exit (0);
      }

      while (bzerror == BZ_OK)
      {
        nBuf = BZ2_bzRead (&bzerror, b, ptr, BZ_BUF_SIZE);
        ptr += nBuf;
      }
    } while (!feof(f));
    BZ2_bzReadClose (&bzerror, b);
  }
  fclose (f);
}
