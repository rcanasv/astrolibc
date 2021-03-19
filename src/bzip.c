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
  fclose(f);

  printf ("%s\n", buf);
  printf ("filesize  %d\n", fsize);



  /*
  // Read the whole file
  if (nbytes == 0)
  {
    bzerror = BZ_OK;
    fsize = 0;
    while (bzerror == BZ_OK)
    {
      nBuf = BZ2_bzRead (&bzerror, b, buf, 1000);
      fsize += nBuf;
      if (bzerror == BZ_OK)
      {
        //printf ("%s\n", buf);
        //exit(0);
      }
    }
    if (bzerror == BZ_STREAM_END)
      printf ("End of stream\n");
    BZ2_bzReadClose (&bzerror, b);

    rewind(f);
    *(fbuff) = (char *) malloc (fsize);
    char * ptr = *(fbuff);
    b = BZ2_bzReadOpen (&bzerror, f, 0, 0, NULL, 0);
    nBuf = BZ2_bzRead (&bzerror, b, ptr, fsize);
    printf ("%s\n", &ptr[fsize-500]);

    if (bzerror == BZ_STREAM_END)
      printf ("End of stream\n");
    printf ("BZ_MAX_UNUSED  %d\n", BZ_MAX_UNUSED);

    void * unused;
    int    nUnused;
    BZ2_bzReadGetUnused (&bzerror, b, &unused, &nUnused);
    printf ("bzerror %d  nunused  %d\n", bzerror, nUnused);
    b = BZ2_bzReadOpen (&bzerror, f, 0, 0, unused, nUnused);
    nBuf = BZ2_bzRead (&bzerror, b, buf, 1000);
    printf ("%s\n", buf);
    //int i;
    //for (i = 0; i < 15; i++)
    //  printf ("%c", ptr[i]);
    //printf ("%c", ptr[fsize-500+i]);
    //for (i = 0; i < 500; i++)
    //printf ("%s\n", &ptr[fsize-500]);
    //printf ("\n");
    BZ2_bzReadClose (&bzerror, b);
    fclose(f);
  }
  exit(0);
  */




  /*
  if (bzerror != BZ_STREAM_END)
  {
    BZ2_bzReadClose (&bzerror, b);
    // HANDlE ERROR
  }
  else
  {
    BZ2_bzReadClose (&bzerror, b);
  }*/

}
