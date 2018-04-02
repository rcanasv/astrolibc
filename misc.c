/*
 *  \file misc.c
 *  \brief This file contains misc functions.
 *
 */


#include "misc.h"


int get_n_num_from_string (char * strng, int n_num, int ** nums)
{
  int  i, j, k;
  int  ndigits = 10;
  char tmpbuff [ndigits];
  int  length = strlen (strng);

  for (i = 0; i < ndigits; i++)
    tmpbuff[i] = '\0';

  for (i = 0, k = 0; i < length; i++)
    if (strng[i] == ' ')
      k++;

  if (k != n_num)
  {
    printf ("ERROR: there should be %d numbers, but there are only %d\n", n_num, k);
    printf ("Exiting...\n");
  }

  *nums = (int *) malloc (k * sizeof(int));

  for (i = 0, j = 0, k = 0; i < length; i++)
  {
    if (strng[i] != ' ')
    {
      tmpbuff[j] = strng[i];
      tmpbuff[j+1] = '\0';
      j++;
    }
    else
    {
      nums[0][k] = atoi (tmpbuff);
      k++;
      j = 0;
    }
  }
  return 0;
}
