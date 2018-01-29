/*
 *  \file   catalog.c
 *  \brief  This file contains functions for Catalog struct type.
 *
 *
 *
 */

#include "catalog.h"
#include "stf.h"
#include "halomaker.h"



void Catalog_init (Catalog * catalog)
{
  catalog->nprocs     = 0;
  catalog->nstruct    = 0;
  catalog->iprops     = 0;
  catalog->iparts     = 0;
  catalog->strctProps = NULL;
  catalog->strctParts = NULL;

  if ((!strcmp(catalog->archive.format, "stf"))          ||  \
      (!strcmp(catalog->archive.format, "VELOCIraptor")) ||  \
      (!strcmp(catalog->archive.format, "velociraptor")))
    catalog->format = STF;
  else
    if ((!strcmp(catalog->archive.format, "hmkr"))            ||  \
        (!strcmp(catalog->archive.format, "hmkr_treebricks")) ||  \
        (!strcmp(catalog->archive.format, "HaloMaker"))       ||  \
        (!strcmp(catalog->archive.format, "halomaker")))
      catalog->format = HALOMAKER;
    else
    {
      printf ("Format %s not supported\n", catalog->archive.format);
      printf ("Exiting...\n");
      exit (0);
    }
}



void Catalog_load (Catalog * catalog)
{
  Catalog_load_properties (catalog);
  Catalog_load_particles  (catalog);
}



void Catalog_load_properties (Catalog * ctlg)
{
  if (ctlg->format == STF)         stf_read_properties (ctlg);
  if (ctlg->format == HALOMAKER)   halomaker_read_properties (ctlg);
}



void Catalog_load_particles (Catalog * ctlg)
{
  if (!ctlg->iprops)
    Catalog_load_properties (ctlg);

  //if (catalog->format == STF)         stf_read_particles (catalog);
  if (ctlg->format == HALOMAKER)   halomaker_read_particles (catalog);
}



void Catalog_get_particle_properties (Catalog * ctlg, Archive * arx)
{
  if (ctlg->format == STF)        stf_catalog_get_particle_properties (ctlg, arx);
  if (ctlg->format == HALOMAKER)  halomaker_catalog_get_particle_properties (ctlg, arx);
}



void Catalog_free (Catalog * catalog)
{
  int i;

  if (catalog->iprops)
  {
    for (i = 1; i <= catalog->nstruct; i++)
    {
      if (catalog->strctProps[i].SubIDs != NULL)
        free (catalog->strctProps[i].SubIDs);

      if (catalog->strctProps[i].MatchIDs != NULL)
        free (catalog->strctProps[i].MatchIDs);

      if (catalog->strctProps[i].MatchMrrts != NULL)
        free (catalog->strctProps[i].MatchMrrts);
    }
    free (catalog->strctProps);
  }


  if (catalog->iparts)
  {
    for (i = 1; i <= catalog->nstruct; i++)
      free (catalog->strctParts[i]);
    free (catalog->strctParts);
  }
}



void Catalog_fill_SubIDS (Catalog * catalog)
{
  int i, j;
  int bob;
  int tmp;


  for (i = 1; i <= catalog->nstruct; i++)
    if (catalog->strctProps[i].HostID == -1)
    {
      bob = catalog->strctProps[i].NumSubs;
      catalog->strctProps[i].SubIDs = (int *) malloc (bob * sizeof(int));
      catalog->strctProps[i].dummy = 0;
      for (j = 0; j < bob; j++)
        catalog->strctProps[i].SubIDs[j] = 0;
    }

  for (i = 1; i <= catalog->nstruct; i++)
  {
    if (catalog->strctProps[i].HostID != -1)
    {
      bob = i;
      while (catalog->strctProps[bob].HostID != -1)
      {
        bob = catalog->strctProps[bob].DirectHostID;
      }
      tmp = catalog->strctProps[bob].dummy++;
      catalog->strctProps[bob].SubIDs[tmp] = i;
    }
  }
}


/*
void Catalog_fill_ProgIDs (Catalog * catalog, char * tffile)
{
  int i, j;
  int bob;
  int tmpid ;
  int nprogs;

  FILE * f = fopen (tffile, "r");

  fgets (buffer, LONG_LENGTH, f);
  fgets (buffer, LONG_LENGTH, f);
  fgets (buffer, LONG_LENGTH, f);
  fgets (buffer, LONG_LENGTH, f);

  for (i = 1; i <= catalog->nstruct; i++)
  {
    fgets (buffer, LONG_LENGTH, f);
    sscanf (buffer, "%d  %d", &tmpid, &nprogs);

    catalog->strctProps[i].NumProg = nprogs;
    if (nprogs > 0)
    {
      catalog->strctProps[i].ProgIDs   = (int *)    malloc (nprogs * sizeof(int));
      catalog->strctProps[i].ProgMrrts = (double *) malloc (nprogs * sizeof(double));
      for (j = 0; j < nprogs; j++)
      {
        fgets  (buffer, LONG_LENGTH, f);
        sscanf (buffer, "%d  %lf", &(catalog->strctProps[i].ProgIDs[j]), &(catalog->strctProps[i].ProgMrrts[j]));
      }
    }
  }

  fclose (f);
}
*/
