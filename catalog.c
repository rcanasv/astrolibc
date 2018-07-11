/*
 *  \file   catalog.c
 *  \brief  This file contains functions for Catalog struct type.
 *
 *
 *
 */

#include "catalog.h"


void Catalog_init (Catalog * ctlg)
{
  ctlg->nprocs     = 0;
  ctlg->nstruct    = 0;
  ctlg->iprops     = 0;
  ctlg->iparts     = 0;
  ctlg->strctProps = NULL;
  ctlg->strctParts = NULL;

  if ((!strcmp(ctlg->archive.format, "stf"))          ||  \
      (!strcmp(ctlg->archive.format, "STF"))          ||  \
      (!strcmp(ctlg->archive.format, "VELOCIraptor")) ||  \
      (!strcmp(ctlg->archive.format, "velociraptor")))
    ctlg->format = STF;
  else
    if ((!strcmp(ctlg->archive.format, "hmkr"))            ||  \
        (!strcmp(ctlg->archive.format, "hmkr_treebricks")) ||  \
        (!strcmp(ctlg->archive.format, "HaloMaker"))       ||  \
        (!strcmp(ctlg->archive.format, "halomaker")))
      ctlg->format = HALOMAKER;
    else
    {
      printf ("Format %s not supported\n", ctlg->archive.format);
      printf ("Exiting...\n");
      exit (0);
    }
}



void Catalog_load (Catalog * ctlg)
{
  Catalog_load_properties (ctlg);
  Catalog_load_particles  (ctlg);
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
  if (ctlg->format == HALOMAKER)   halomaker_read_particles (ctlg);
}



void Catalog_get_particle_properties (Catalog * ctlg, Simulation * sim)
{
  if (ctlg->format == STF)        stf_catalog_get_particle_properties       (ctlg, sim);
  if (ctlg->format == HALOMAKER)  halomaker_catalog_get_particle_properties (ctlg, sim);
}



void Catalog_free (Catalog * ctlg)
{
  int i;

  if (ctlg->iprops)
  {
    for (i = 1; i <= ctlg->nstruct; i++)
    {
      if (ctlg->strctProps[i].SubIDs != NULL)
        free (ctlg->strctProps[i].SubIDs);

      if (ctlg->strctProps[i].MatchIDs != NULL)
        free (ctlg->strctProps[i].MatchIDs);

      if (ctlg->strctProps[i].MatchMrrts != NULL)
        free (ctlg->strctProps[i].MatchMrrts);
    }
    free (ctlg->strctProps);
  }


  if (ctlg->iparts)
  {
    for (i = 1; i <= ctlg->nstruct; i++)
      free (ctlg->strctParts[i]);
    free (ctlg->strctParts);
  }
}



void Catalog_fill_SubIDS (Catalog * ctlg)
{
  int i, j;
  int bob;
  int tmp;

  for (i = 1; i <= ctlg->nstruct; i++)
    if (ctlg->strctProps[i].HostID == -1)
    {
      bob = ctlg->strctProps[i].NumSubs;
      ctlg->strctProps[i].SubIDs = (int *) malloc (bob * sizeof(int));
      ctlg->strctProps[i].dummy = 0;
      for (j = 0; j < bob; j++)
        ctlg->strctProps[i].SubIDs[j] = 0;
    }

  for (i = 1; i <= ctlg->nstruct; i++)
  {
    if (ctlg->strctProps[i].HostID != -1 && \
        ctlg->strctProps[i].HostID != ctlg->strctProps[i].ID)
    {
      bob = i;
      while (ctlg->strctProps[bob].HostID != -1)
      {
        bob = ctlg->strctProps[bob].DirectHostID;
      }
      tmp = ctlg->strctProps[bob].dummy++;
      ctlg->strctProps[bob].SubIDs[tmp] = i;
    }
  }
}



void Catalog_fill_isolated (Catalog * ctlg)
{
  if (ctlg->format == STF)  stf_catalog_fill_isolated (ctlg);
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
