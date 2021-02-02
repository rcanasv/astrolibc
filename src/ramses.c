/*
 *  \file   ramses.c
 *  \brief  This file contains RAMSES input-output routines
 *
 *
 */

#include "base.h"
#include "ramses.h"


void ramses_hydro_read (Simulation * ramses, int filenum, Grid * grid)
{
  FILE  * f;
  int     ncpu;
  int     ndim;
  int     nvarh;
  int     nlevelmax;
  int     ngridmax;
  int     nboundary;
  int     i, j, k, l, n;
  int     dummy;
  int     dummyi;
  double  dummyd;
  double  gamma;
  int     twondim, twotondim;
  double  xc[8][3];
  int     ix, iy, iz;
  double  dx;
  int     tmplvl;
  int     tmpng;

  char    fname    [NAME_LENGTH];
  char    dummys   [NAME_LENGTH];
  char    buffer   [NAME_LENGTH];
  char    ordering [NAME_LENGTH];

  if (!(grid->alloc_ngrid && grid->alloc_level))
  {
    printf ("Can't read hydro without amr data\n");
    printf ("Load amr data first\n");
    printf ("Exiting\n");
    exit (0);
  }

  sprintf (fname, "%s/hydro_%s.out%05d", ramses->archive.path, ramses->archive.prefix, filenum+1);
  if ((f = fopen (fname, "r")) == NULL)
  {
    printf ("Couldn't open file %s\n", fname);
    exit (0);
  }

  RMSSSKIP  fread(&ncpu,  sizeof(int), 1, f);   RMSSSKIP   // ncpu
  RMSSSKIP  fread(&nvarh, sizeof(int), 1, f);   RMSSSKIP   // nvarh
  RMSSSKIP  fread(&ndim,  sizeof(int), 1, f);   RMSSSKIP   // ndim

  twondim   = 2 * ndim;
  twotondim = (int) pow (2, ndim);

  RMSSSKIP  fread(&nlevelmax, sizeof(int),    1, f);   RMSSSKIP   // nlevelmax
  RMSSSKIP  fread(&nboundary, sizeof(int),    1, f);   RMSSSKIP   // nboundary
  RMSSSKIP  fread(&gamma,     sizeof(double), 1, f);   RMSSSKIP   // gamma

  for (i = 0; i < nlevelmax; i++)
  {
    // Geometry
    dx = pow(0.5,i+1);
    for (j = 0; j < twotondim; j++)
    {
      iz =  j / 4;
      iy = (j -4*iz) / 2;
      ix = (j -4*iz -2*iy);
      xc[j][0] = ((double)ix - 0.5)*dx;
      xc[j][1] = ((double)iy - 0.5)*dx;
      xc[j][2] = ((double)iz - 0.5)*dx;
    }

    // Loop over cpus
    for (j = 0; j < ncpu; j++)
    {
      // Level
      RMSSSKIP
      fread (&tmplvl, sizeof(int), 1, f);
      RMSSSKIP

      // Ngrids
      RMSSSKIP
      fread (&tmpng, sizeof(int), 1, f);
      RMSSSKIP

      if (j != filenum)
      {
        if (tmpng > 0)
        {
          for (k = 0; k < twotondim; k++)
          {
            for (l = 0; l < nvarh; l++)
            {
              RMSSSKIP
              fseek (f, dummy, SEEK_CUR);
              RMSSSKIP
            }
          }
        }
      }
      else
      {
        if (tmpng > 0)
        {
          if (tmpng != grid->ngrid[i][filenum])
          {
            printf ("ngrid doesn't match, aborting\n");
            exit (0);
          }

          for (n = 0; n < tmpng; n++)
            grid->level[i].cell[n].dx = dx;

          for (k = 0; k < twotondim; k++)
          {
            // Compute oct center
            for (n = 0; n < tmpng; n++)
            {
              grid->level[i].cell[n].octPos[k][0] = grid->level[i].cell[n].Pos[0] + xc[k][0];
              grid->level[i].cell[n].octPos[k][1] = grid->level[i].cell[n].Pos[1] + xc[k][1];
              grid->level[i].cell[n].octPos[k][2] = grid->level[i].cell[n].Pos[2] + xc[k][2];
            }

            // Check if cell is refined
            for (n = 0; n < tmpng; n++)
              grid->level[i].cell[n].okOct[k] = ((grid->level[i].cell[n].sonIndex[k] == 0) || (tmplvl == nlevelmax));

            for (l = 0; l < nvarh; l++)
            {
              //
              // Only read density for the moment
              // Skipping the rest of the properties
              //
              RMSSSKIP
              if (l == 0)
                for (n = 0; n < tmpng; n++)
                  fread (&grid->level[i].cell[n].octRho[k], sizeof(double), 1, f);
              else
                fseek (f, dummy, SEEK_CUR);
              RMSSSKIP
            }
          } // k loop
        } // if tmpng > 0
      } // else - Actually reading stuff
    } // cpus loop
  } // nlevelmax loop
}




void ramses_amr_init (Grid * grid)
{
  grid->alloc_ngrid = 0;
  grid->alloc_level = 0;
  grid->level = NULL;
  grid->ngrid = NULL;
}


void ramses_amr_free (Grid * grid)
{
  int k;

  if (grid->alloc_level)
  {
    for (k = 0; k < grid->nlevelmax; k++)
      if (grid->level[k].num)
        free (grid->level[k].cell);
    free (grid->level);
  }

  if (grid->alloc_ngrid)
  {
    for (k = 0; k < grid->nlevelmax; k++)
      free (grid->ngrid[k]);
    free (grid->ngrid);
  }
}


void ramses_amr_load (Simulation * ramses, int filenum, Grid * grid)
{
  if (grid->alloc_level || grid->alloc_ngrid)
    return;

  FILE  * f;
  int     mycpu;
  int     bsize;
  int     ncpu;
  int     ndim;
  int     nx [3];
  int     nlevelmax;
  int     ngridmax;
  int     nboundary;
  int     ngrid_current;
  double  boxlen;
  int     i, n, j, k;
  int     dummy;
  int     dummyi;
  double  dummyd;
  char    name[100];
  int     twondim, twotondim;

  char    fname    [NAME_LENGTH];
  char    dummys   [NAME_LENGTH];
  char    buffer   [NAME_LENGTH];
  char    ordering [NAME_LENGTH];

  sprintf (fname, "%s/amr_%s.out%05d", ramses->archive.path, ramses->archive.prefix, filenum+1);
  if ((f = fopen (fname, "r")) == NULL)
  {
    printf ("Couldn't open file %s\n", fname);
    exit (0);
  }

  RMSSSKIP   fread (&ncpu,  sizeof(int), 1, f);   RMSSSKIP  // ncpu
  RMSSSKIP   fread (&ndim,  sizeof(int), 1, f);   RMSSSKIP  // ndim

  twondim   = 2 * ndim;
  twotondim = (int) pow (2, ndim);

  RMSSSKIP   fread (&nx,            sizeof(nx),     1, f);   RMSSSKIP  // nx, ny, nz
  RMSSSKIP   fread (&nlevelmax,     sizeof(int),    1, f);   RMSSSKIP  // nlevelmax
  RMSSSKIP   fread (&ngridmax,      sizeof(int),    1, f);   RMSSSKIP  // ngridmax
  RMSSSKIP   fread (&nboundary,     sizeof(int),    1, f);   RMSSSKIP  // nboundary
  RMSSSKIP   fread (&ngrid_current, sizeof(int),    1, f);   RMSSSKIP  // ngrid_current
  RMSSSKIP   fread (&boxlen,        sizeof(double), 1, f);   RMSSSKIP  // boxlen
  RMSSSKIP   fread (&nx,            sizeof(nx),     1, f);   RMSSSKIP  // noutput, iout, ifout
  RMSSSKIP   fseek (f, dummy, SEEK_CUR);                     RMSSSKIP  // tout
  RMSSSKIP   fseek (f, dummy, SEEK_CUR);                     RMSSSKIP  // aout
  RMSSSKIP   fseek (f, dummy, SEEK_CUR);                     RMSSSKIP  // t
  RMSSSKIP   fseek (f, dummy, SEEK_CUR);                     RMSSSKIP  // dtold
  RMSSSKIP   fseek (f, dummy, SEEK_CUR);                     RMSSSKIP  // dtnew
  RMSSSKIP   fseek (f, dummy, SEEK_CUR);                     RMSSSKIP  // nstep, nstep_coarse
  RMSSSKIP   fseek (f, dummy, SEEK_CUR);                     RMSSSKIP  // einit, mass_tot_0, rho_tot
  RMSSSKIP   fseek (f, dummy, SEEK_CUR);                     RMSSSKIP  // omega_m, omega_l, omega_k, omega_b, h0, aexp_ini, boxlen_ini
  RMSSSKIP   fseek (f, dummy, SEEK_CUR);                     RMSSSKIP  // aexp, hexp, aexp_old, epot_tot_init, epot_tot_old
  RMSSSKIP   fseek (f, dummy, SEEK_CUR);                     RMSSSKIP  // mass_sph
  RMSSSKIP   fseek (f, dummy, SEEK_CUR);                     RMSSSKIP  // headl Nproc x Nlevel  array
  RMSSSKIP   fseek (f, dummy, SEEK_CUR);                     RMSSSKIP  // taill Nproc x Nlevel  array

  // Allocate memory
  grid->ncpu      = ncpu;
  grid->nlevelmax = nlevelmax;
  grid->ngrid     = (int **) malloc (nlevelmax * (sizeof(int *)));
  for (k = 0; k < nlevelmax; k++)
    grid->ngrid[k] = (int *) malloc ((ncpu+nboundary) * sizeof(int));
  grid->alloc_ngrid = 1;

  // Numbl
  RMSSSKIP
  for (k = 0; k < nlevelmax; k++)
    for (j = 0; j < ncpu; j++)
      fread (&grid->ngrid[k][j], sizeof(int), 1, f);
  RMSSSKIP

  RMSSSKIP   fseek (f, dummy, SEEK_CUR);   RMSSSKIP  // numbtot

  if (nboundary > 0)
  {
    RMSSSKIP   fseek (f, dummy, SEEK_CUR);   RMSSSKIP  // headb
    RMSSSKIP   fseek (f, dummy, SEEK_CUR);   RMSSSKIP  // tailb
    // numb
    RMSSSKIP
    for (k = 0; k < nlevelmax; k++)
      for (j = ncpu; j < nboundary; j++)
        fread (&grid->ngrid[k][j], sizeof(int), 1, f);
    RMSSSKIP
  }

  // Free memory
  RMSSSKIP   fseek (f, dummy, SEEK_CUR);   RMSSSKIP    // headf, tailf, numb, used_mem, used_mem_tot

  // Cpu boundaries
  //
  // Ordering
  //
  RMSSSKIP   fread (&ordering, dummy, 1, f);   RMSSSKIP
  if (!strncmp("bisection", ordering, 9))
  {
    RMSSSKIP   fseek (f, dummy, SEEK_CUR);   RMSSSKIP
    RMSSSKIP   fseek (f, dummy, SEEK_CUR);   RMSSSKIP
    RMSSSKIP   fseek (f, dummy, SEEK_CUR);   RMSSSKIP
    RMSSSKIP   fseek (f, dummy, SEEK_CUR);   RMSSSKIP
    RMSSSKIP   fseek (f, dummy, SEEK_CUR);   RMSSSKIP
  }
  else
  {
    RMSSSKIP   fseek (f, dummy, SEEK_CUR);   RMSSSKIP  // Indices as info.txt
  }

  RMSSSKIP   fseek (f, dummy, SEEK_CUR);   RMSSSKIP   // Coarse level son
  RMSSSKIP   fseek (f, dummy, SEEK_CUR);   RMSSSKIP   // Coarse level flag1
  RMSSSKIP   fseek (f, dummy, SEEK_CUR);   RMSSSKIP   // Coarse level cpu_map

  // Skip for cells that this cpu is not reading
  int toskip = 3 + ndim + 1 + twondim + 3*twotondim;

  // Allocate memory
  grid->level = (Level *) malloc (grid->nlevelmax * sizeof(Level));
  for (k = 0; k < grid->nlevelmax; k++)
  {
    grid->level[k].num = grid->ngrid[k][filenum];
    if (grid->level[k].num)
      grid->level[k].cell = (Cell *) malloc (grid->level[k].num * sizeof(Cell));
  }
  grid->alloc_level = 1;

  // Load data
  for (k = 0; k < grid->nlevelmax; k++)
  {
    for (j = 0; j < grid->ncpu; j++)
    {
      if (grid->ngrid[k][j])
      {
        // If not cpu skip
        if (j != filenum)
        {
          for (n = 0; n < toskip; n++)
          {
            RMSSSKIP
            fseek (f, dummy, SEEK_CUR);
            RMSSSKIP
          }
        }
        else
        {
          // Grid Index
          RMSSSKIP
          for (n = 0; n < grid->level[k].num; n++)
            fread (&grid->level[k].cell[n].myIndex, sizeof(int), 1, f);
          RMSSSKIP

          // Next Index
          RMSSSKIP
          for (n = 0; n < grid->level[k].num; n++)
            fread (&grid->level[k].cell[n].nextIndex, sizeof(int), 1, f);
          RMSSSKIP

          // Prev Index
          RMSSSKIP
          for (n = 0; n < grid->level[k].num; n++)
            fread (&grid->level[k].cell[n].prevIndex, sizeof(int), 1, f);
          RMSSSKIP

          // Grid Position
          for (i = 0; i < ndim; i++)
          {
            RMSSSKIP
            for (n = 0; n < grid->level[k].num; n++)
              fread (&grid->level[k].cell[n].Pos[i], sizeof(double), 1, f);
            RMSSSKIP
          }

          // Father Index
          RMSSSKIP
          for (n = 0; n < grid->level[k].num; n++)
            fread (&grid->level[k].cell[n].fatherIndex, sizeof(int), 1, f);
          RMSSSKIP

          // Neighbour
          for (i = 0; i < twondim; i++)
          {
            RMSSSKIP
            for (n = 0; n < grid->level[k].num; n++)
              fread (&grid->level[k].cell[n].nborIndex[i], sizeof(int), 1, f);
            RMSSSKIP
          }

          // Son
          for (i = 0; i < twotondim; i++)
          {
            RMSSSKIP
            for (n = 0; n < grid->level[k].num; n++)
              fread (&grid->level[k].cell[n].sonIndex[i], sizeof(int), 1, f);
            RMSSSKIP
          }

          // Cpu map
          for (i = 0; i < twotondim; i++)
          {
            RMSSSKIP
            for (n = 0; n < grid->level[k].num; n++)
              fread (&grid->level[k].cell[n].cpuMap[i], sizeof(int), 1, f);
            RMSSSKIP
          }

          // Ref Map
          for (i = 0; i < twotondim; i++)
          {
            RMSSSKIP
            for (n = 0; n < grid->level[k].num; n++)
              fread (&grid->level[k].cell[n].refMap[i], sizeof(int), 1, f);
            RMSSSKIP
          }
        }
      }
    } // ncpu
  } // nlevelmax

  fclose (f);
}



void ramses_init (Simulation * ramses)
{
  int     i;
  FILE  * f;
  char    fname  [NAME_LENGTH];
  char    dummys [NAME_LENGTH];
  char    buffer [NAME_LENGTH];


  sprintf (fname, "%s/info_%s.txt", ramses->archive.path, ramses->archive.prefix);
  if ((f = fopen (fname, "r")) == NULL)
  {
    printf ("Couldn't open file %s\n", fname);
    exit (0);
  }
  for (i = 0; i < 7; i++)
    fgets (buffer, 100, f);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->Lbox);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->Time);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->a);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->cosmology.HubbleParam);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->cosmology.OmegaM);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->cosmology.OmegaL);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %s ", dummys, dummys, dummys);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->cosmology.OmegaB);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->unit_l);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->unit_d);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->unit_t);
  fclose (f);

  ramses->unit_m = ramses->unit_d * ramses->unit_l * ramses->unit_l * ramses->unit_l;  // in grams
  ramses->unit_m = ramses->unit_m / 1.989e+33;                                         // in solar masses

  ramses->unit_v = ramses->unit_l / ramses->unit_t;                                    // in cm / s
  ramses->unit_v = ramses->unit_v / 100000.0;                                          // in km / s

  ramses->unit_l = ramses->unit_l / 3.08e+21;                                          // in kpc
  ramses->h = ramses->cosmology.HubbleParam / 100.0;

  // Box is now in kpc
  ramses->Lbox *= ramses->unit_l;

  ramses->npartinfile = NULL;
}


//
//  Read RAMSES Particle FILE
//
void ramses_load_particles (Simulation * ramses, int filenum, Particle ** part)
{

  int     i, j;
  char    fname  [NAME_LENGTH];
  char    buffer [NAME_LENGTH];

  int     dummy;
  int     dummyi;
  float   dummyf;
  double  dummyd;
  char    dummys [NAME_LENGTH];

  FILE * f;

  Particle * P;

  //
  //  Read Info file to get Simulation info
  //
  sprintf (fname, "%s/info_%s.txt", ramses->archive.path, ramses->archive.prefix);
  if ((f = fopen (fname, "r")) == NULL)
  {
    printf ("Couldn't open file %s\n", fname);
    exit (0);
  }
  for (i = 0; i < 7; i++)
    fgets (buffer, 100, f);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->Lbox);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->Time);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->a);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->cosmology.HubbleParam);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->cosmology.OmegaM);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->cosmology.OmegaL);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %s ", dummys, dummys, dummys);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->cosmology.OmegaB);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->unit_l);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->unit_d);
  fgets(buffer, 100, f);   sscanf (buffer, "%s %s %lf", dummys, dummys, &ramses->unit_t);
  fclose (f);

  ramses->unit_m = ramses->unit_d * ramses->unit_l * ramses->unit_l * ramses->unit_l;  // in grams
  ramses->unit_m = ramses->unit_m / 1.989e+33;                                         // in solar masses

  ramses->unit_v = ramses->unit_l / ramses->unit_t;                                    // in cm / s
  ramses->unit_v = ramses->unit_v / 100000.0;                                          // in km / s

  ramses->unit_l = ramses->unit_l / 3.08e+21;                                          // in kpc

  // Box is now in kpc
  ramses->Lbox *= ramses->unit_l;

  // H0 -> h
  ramses->h = ramses->cosmology.HubbleParam / 100.0;

  //
  //  Read Particle file to get Simulation info
  //
  if (ramses->format == RAMSES)       sprintf (fname, "%s/part_%s.out%05d", ramses->archive.path, ramses->archive.prefix, filenum+1);
  if (ramses->format == RAMSES_STAR)  sprintf (fname, "%s/star_%s.out%05d", ramses->archive.name, ramses->archive.prefix, filenum+1);
  f = fopen(fname,"r");
  if ((f = fopen (fname, "r")) == NULL)
  {
    printf ("Couldn't open file %s\n", fname);
    exit (0);
  }

  //!--- Header
  RMSSSKIP  fread(&ramses->ncpu,     sizeof(int),    1, f);  RMSSSKIP
  RMSSSKIP  fread(&ramses->ndim,     sizeof(int),    1, f);  RMSSSKIP
  RMSSSKIP  fread(&ramses->npart,    sizeof(int),    1, f);  RMSSSKIP

  if (ramses->npartinfile == NULL)
    ramses->npartinfile = (int *) malloc (ramses->ncpu * sizeof(int));
  ramses->npartinfile[filenum] = ramses->npart;

  if (ramses->format == RAMSES)
  {
    RMSSSKIP  fread(&ramses->seed[0],  sizeof(int),    4, f);  RMSSSKIP
  }

  RMSSSKIP  fread(&ramses->nstarTot, sizeof(int),    1, f);  RMSSSKIP

  if (ramses->format == RAMSES)
  {
    RMSSSKIP  fread(&ramses->mstarTot, sizeof(double), 1, f);  RMSSSKIP
    RMSSSKIP  fread(&ramses->mstarLst, sizeof(double), 1, f);  RMSSSKIP
    RMSSSKIP  fread(&ramses->nsink,    sizeof(int),    1, f);  RMSSSKIP
  }

/*
   printf ("NumProcs        %d\n", ramses->ncpu);
   printf ("Num Dims        %d\n", ramses->ndim);
   printf ("Npart           %d\n", ramses->npart);
   for (i = 0; i < 4; i++)
     printf ("LocalSeed[%d]    %d\n", i, ramses->seed[i]);
   printf ("NstarTot        %d\n", ramses->nstarTot);
   printf ("Mstar_tot       %g\n", ramses->mstarTot);
   printf ("Mstar_lost      %g\n", ramses->mstarLst);
   printf ("NumSink         %d\n", ramses->nsink);
*/

  if ((*(part) = (Particle *) malloc (ramses->npart * sizeof(Particle))) == NULL)
  {
    printf ("Couldn't allocate memory for Particle array\n");
    exit(0);
  }

  P = *(part);

  //--- Pos
  for (i = 0; i < ramses->ndim; i++)
  {
    RMSSSKIP
    for (j = 0; j < ramses->npart; j++)
    {
      fread(&dummyd, sizeof(double), 1, f);
      P[j].Pos[i] = dummyd;
    }
    RMSSSKIP
  }

  //--- Vel
  for (i = 0; i < ramses->ndim; i++)
  {
    RMSSSKIP
    for (j = 0; j < ramses->npart; j++)
    {
      fread(&dummyd, sizeof(double), 1, f);
      P[j].Vel[i] = dummyd;
    }
    RMSSSKIP
  }

  //--- Mass
  RMSSSKIP
  for (j = 0; j < ramses->npart; j++)
  {
    fread(&dummyd, sizeof(double), 1, f);
    P[j].Mass = dummyd;
  }
  RMSSSKIP

  //--- Id
  RMSSSKIP
  for (j = 0; j < ramses->npart; j++)
  {
    fread(&dummyi, sizeof(int), 1, f);
    P[j].Id = dummyi;
  }
  RMSSSKIP

  //--- Level
  RMSSSKIP
  for (j = 0; j < ramses->npart; j++)
  {
    fread(&dummyi, sizeof(int), 1, f);
    P[j].Level = dummyi;
  }
  RMSSSKIP

  //--- Birth Epoch
  RMSSSKIP
  for (j = 0; j < ramses->npart; j++)
  {
    fread(&dummyd, sizeof(double), 1, f);
    P[j].Age = dummyd;
  }
  RMSSSKIP

  //--- Metallicity if ((STAR || SINK) && (METAL))
  //
  //  Skip for the moment.
  //
  //for (i = 0; i < 11; i++)
  //{
  //    SKIP  fread(&ramses_met[i][0], sizeof(double), npart, f);  SKIP
  //}
  fclose (f);

  //
  // Convert to human readable units
  //
  for (i = 0; i < ramses->npart; i++)
  {
    P[i].Pos[0] *= ramses->unit_l;
    P[i].Pos[1] *= ramses->unit_l;
    P[i].Pos[2] *= ramses->unit_l;

    P[i].Vel[0] *= ramses->unit_v;
    P[i].Vel[1] *= ramses->unit_v;
    P[i].Vel[2] *= ramses->unit_v;

    P[i].Mass   *= ramses->unit_m;

    P[i].Type    = 1;
  }
}


double friedman(double Omega0, double OmegaL, double OmegaK, double alpha, double axp_min, double ** axp_out, double ** hexp_out, double ** tau_out, double ** t_out, int ntable)
{
  /*
  ! ######################################################!
  ! This subroutine assumes that axp = 1 at z = 0 (today) !
  ! and that t and tau = 0 at z = 0 (today).              !
  ! axp is the expansion factor, hexp the Hubble constant !
  ! defined as hexp=1/axp*daxp/dtau, tau the conformal    !
  ! time, and t the look-back time, both in unit of 1/H0. !
  ! alpha is the required accuracy and axp_min is the     !
  ! starting expansion factor of the look-up table.       !
  ! ntable is the required size of the look-up table.     !
  ! ######################################################!
  */
  double axp_tau = 1.0;
  double axp_t   = 1.0;
  double tau     = 0.0;
  double t       = 0.0;
  double age_tot;

  double dtau;
  double dt;
  double axp_tau_pre;
  double axp_t_pre;

  int nstep = 0;
  int nskip;
  int nout;

  *axp_out  = (double *) malloc (ntable * sizeof(double));
  *hexp_out = (double *) malloc (ntable * sizeof(double));
  *tau_out  = (double *) malloc (ntable * sizeof(double));
  *t_out    = (double *) malloc (ntable * sizeof(double));

  while ((axp_tau >= axp_min) || (axp_t >= axp_min))
  {
    nstep++;
    dtau        = alpha * axp_tau / dadtau(axp_tau, Omega0, OmegaL, OmegaK);
    axp_tau_pre = axp_tau - dadtau(axp_tau, Omega0, OmegaL, OmegaK) * dtau / 2.0;
    axp_tau     = axp_tau - dadtau(axp_tau_pre, Omega0, OmegaL, OmegaK) * dtau;
    tau         = tau - dtau;

    dt          = alpha * axp_t / dadt(axp_t, Omega0, OmegaL, OmegaK);
    axp_t_pre   = axp_t - dadt(axp_t, Omega0, OmegaL, OmegaK) * dt / 2.0;
    axp_t       = axp_t - dadt(axp_t_pre, Omega0, OmegaL, OmegaK) * dt;
    t           = t - dt;
  }

  age_tot =-t;
//   printf ("Age of the universe (in unit of 1/H0)=%e \n", -t);

  nskip = nstep / ntable;

  axp_t   = 1.0;
  axp_tau = 1.0;
  tau     = 0.0;
  t       = 0.0;

  nstep = 0;
  nout  = 0;

  t_out   [0][nout] = t;
  tau_out [0][nout] = tau;
  axp_out [0][nout] = axp_tau;
  hexp_out[0][nout] = dadtau (axp_tau, Omega0, OmegaL, OmegaK) / axp_tau;


  while ((axp_tau >= axp_min) || (axp_t >= axp_min))
  {
    nstep++;
    dtau        = alpha * axp_tau / dadtau (axp_tau, Omega0, OmegaL, OmegaK);
    axp_tau_pre = axp_tau - dadtau(axp_tau, Omega0, OmegaL, OmegaK) * dtau/2.0;
    axp_tau     = axp_tau - dadtau(axp_tau_pre, Omega0, OmegaL, OmegaK) * dtau;
    tau         = tau - dtau;

    dt          = alpha * axp_t / dadt(axp_t, Omega0, OmegaL, OmegaK);
    axp_t_pre   = axp_t - dadt(axp_t, Omega0, OmegaL, OmegaK) * dt / 2.0;
    axp_t       = axp_t - dadt(axp_t_pre, Omega0, OmegaL, OmegaK) * dt;
    t           = t -dt;

    if ((nstep%nskip) == 0)
    {
      nout = nout + 1;
      t_out   [0][nout] = t;
      tau_out [0][nout] = tau;
      axp_out [0][nout] = axp_tau;
      hexp_out[0][nout] = dadtau(axp_tau, Omega0, OmegaL, OmegaK) / axp_tau;
    }
  }

  t_out   [0][ntable-1] = t;
  tau_out [0][ntable-1] = tau;
  axp_out [0][ntable-1] = axp_tau;
  hexp_out[0][ntable-1] = dadtau(axp_tau, Omega0, OmegaL, OmegaK) / axp_tau;

  return age_tot;
}



double dadtau (double axp_tau, double Omega0, double OmegaL, double OmegaK)
{
  return sqrt(axp_tau*axp_tau*axp_tau * (Omega0 + OmegaL*axp_tau*axp_tau*axp_tau + OmegaK*axp_tau));
}



double dadt (double axp_t, double Omega0, double OmegaL, double OmegaK)
{
  return sqrt((1.0/axp_t) * (Omega0 + OmegaL*axp_t*axp_t*axp_t + OmegaK*axp_t));
}



void  ramses_structure_calculate_star_age (Simulation * ramses, Catalog * ctlg, int * strct_to_get)
{
  int      i, j, k;

  double * axp_frw;
  double * hexp_frw;
  double * tau_frw;
  double * t_frw;
  int      n_frw = 1000;
  double   t;

  double   time_tot;
  double   time_simu;

  FILE   * f;
  char     fname  [NAME_LENGTH];
  char     dummys [NAME_LENGTH];
  char     buffer [NAME_LENGTH];
  double   dummyd;

  Structure * strct;

  //
  //  Fill frw arrays
  //
  time_tot = friedman(ramses->cosmology.OmegaM, ramses->cosmology.OmegaL, 0.0, 1e-6, 1e-3, &axp_frw, &hexp_frw, &tau_frw, &t_frw, n_frw);

  // Find neighbouring conformal time
  i = 1;
  while (tau_frw[i] > ramses->Time  && i < n_frw)
    i = i+1;

  // Interpolate time
  time_simu = t_frw[i]   * (ramses->Time - tau_frw[i-1]) / (tau_frw[i]   - tau_frw[i-1]) + \
              t_frw[i-1] * (ramses->Time - tau_frw[i])   / (tau_frw[i-1] - tau_frw[i]);


  ramses->LookBackTime = (time_tot + time_simu) / (ramses->cosmology.HubbleParam*1e5/3.08e24) / (365*24*3600*1e9);
  printf ("Time simu      %lf\n", (time_tot + time_simu) / (ramses->cosmology.HubbleParam*1e5/3.08e24) / (365*24*3600*1e9));
  printf ("Hubble time    %lf\n", time_tot / (ramses->cosmology.HubbleParam*1e5/3.08e24) / (365*24*3600*1e9));
  printf ("i               %d\n", i);
  printf ("time            %e\n", ramses->Time);
  printf ("time_tot        %e\n", time_tot);
  printf ("time_simu       %e\n", time_simu);
  printf ("t_tot + t_simu  %e\n", time_tot + time_simu);


  for (i = 1; i <= ctlg->nstruct; i++)
  {
    if (strct_to_get[i])
    {
      if (ctlg->strctProps[i].iPart)
      {
        strct = &ctlg->strctProps[i];
        for (j = 0; j < strct->NumPart; j++)
        {
          if (strct->Part[j].Age != 0)
          {
            k = 1;
            while (tau_frw[k] > strct->Part[j].Age  && k < n_frw)
              k++;

            t = t_frw[k]   * (strct->Part[j].Age - tau_frw[k-1]) / (tau_frw[k]   - tau_frw[k-1]) + \
                t_frw[k-1] * (strct->Part[j].Age - tau_frw[k])   / (tau_frw[k-1] - tau_frw[k]);

           // Age in years
            strct->Part[j].Age = (time_simu - t) / (ramses->cosmology.HubbleParam*1e5/3.08e24) / (365*24*3600.0);
          }
        }
      }
      else
      {
        printf ("Somethin weird happened. Particle array for Structure has not been initialized\n");
        exit (0);
      }
    }
  }

  free (axp_frw);
  free (hexp_frw);
  free (tau_frw);
  free (t_frw);
}



void  ramses_catalog_calculate_star_age (Simulation * ramses, Catalog * ctlg)
{
  int      i, j, k;

  double * axp_frw;
  double * hexp_frw;
  double * tau_frw;
  double * t_frw;
  int      n_frw = 1000;
  double   t;

  double   time_tot;
  double   time_simu;

  FILE   * f;
  char     fname  [NAME_LENGTH];
  char     dummys [NAME_LENGTH];
  char     buffer [NAME_LENGTH];
  double   dummyd;

  Structure * strct;

  //
  //  Fill frw arrays
  //
  time_tot = friedman(ramses->cosmology.OmegaM, ramses->cosmology.OmegaL, 0.0, 1e-6, 1e-3, &axp_frw, &hexp_frw, &tau_frw, &t_frw, n_frw);

  // Find neighbouring conformal time
  i = 1;
  while (tau_frw[i] > ramses->Time  && i < n_frw)
    i = i+1;

  // Interpolate time
  time_simu = t_frw[i]   * (ramses->Time - tau_frw[i-1]) / (tau_frw[i]   - tau_frw[i-1]) + \
              t_frw[i-1] * (ramses->Time - tau_frw[i])   / (tau_frw[i-1] - tau_frw[i]);


  ramses->LookBackTime = (time_tot + time_simu) / (ramses->cosmology.HubbleParam*1e5/3.08e24) / (365*24*3600*1e9);
  printf ("Time simu    %lf\n", (time_tot + time_simu) / (ramses->cosmology.HubbleParam*1e5/3.08e24) / (365*24*3600*1e9));
  printf ("Hubble time  %lf\n", time_tot / (ramses->cosmology.HubbleParam*1e5/3.08e24) / (365*24*3600*1e9));
  printf ("i               %d\n", i);
  printf ("time            %e\n", ramses->Time);
  printf ("time_tot        %e\n", time_tot);
  printf ("time_simu       %e\n", time_simu);
  printf ("t_tot + t_simu  %e\n", time_tot + time_simu);


  for (i = 1; i <= ctlg->nstruct; i++)
  {
    if (ctlg->strctProps[i].iPart)
    {
      strct = &ctlg->strctProps[i];
      for (j = 0; j < strct->NumPart; j++)
      {
        if (strct->Part[j].Age != 0)
        {
          k = 1;
          while (tau_frw[k] > strct->Part[j].Age  && k < n_frw)
            k++;

          t = t_frw[k]   * (strct->Part[j].Age - tau_frw[k-1]) / (tau_frw[k]   - tau_frw[k-1]) + \
              t_frw[k-1] * (strct->Part[j].Age - tau_frw[k])   / (tau_frw[k-1] - tau_frw[k]);

          // Age in years
          strct->Part[j].Age = (time_simu - t) / (ramses->cosmology.HubbleParam*1e5/3.08e24) / (365*24*3600.0);
        }
      }
    }
    else
    {
      printf ("Somethin weird happened. Particle array for Structure has not been initialized\n");
      exit (0);
    }
  }


  free (axp_frw);
  free (hexp_frw);
  free (tau_frw);
  free (t_frw);
}
