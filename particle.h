/*
 *  \file particle.h
 *  \brief  header file for particle properties structure
 *
 */


typedef struct pdata_d
{
  double   Pos[3];
  double   Vel[3];
  double   Mass;
  long     Id;
  double   Age;
  double   Metal;
  double   Chem;
  int      Type;
} pdata_d;


typedef struct pdata_s
{
  float   Pos[3];
  float   Vel[3];
  float   Mass;
  int     Id;
  float   Age;
  float   Metal;
  float   Chem;
  int     Type;
} pdata_s;
