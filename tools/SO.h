#ifndef SO_H
#define SO_H

#include "../src/base.h"
#include "../src/typedef.h"
#include "../src/archive.h"
#include "../src/catalog.h"
#include "../src/simulation.h"
//#include <mpi.h>
#include <time.h>

void get_structure_SO (Catalog * ctlg, Simulation * sim, int * tasks);

#endif    /*  SO_H  */
