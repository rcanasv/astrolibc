CC     = gcc
INC    = -I/opt/gcc/5.4.0/hdf5/1.8.18/include
LIB    = -L/opt/gcc/5.4.0/hdf5/1.8.18/lib

analyze_galaxy_catalog: analyze_galaxy_catalog.c
	gcc analyze_galaxy_catalog.c -lm -o bin/analyze_galaxy_catalog

ctlgMatch: ctlgmatch.c archive.c catalog.c stf.c halomaker.c simulation.c particle.c ramses.c structure.c
	gcc ctlgmatch.c archive.c catalog.c stf.c  halomaker.c simulation.c particle.c ramses.c structure.c -lm -o bin/ctlgMatch

convert: convert.c archive.c catalog.c stf.c halomaker.c
	gcc convert.c archive.c catalog.c stf.c  halomaker.c -lm -o bin/convert

test: test.c archive.c catalog.c stf.c simulation.c particle.c structure.c hdf5.c hdf5sim.c
	$(CC) $(INC) $(LIB) test.c archive.c catalog.c stf.c simulation.c particle.c structure.c hdf5.c hdf5sim.c -lm -o bin/test
