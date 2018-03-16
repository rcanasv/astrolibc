CC          = gcc

HDF5_INCL   = -I/opt/gcc/5.4.0/zlib/1.2.11/include  -I/opt/gcc/5.4.0/hdf5/1.8.18/include
HDF5_LIB    = -L/opt/gcc/5.4.0/hdf5/1.8.18/lib  /opt/gcc/5.4.0/hdf5/1.8.18/lib/libhdf5_hl.a /opt/gcc/5.4.0/hdf5/1.8.18/lib/libhdf5.a -L/opt/gcc/5.4.0/zlib/1.2.11/lib
HDF5_FLAGS  = -lhdf5 -lhdf5_hl -lz -ldl -Wl,-rpath -Wl,/opt/gcc/5.4.0/hdf5/1.8.18/lib

INC         = $(HDF5_INCL)
LIB         = $(HDF5_LIB)
FLAGS       = -lm $(HDF5_FLAGS)

analyze_galaxy_catalog: analyze_galaxy_catalog.c
	gcc analyze_galaxy_catalog.c -lm -o bin/analyze_galaxy_catalog

ctlgMatch: ctlgmatch.c archive.c catalog.c stf.c halomaker.c simulation.c particle.c ramses.c structure.c
	gcc ctlgmatch.c archive.c catalog.c stf.c  halomaker.c simulation.c particle.c ramses.c structure.c -lm -o bin/ctlgMatch

convert: convert.c archive.c catalog.c stf.c halomaker.c
	gcc convert.c archive.c catalog.c stf.c  halomaker.c -lm -o bin/convert

test: test.c archive.c catalog.c stf.c simulation.c halomaker.c particle.c ramses.c structure.c hdf5routines.c hdf5sim.c
	$(CC) $(INC) $(LIB) test.c archive.c catalog.c stf.c halomaker.c ramses.c simulation.c particle.c structure.c hdf5routines.c hdf5sim.c -o bin/test $(FLAGS)
