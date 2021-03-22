ifeq ($(MPI),yes)
CC          = mpicc
else
CC          = gcc
MPI_INCL    =
MPI_LIB     =
endif

ifeq ($(MACHINE),horizon)
HDF5_INCL   = -I//softs/hdf5/1.8.18-gcc4/include
HDF5_LIB    = -L/softs/hdf5/1.8.18-gcc4/lib  /softs/hdf5/1.8.18-gcc4/lib/libhdf5_hl.a /softs/hdf5/1.8.18-gcc4/lib/libhdf5.a
HDF5_FLAGS  = -lhdf5 -lhdf5_hl -lz -ldl -Wl,-rpath -Wl,/softs/hdf5/1.8.18-gcc4/lib

GSL_INCL    = -I/softs/gsl/2.3/include
GSL_LIB     = -L/softs/gsl/2.3/lib
GSL_FLAGS   = -lgsl -lgslcblas

ifeq ($(MPI),yes)
MPI_INCL    = -I//softs/openmpi/1.6.3-ifort-12.1-torque-CentOS6/include
MPI_LIB     = -L/softs/openmpi/1.6.3-ifort-12.1-torque-CentOS6/lib
endif
endif

ifeq ($(MACHINE),raijin)
HDF5_INCL   = -I/home/571/rac571/opt/gcc-5.2.0/hdf5-1.8.18/include
HDF5_LIB    = -L/home/571/rac571/opt/gcc-5.2.0/hdf5-1.8.18/lib  /home/571/rac571/opt/gcc-5.2.0/hdf5-1.8.18/lib/libhdf5_hl.a /home/571/rac571/opt/gcc-5.2.0/hdf5-1.8.18/lib/libhdf5.a
HDF5_FLAGS  = -lhdf5 -lhdf5_hl -lz -ldl -Wl,-rpath -Wl,//home/571/rac571/opt/gcc-5.2.0/hdf5-1.8.18/lib

GSL_INCL    = -I//apps/gsl/2.3/include
GSL_LIB     = -L/apps/gsl/2.3/lib
GSL_FLAGS   = -lgsl -lgslcblas
endif

ifeq ($(MACHINE),icrar)
HDF5_INCL   = -I/opt/gcc/5.4.0/zlib/1.2.11/include  -I/opt/gcc/5.4.0/hdf5/1.8.18/include
HDF5_LIB    = -L/opt/gcc/5.4.0/hdf5/1.8.18/lib  /opt/gcc/5.4.0/hdf5/1.8.18/lib/libhdf5_hl.a /opt/gcc/5.4.0/hdf5/1.8.18/lib/libhdf5.a -L/opt/gcc/5.4.0/zlib/1.2.11/lib
HDF5_FLAGS  = -lhdf5 -lhdf5_hl -lz -ldl -Wl,-rpath -Wl,/opt/gcc/5.4.0/hdf5/1.8.18/lib

GSL_INCL    = -I/opt/gcc/5.4.0/gsl/2.2/include
GSL_LIB     = -L/opt/gcc/5.4.0/gsl/2.2/lib
GSL_FLAGS   = -lgsl -lgslcblas

endif

ifeq ($(MACHINE),pulsar)
HDF5_INCL   = -I/opt/gcc/5.4.0/zlib/1.2.11/include  -I/opt/gcc/5.4.0/hdf5/1.8.18/include
HDF5_LIB    = -L/opt/gcc/5.4.0/hdf5/1.8.18/lib  /opt/gcc/5.4.0/hdf5/1.8.18/lib/libhdf5_hl.a /opt/gcc/5.4.0/hdf5/1.8.18/lib/libhdf5.a -L/opt/gcc/5.4.0/zlib/1.2.11/lib
HDF5_FLAGS  = -lhdf5 -lhdf5_hl -lz -ldl -Wl,-rpath -Wl,/opt/gcc/5.4.0/hdf5/1.8.18/lib

GSL_INCL    = -I/opt/gsl/gsl-2.3_gcc-5.4.0/include
GSL_LIB     = -L/opt/gsl/gsl-2.3_gcc-5.4.0/lib
GSL_FLAGS   = -lgsl -lgslcblas

BZL_INCL    = -I/usr/include
BZL_LIB     = -L/usr/lib/x86_64-linux-gnu/
BZL_FLAGS   = -lbz2
endif

ifeq ($(MACHINE),hyades)
HDF5_INCL   = -I/home/rcanas/opt/gcc/6.3.0/hdf5/1.8.18/include
HDF5_LIB    = -L/home/rcanas/opt/gcc/6.3.0/hdf5/1.8.18/lib /home/rcanas/opt/gcc/6.3.0/hdf5/1.8.18/lib/libhdf5_hl.a /home/rcanas/opt/gcc/6.3.0/hdf5/1.8.18/lib/libhdf5.a
HDF5_FLAGS  = -lhdf5 -lhdf5_hl -lz -ldl -Wl,-rpath -Wl,/home/rcanas/opt/gcc/6.3.0/hdf5/1.8.18/lib
endif

ifeq ($(MACHINE),cosma)
HDF5_INCL   = -I/cosma/local/hdf5/intel_2018/1.8.20/include
HDF5_LIB    = -L/cosma/local/hdf5/intel_2018/1.8.20/lib /cosma/local/hdf5/intel_2018/1.8.20/lib/libhdf5_hl.a /cosma/local/hdf5/intel_2018/1.8.20/lib/libhdf5.a
HDF5_FLAGS  = -lhdf5 -lhdf5_hl -lz -ldl -Wl,-rpath -Wl,/cosma/local/hdf5/intel_2018/1.8.20/lib

GSL_INCL    = -I/cosma/local/gsl/2.4/include
GSL_LIB     = -L/cosma/local/gsl/2.4/lib
GSL_FLAGS   = -lgsl -lgslcblas
endif

ifeq ($(MACHINE),taurus)
HDF5_INCL   = -I/home2/rodrigo/opt/gcc/9.3.0/hdf5/1.8.18/include
HDF5_LIB    = -L/home2/rodrigo/opt/gcc/9.3.0/hdf5/1.8.18/lib /home2/rodrigo/opt/gcc/9.3.0/hdf5/1.8.18/lib/libhdf5_hl.a  /home2/rodrigo/opt/gcc/9.3.0/hdf5/1.8.18/lib/libhdf5.a
HDF5_FLAGS  = -lhdf5 -lhdf5_hl -lz -ldl

GSL_INCL    = -I/home2/rodrigo/opt/gcc/9.3.0/gsl/2.3/include
GSL_LIB     = -L/home2/rodrigo/opt/gcc/9.3.0/gsl/2.3/lib
GSL_FLAGS   = -lgsl -lgslcblas

BZL_INCL    = -I/usr/include
BZL_LIB     = -L/usr/lib/x86_64-linux-gnu/
BZL_FLAGS   = -lbz2
endif

INC         = $(HDF5_INCL) $(GSL_INCL) $(MPI_INCL) $(BZL_INCL)
LIB         = $(HDF5_LIB) $(GSL_LIB) $(MPI_LIB) $(BZL_LIB)
FLAGS       = -lm -O2 $(HDF5_FLAGS) $(GSL_FLAGS) $(BZL_FLAGS)


SRC   = $(wildcard src/*.c)
TOOLS = $(wildcard tools/*.c)

#
#  get_structure
#    Extract structure(s) from simulations as a Gadget-2 binay file.
#    Either to visualize only or to do postprocessing.
#
get_structure: apps/get_structure.c $(SRC)
		$(CC) $(INC) $(LIB) $^ -o bin/$@ $(FLAGS)


#
#  ihsc
#    Calculate IHSC mass fraction, and other relevant properties of 3DFOF
#    objects.
#
ihsc: apps/ihsc.c $(SRC) $(TOOLS)
		$(CC) $(INC) $(LIB) $^ -o bin/$@ $(FLAGS)

ihsc_AHF_SO: apps/ihsc_AHF_SO.c $(SRC) $(TOOLS)
		$(CC) $(INC) $(LIB) $^ -o bin/$@ $(FLAGS)

#
#		check_gagetio
#
check_gadgetio: tools/check_gadgetio.c $(SRC)
		$(CC) $(INC) $(LIB) $^ -o bin/$@ $(FLAGS)


#---------------------------------------------
#---------------------------------------------
#---------------------------------------------
#---------------------------------------------
sizemass_eagle: sizemass_eagle.c archive.c catalog.c stf.c simulation.c misc.c halomaker.c particle.c ramses.c structure.c hdf5routines.c hdf5sim.c
		$(CC) $(INC) $(LIB) sizemass_eagle.c archive.c catalog.c stf.c halomaker.c misc.c ramses.c simulation.c particle.c structure.c hdf5routines.c hdf5sim.c -o bin/sizemass_eagle $(FLAGS)

surface_density: test.c archive.c catalog.c misc.c stf.c simulation.c gadget.c halomaker.c particle.c ramses.c structure.c hdf5routines.c hdf5sim.c
		$(CC) $(INC) $(LIB) test.c archive.c catalog.c misc.c stf.c halomaker.c gadget.c ramses.c simulation.c particle.c structure.c hdf5routines.c hdf5sim.c -o bin/surface_density $(FLAGS)

density: density.c archive.c catalog.c misc.c stf.c simulation.c gadget.c halomaker.c particle.c ramses.c structure.c hdf5routines.c hdf5sim.c
		$(CC) $(INC) $(LIB) density.c archive.c catalog.c misc.c stf.c halomaker.c gadget.c ramses.c simulation.c particle.c structure.c hdf5routines.c hdf5sim.c -o bin/density $(FLAGS)


isgal: isgal.c archive.c catalog.c misc.c stf.c simulation.c gadget.c halomaker.c particle.c ramses.c structure.c hdf5routines.c hdf5sim.c
	$(CC) $(INC) $(LIB) isgal.c archive.c catalog.c misc.c stf.c halomaker.c gadget.c ramses.c simulation.c particle.c structure.c hdf5routines.c hdf5sim.c -o bin/isgal $(FLAGS)

calcSO: calcSO.c archive.c catalog.c misc.c stf.c simulation.c gadget.c halomaker.c particle.c ramses.c structure.c hdf5routines.c hdf5sim.c
	$(CC) $(INC) $(LIB) calcSO.c archive.c catalog.c misc.c stf.c halomaker.c gadget.c ramses.c simulation.c particle.c structure.c hdf5routines.c hdf5sim.c -o bin/calcSO $(FLAGS)


analyze_galaxy_catalog: analyze_galaxy_catalog.c
		gcc analyze_galaxy_catalog.c -lm -o bin/analyze_galaxy_catalog

ctlgMatch: ctlgmatch.c archive.c gadget.c catalog.c stf.c halomaker.c simulation.c particle.c ramses.c structure.c hdf5routines.c hdf5sim.c misc.c
		$(CC) $(INC) $(LIB) ctlgmatch.c archive.c gadget.c catalog.c stf.c  halomaker.c simulation.c particle.c ramses.c structure.c hdf5routines.c hdf5sim.c misc.c -o bin/ctlgMatch $(FLAGS)

find: find.c archive.c catalog.c stf.c gadget.c halomaker.c simulation.c particle.c ramses.c structure.c hdf5routines.c hdf5sim.c misc.c
		$(CC) $(INC) $(LIB) find.c archive.c catalog.c stf.c gadget.c halomaker.c simulation.c particle.c ramses.c structure.c hdf5routines.c hdf5sim.c misc.c -o bin/find $(FLAGS)

convert: convert.c archive.c catalog.c stf.c halomaker.c
		gcc convert.c archive.c catalog.c stf.c  halomaker.c -lm -o bin/convert

test: test.c archive.c catalog.c misc.c stf.c simulation.c gadget.c halomaker.c particle.c ramses.c structure.c hdf5routines.c hdf5sim.c
		$(CC) $(INC) $(LIB) test.c archive.c catalog.c misc.c stf.c halomaker.c gadget.c ramses.c simulation.c particle.c structure.c hdf5routines.c hdf5sim.c -o bin/test $(FLAGS)
