#CC     = gcc
#CFLAGS = -I

analyze_galaxy_catalog: analyze_galaxy_catalog.c
	gcc analyze_galaxy_catalog.c -lm -o bin/analyze_galaxy_catalog

ctlgMatch: ctlgmatch.c archive.c catalog.c stf.c halomaker.c simulation.c particle.c ramses.c
	gcc ctlgmatch.c archive.c catalog.c stf.c  halomaker.c simulation.c particle.c ramses.c -lm -o bin/ctlgMatch

convert: convert.c archive.c catalog.c stf.c halomaker.c
	gcc convert.c archive.c catalog.c stf.c  halomaker.c -lm -o bin/convert
