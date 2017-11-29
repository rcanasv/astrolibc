#CC     = gcc
#CFLAGS = -I

analyze_galaxy_catalog: analyze_galaxy_catalog.c
	gcc analyze_galaxy_catalog.c -lm -o bin/analyze_galaxy_catalog

ctlgMatch: ctlgmatch.c archive.c catalog.c stf.c
	gcc ctlgmatch.c archive.c catalog.c stf.c -lm -o bin/ctlgMatch
