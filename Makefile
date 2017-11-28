#CC     = gcc
#CFLAGS = -I

analyze_galaxy_catalog: analyze_galaxy_catalog.c
	gcc analyze_galaxy_catalog.c -lm -o bin/analyze_galaxy_catalog

cross_match_catalogs: cross_match_catalogs.c
	gcc cross_match_catalogs.c -lm -o bin/cross_match_catalogs
