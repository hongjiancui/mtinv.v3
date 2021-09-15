#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hist.h"

char progname[128];

int main( int ac, char **av )
{
	Histogram h;
	float *data, tmp;
	int npts;
	int i;

	float xstart, xstop, binwidth;
	int verbose = 0;
	int ilog = 0;
	int retv;
	FILE *fp;
	char filename[128];
	
	void makehist( float *data, int npts, float xstart, float xstop, float binwidth, Histogram *h, int verbose );
	int setpar(int ac, char **av);
	int getpar(), mstpar();
	void endpar();
	
	strcpy( progname, av[0] );

	setpar( ac, av );
	mstpar( "f",        "s", filename );
	mstpar( "xstart",   "f", &xstart );
	mstpar( "xstop",    "f", &xstop );
	mstpar( "binwidth", "f", &binwidth );
	getpar( "log",      "b", &ilog );
	getpar( "verbose",  "b", &verbose );
	endpar();

	if( (fp = fopen( filename, "r" )) == NULL )
	{
		fprintf( stderr, "%s: error opening file %s\n",
			progname, filename );
		exit(-1);
	}
	fprintf( stderr, "%s: opening file %s\n",  progname, filename );

	data = (float *)malloc( 1000 * sizeof(float) );
	npts = 0;
	do {
		data = (float *)realloc( data, (npts+1)*sizeof(float) );
		retv = fscanf( fp, "%f", &data[npts] );
		if( ilog )
		{
			if( data[npts] <= 0  ) 
				tmp = log10( 1.0 );
			else
				tmp = log10( data[npts] );
			data[npts] = tmp;
		}
		npts++;
	} while( retv == 1 );
	npts--;

	fprintf( stderr, "%s: read %d points from %s\n", progname, npts, filename );

	makehist( data, npts, xstart, xstop, binwidth, &h, verbose );

	for( i = 0; i < h.nbin; i++ )
	{
		fprintf( stdout, "%g %g\n", 
			h.b[i].start, 
			h.b[i].counts );
	}
	free(data);
	exit(0);
}
