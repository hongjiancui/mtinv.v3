#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sac.h"

char progname[256];

int main( int ac, char **av)
{
	Sac_Header s;
	int i, nt, xreset=FALSE;
	float *x, *y;
	FILE *fp;
	int norm = 0, xdeg = 0;
	float ymax, scale = 1.0, ydist;
	int setpar(int, char **);
	int getpar();
	void endpar();
	int verbose = 0;

	sprintf( progname, "%s", av[0] );

	setpar( ac, av );
	getpar( "xreset", "b", &xreset);
	getpar( "norm",   "b", &norm );
	getpar( "sc",     "f", &scale );
	getpar( "xdeg",   "b", &xdeg );
	getpar( "verbose", "b", &verbose );
	endpar();

	if(verbose)
	  fprintf( stderr, "%s: xreset=%d norm=%d sc=%g xdeg=%d\n",
		progname, xreset, norm, scale, xdeg );

	fread(&s, sizeof(Sac_Header), 1, stdin);
	x = (float *)calloc(s.npts, sizeof(float));
        y = (float *)calloc(s.npts, sizeof(float));

        fread(&y[0], s.npts*sizeof(float), 1, stdin);

	if(verbose) fprintf( stderr, "%s: npts=%d delta=%g b=%g sta=%s\n",
		progname, s.npts, s.delta, s.b, s.kstnm );

	/* if( s.iftype == ITIME ) s.b=0; */
	if( xreset == TRUE ) s.b=0;

	if( norm )
	{
		ymax = fabs(y[0]);

		for( i=1; i<s.npts; i++ )
		{
			if( fabs(y[i]) > ymax ) ymax = fabs(y[i]);
		}
	}

	fprintf( stdout, ">  %s %g\n", s.kstnm, s.dist );

        for( i=0; i<s.npts; i++)
	{
		x[i]=s.b + (float)i*s.delta;

	/*** prevent division by zero ***/
		if( x[i] == 0 ) x[i]=1.0e-9;
		if( y[i] == 0 ) y[i]=1.0e-9;
	
		if( xdeg )
		{
			ydist = s.dist;
		}
		else
		{
			ydist = 0;
		}

		if( norm )
			fprintf( stdout, "%g %g\n", x[i], scale*(y[i]/ymax) + ydist );
		else
                	fprintf( stdout, "%g %g\n", x[i], scale*y[i] + ydist );
	}

	if(verbose) fprintf( stderr, "%s: exiting \n", progname );
}
