#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sac.h"

int main( int ac, char **av )
{
	char SacFileName[256];
	char sacoutfile[256];

/*** original sac header ***/
	Sac_Header s;
	float *data;

/*** function prototypes ***/
	void envelope( float *, int, float );
	float *readsac( Sac_Header *, char * );
	void set_sac_minmax( Sac_Header *, float * );
	void wrtoldsac( char *, Sac_Header *, float * );
        int setpar( int, char ** );
        int mstpar(), getpar();
        void endpar();

/*** get command line args ***/
	setpar(ac,av);
        mstpar( "sacf", "s", SacFileName );
        endpar();

/*** read the sac data and allocate memory for output ***/
	data = calloc( 10000, sizeof(float) );
	data = readsac( &s, SacFileName );

	/* env = calloc( s.npts+1, sizeof(float) ); */

	fprintf( stdout, "%s: nt=%d dt=%g b=%g e=%g min=%g max=%g amp0=%g\n",
	    SacFileName, s.npts, s.delta, s.b, s.e, s.depmin, s.depmax, data[0] );

/*** calculate the envelope ***/

	envelope( &data[0], s.npts, s.delta );

/*** output ***/
	set_sac_minmax( &s, &data[0] );

	sprintf( sacoutfile, "env.sac" ); 

	wrtoldsac( sacoutfile, &s, &data[0] );

	return 0;
}
