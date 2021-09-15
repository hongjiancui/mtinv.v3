#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../include/sac.h"

int main( int ac, char **av )
{
	float *data;
	Sac_Header s;
	char FileName[128];
	int i;
	float new_dt;
	int new_nt;
	float *z;

	void wrtoldsac( char *, Sac_Header *, float * );
	float *readsac( Sac_Header *, char * );
	void interpolate_wiggins( float *, int, float, float, float *, int, float );
	void set_sac_minmax( Sac_Header *, float * );

	setpar(ac,av);
	mstpar("f",  "s", FileName);
	mstpar("nt", "d", &new_nt);
	mstpar("dt", "f", &new_dt);
	endpar();

	data = (float *)readsac( &s, FileName );
	fprintf( stdout, "%s: nt=%d dt=%g\n", FileName, s.npts, s.delta );

	z = (float *)calloc( new_nt, sizeof(float) );

	interpolate_wiggins( data, s.npts, s.delta, s.b, z, new_nt, new_dt );

	s.npts = new_nt;
	s.delta = new_dt;
	
	set_sac_minmax( &s, data );

	wrtoldsac( "xint.sac", &s, z );

	return 0;
}
