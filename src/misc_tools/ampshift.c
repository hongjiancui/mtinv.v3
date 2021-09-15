#include <stdio.h>

void ampshift( float *x, int n, int verbose )
{
	int i;
	float x0,xtmp;
	x0 = x[0];
	if( verbose ) printf( "ampshift(): n=%d x0=%g\n", n, x0 );

	for( i=0; i<n; i++ )
	{
		xtmp = x[i] - x0;
		x[i] = xtmp;
	}
}
