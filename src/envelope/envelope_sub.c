#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "complex.h"
#include "sac.h"

void envelope( float *y, int npts, float dt )
{
        Complex *z, *v;
        int i, nt, n2;
        void fft1d( Complex *, float, float, int );
	int GetPow2( int );

        nt = GetPow2( npts );
        n2 = (int)(nt/2);

        z = (Complex *)calloc(nt, sizeof(Complex));
        v = (Complex *)calloc(nt, sizeof(Complex));

        for(i=0; i<npts; i++) 
	{ 
                z[i].re = y[i];
                z[i].im = 0.0;
                v[i].re = y[i];
                v[i].im = 0.0;
        }
        for(i=npts; i<nt; i++) 
	{
                z[i].re = 0.0;
                z[i].im = 0.0;
                v[i].re = 0.0;
                v[i].im = 0.0;
        }

        fft1d( &z[0], FORWARD, dt, nt );

        for( i=0; i<n2; i++) z[i].im = -1 * z[i].im;
        for( i=0; i<n2; i++) z[i].re = -1 * z[i].re;

        fft1d( &z[0], INVERSE, dt, nt);

        for( i=0; i<npts; i++)
	{
		y[i]=sqrt(v[i].re*v[i].re + z[i].im*z[i].im);
	}

        free(v);
	free(z);
}
