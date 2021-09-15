/*      Fft.c
        1-d fast Fourier transform subroutine.
        From Claerbout (1985) p. 70.
                        fft(lx, cx, signi, scale);
        Arguments:      lx      number of elements of vector to transform
                                MUST BE A POWER OF 2
                        cx      pointer to vector of complex structures
                        signi   sign of transform- +1 forward to Fourier domain
                                                   -1 inverse to real domain
       fft( npts, &samp[0], FORWARD, scale)
**/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "complex.h"

void fft1d( Complex *cx, float signi, float dt, int lx )
{
	int i, j, k, m, istep;
	float arg, df, scale;
	Complex cw, ct;
	int GetPow2( int );

	if( lx != GetPow2( lx ) )
	{
		fprintf( stderr, "\nfft: ERROR!! lx must be power of 2 lx=%d\n", lx );
		exit(-1);
	}

/*** set constants ***/
	df = 1.0/(lx*dt);
	if( signi == FORWARD ) scale = sqrt( (1.0/lx) * (dt/df) );
	if( signi == INVERSE ) scale = sqrt( (1.0/lx) * (df/dt) );

	j = 0;
	for( i = 0; i < lx; i++) 
	{
		if ( i <= j ) 
		{
			ct.re    = scale * cx[j].re;
			ct.im    = scale * cx[j].im;
			cx[j].re = scale * cx[i].re;
			cx[j].im = scale * cx[i].im;
			cx[i].re = ct.re;
			cx[i].im = ct.im;
		}
		m = lx / 2;
		while ( j > (m - 1) && (m > 1) ) 
		{
			j = j - m;
			m = m / 2;
		}
		j = j + m;
	}
	k = 1;

	do {
		istep = 2 * k;
		for ( m = 0; m < k; m++ ) 
		{
			arg = M_PI * signi * m/k;
			cw.re = cos(arg);
			cw.im = sin(arg);

			for ( i = m; i < lx; i += istep ) 
			{
				ct.re      = (cw.re * cx[i+k].re) - (cw.im * cx[i+k].im);
				ct.im      = (cw.im * cx[i+k].re) + (cw.re * cx[i+k].im);
				cx[i+k].re = cx[i].re - ct.re;
				cx[i+k].im = cx[i].im - ct.im;
				cx[i].re   = cx[i].re + ct.re;
				cx[i].im   = cx[i].im + ct.im;
			}
		}
		k = istep;
	} while (k < lx);
	return;
}

int GetPow2( int inum )
{
        int j=1;
	while( j<inum ) j=j*2;
        return j;
}
