/**************************************************************************************/
/*** interpolate_wiggins.c - Interpolates evenly or unevenly spaced data.           ***/
/***    Wiggins, 1976, BSSA, 66, p.2077.                                            ***/
/*** INPUT                                                                          ***/
/***  y       is the input array                                                    ***/
/***  npts    is the number of points                                               ***/
/***  delta   is the sampling rate                                                  ***/
/***  b       is the start time                                                     ***/
/***  new_nt  new number of points                                                  ***/
/***  new_dt  new sampling rate                                                     ***/
/*** OUTPUT                                                                         ***/
/***  z       new array with new nt and dt                                          ***/
/**************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

char progname[128];

void interpolate_wiggins2( float *data, int npts, float delta, float b, int new_nt, float new_dt, int verbose )
{
	int i;
	float xnew;
	float eps=0.000001;
	float *z, *x, *ytmp;

	float wigint( float *, float *, int, float, float, float );

	if(verbose) 
	{
	  fprintf( stdout, "%s: interpolate_wiggins2(): from nt=%d dt=%g to nt=%d dt=%g\n",
		progname, npts, delta, new_nt, new_dt );
	}

/*** allocate memory for temporary spaces ***/
	x    = (float *)calloc(npts+1, sizeof(float) );
	ytmp = (float *)calloc(npts+1, sizeof(float) );
	z    = (float *)calloc(new_nt+1, sizeof(float) );

/*** do the interpolation ***/
	for( i=1; i<= npts; i++ )
	{
		x[i] = b + (delta*(i-1));
		ytmp[i] = data[i-1];
	}
	for( i=1; i <= new_nt; i++ )
	{
		xnew = b + (new_dt*(i-1));
		z[i-1] = wigint( x, ytmp, npts, delta, eps, xnew );
	}

/*** resize and replace input Sac_File data vector ***/
	for( i=0; i<new_nt; i++ )
		data[i] = z[i];

/*** free memory and return ***/
	free(z);
	free(ytmp);
	free(x);
}

void interpolate_wiggins( float *y, int npts, float delta, float b, float *z, int new_nt, float new_dt )
{
	int i;
	float xnew;
	float eps=0.0001;
	float *x, *ytmp;
	float wigint( float *, float *, int, float, float, float );

	// fprintf( stderr, "interpolate from nt=%d dt=%g to nt=%d dt=%g\n",
	//	npts, delta, new_nt, new_dt );

	x    = (float *)calloc(npts+1, sizeof(float) );
	ytmp = (float *)calloc(npts+1, sizeof(float) );

	for( i=1; i<= npts; i++ )
	{
		x[i] = b + (delta*(i-1));
		ytmp[i] = y[i-1];
	}

	for( i=1; i <= new_nt; i++ )
	{
		xnew = b + (new_dt*(i-1));
		z[i-1] = wigint( x, ytmp, npts, delta, eps, xnew );
		// fprintf( stdout, "%d %g %g\n", i, xnew, z[i-1] );
	}

	free(ytmp);
	free(x);
}

/**************************************************************************************/
/*** wigint.c                                                                       ***/
/*** INPUT                                                                          ***/
/*** X:       X array if unevenly spaced, first x if evenly spaced.                 ***/
/*** Y:       Y array.                                                              ***/
/*** NPTS:    Length of (X and) Y arrays.                                           ***/
/*** DX:      Set to 0.0 if unevenly spaced, to sampling interval if evenly spaced. ***/
/*** EPS:     Interpolation factor.                                                 ***/
/*** T:       Time value to interpolate to.                                         ***/
/*** OUTPUT                                                                         ***/
/*** F:       Interpolated y value.                                                 ***/
/**************************************************************************************/

float wigint( float *x, float *y, int npts, float dx, float eps, float t )
{
	int i, j, n1;
	float fout;
	float dxj1, dxj, a, h, hs, hc, dxjs, dxj1s, dy;
	float am, amd, amu;
	float dxd, dxu, dyd, dyu, wd, w, wu, sp, sp1;
	float t1, t2, t3, t4;
	float epsi;
	float amax1( float, float );

	epsi = 0.0001;
	if( eps > 0 ) epsi = eps;

	if( dx == 0 ) 
	{
		for( j =1; j <= npts; j++ )
		{
			a = x[j] - t;
			if( a > 0 ) break;
		}
		j = j - 1;
		dxj = t - x[j];
		if( dxj == 0 ) return y[j];
		h = x[j+1] - x[j];
		dxj1 = t - x[j+1];
	}	
	else
	{
		j = (t-x[1])/dx;
		dxj = t - x[1] - j * dx;
		j = j + 1;
		if( dxj == 0 ) return y[j];
		h = dx;
		dxj1 = dxj - h;
		dxd = h;
		dxu = h;
	}

	hs = h  * h;
	hc = hs * h;
	dxjs  = dxj  * dxj;
	dxj1s = dxj1 * dxj1;
	dy = y[j+1] - y[j];
	am = dy/h;
	amd = am;
	amu = am;

	if( j != 1 ) 
	{
		if( dx != 0 ) 
			dxd = dx;
		else
			dxd = x[j] - x[j-1];

		dyd = y[j] - y[j-1];
		amd = dyd/dxd;
	}
	n1 = j + 1;
	if( n1 != npts )
	{
		if( dx != 0 )
			dxu = dx;
		else
			dxu = x[j+2] - x[j+1];
		dyu = y[j+2] - y[+1];
		amu = dyu / dxu;
	}
	
	wd = 1.0 / amax1( fabs(amd), epsi );
	w  = 1.0 / amax1( fabs(am),  epsi );
	wu = 1.0 / amax1( fabs(amu), epsi );
	sp  = ( wd * amd + w  * am  )/( wd + w );
	sp1 = ( w  * am  + wu * amu )/( w + wu );
	t1 = y[j]   * ((dxj1s/hs) + (2.0 * dxj  * dxj1s/ hc) );
	t2 = y[j+1] * ((dxjs /hs) - (2.0 * dxj1 * dxjs / hc) );
	t3 = sp * dxj * dxj1s/hs;
	t4 = sp1 * dxjs * dxj1/hs;
	fout  = t1 + t2 + t3 + t4;

	return fout;
}

float amax1( float x, float y )
{
	if( x >= y ) 
		return x;
	else
		return y;
}
