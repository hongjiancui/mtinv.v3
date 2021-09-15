#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hist.h"

char progname[128];
	
void makehist( float *data, int npts, float xstart, float xstop, float binwidth, Histogram *h, int verbose )
{
	int i, j;
	int *idx;
	float xmin, xmax;
	FILE *fp;

/*** begin program ***/

	idx = calloc(npts, sizeof(int));

	h->nbin = (int)( ( xstop - xstart )/binwidth );

	xmin = +1.0E+09;
	xmax = -1.0E+09;

	for( j = 0; j < npts; j++ )
	{
		if( data[j] < xmin ) xmin = data[j];
		if( data[j] > xmax ) xmax = data[j];
	}

	if(verbose)
	{
	  fprintf( stderr, "%s: nbin = %d xmin = %g xmax = %g \n",
		progname,
		h->nbin,
		xmin,
		xmax );
	}

	h->b = malloc( h->nbin * sizeof(Bin) );

	for( i = 0; i < h->nbin; i++ )
	{
		h->b[i].start  = xstart + binwidth * (float)i;
		h->b[i].stop   = xstart + binwidth * (float)(i+1);
		h->b[i].center = h->b[i].start + 0.5 * ( h->b[i].stop - h->b[i].start );
		h->b[i].counts = 0;
	}

/*** debug output ***/
/***
	fp = fopen("tmp.out", "w" );
	for( j = 0; j < npts; j++ )
	{
		fprintf( fp, "%g\n", data[j] );
	}
	fclose(fp);
***/

	for( j = 0; j < npts; j++ )
	{
		idx[j] = (int)floor( (data[j]/(xstop-xstart)) * (float)h->nbin );

		if( idx[j] >= 0 && idx[j] < h->nbin )
		{
			h->b[idx[j]].counts++;
		}
		else
		{
			if(verbose)
			{
			  fprintf( stderr, 
			    "%s: warning j = %d point out of bounds data = %g idx = %d\n",
				progname, j, data[j], idx[j] );
			}
		}
	}

	free(idx);
}

void hist2mode( Histogram *h, int verbose )
{
	int i;
	float mode;
	double max;

	mode = h->b[0].center;
	max = 1.0E-13;
	for( i = 0; i < h->nbin; i++ )
	{
		if( h->b[i].counts > max )
		{
			max  = h->b[i].counts;
			mode = h->b[i].center;
		}
	}
	h->mode = mode;
}

void hist2mean( Histogram *h, int verbose )
{
	int i;
	double nsum = 0;
	float sum = 0;

	for( i = 0; i < h->nbin; i++ )
	{
		sum += h->b[i].center * (float)h->b[i].counts;
		nsum += h->b[i].counts;
	}
	h->mean = sum/(float)nsum;
}

void hist2median( Histogram *h, int verbose )
{
	int i;
	double nsum = 0;
	double dmedian;
	float median;

	for( i = 0; i < h->nbin; i++ )
	{
		nsum += h->b[i].counts;
	}
	dmedian = ( nsum + 1 )/2;

	nsum = 0;
	for( i = 0; i < h->nbin; i++ )
	{
		nsum += h->b[i].counts;
		if( nsum <= dmedian )
		{
			median = h->b[i].center;
		}
	}
	h->median = median;
}

void hist_normalize_percentage( Histogram *h, int verbose )
{
	int i;
	double ntotal = 0;

	for( i = 0; i < h->nbin; i++ )
	{
		ntotal += h->b[i].counts;
	}

	for( i = 0; i < h->nbin; i++ )
	{
		if( ntotal == 0 )
		{
			h->b[i].percent = 0;
		}
		else
		{
			h->b[i].percent = 100 * ( h->b[i].counts / ntotal );
		}
	}
	
/*** set some max non-zero values ***/

	h->xmax_percent =  h->b[0].percent;

	for( i = 1; i < h->nbin; i++ )
	{
		if( h->b[i].percent > 0 )
		{
			h->xmax         = h->b[i].center;
			h->xmax_percent = h->b[i].percent;
		}
	}
}
