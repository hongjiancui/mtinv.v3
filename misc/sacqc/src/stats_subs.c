#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "hist.h"

void stats2( int nsamps, float *y, Statistics *s, Histogram *h, int verbose )
{
	float *x;
	int i, k;
	void median( float *x, int n, float *med, float *min, float *max );
        void moment( float *x, int n, float *ave, float *adev, float *sdev, float *svar, float *skew, float *curt );

        x = (float *)calloc( nsamps+1, sizeof(float) );
        for( k = 1, i = 0; i < nsamps; i++ )
        {
		if( y[i] > ( h->mode - 0.1 ) )
		{
                	x[k] = y[i];
			k++;
		}
        }
        moment( x, k, &(s->ave), &(s->adev), &(s->sdev), &(s->svar), &(s->skew), &(s->curt) );
        median( x, k, &(s->med), &(s->min), &(s->max) );
        free(x);
}

void stats( int nsamps, float *y, Statistics *s, int verbose )
{
	float *x;
	int i;
        void median( float *x, int n, float *med, float *min, float *max );
        void moment( float *x, int n, float *ave, float *adev, float *sdev, float *svar, float *skew, float *curt );

	x = (float *)calloc( nsamps+1, sizeof(float) );
	for( i = 0; i < nsamps; i++ )
	{
		x[i+1] = y[i];
	}
        moment( x, nsamps, &(s->ave), &(s->adev), &(s->sdev), &(s->svar), &(s->skew), &(s->curt) );
        median( x, nsamps, &(s->med), &(s->min), &(s->max) );
	free(x);

	if(verbose)
	{
          fprintf( stderr,
            "n=%d min=%g med=%g max=%g ave=%g adev=%g sdev=%g svar=%g skew=%g curt=%g\n",
                nsamps,
                s->min,
                s->med,
                s->max,
                s->ave,
                s->adev,
                s->sdev,
                s->svar,
                s->skew,
                s->curt);
	}
}

void median( float *x, int n, float *med, float *min, float *max )
{
	int *indx;
	void indexx( int n, float *arrin, int *indx );

	indx = (int *)malloc( (n+1) * sizeof(int));

	indexx( n, x, indx );

	if( (n%2) == 0 ) 
	{
		*med=x[indx[n/2]];
	}
	else 
	{
		*med=x[indx[(n+1)/2]];
	}
	*min = x[indx[1]];
	*max = x[indx[n]];

	free(indx);
}

void moment( float *data, int n, float *ave, float *adev, float *sdev, float *svar, float *skew, float *curt )
{
        int j;
        float s,p;

        if (n <= 1) 
	{
		fprintf(stderr, "n must be at least 2 in MOMENT\n");
		exit(-1);
	}
        s=0.0;
        for (j=1;j<=n;j++) s += data[j];
        *ave=s/n;
        *adev=(*svar)=(*skew)=(*curt)=0.0;
        for (j=1;j<=n;j++) {
                *adev += fabs(s=data[j]-(*ave));
                *svar += (p=s*s);
                *skew += (p *= s);
                *curt += (p *= s);
        }
        *adev /= n;
        *svar /= (n-1);
        *sdev=sqrt(*svar);
        if (*svar) 
	{
                *skew /= (n*(*svar)*(*sdev));
                *curt=(*curt)/(n*(*svar)*(*svar))-3.0;
        } 
	else 
	{
		fprintf(stderr, "No skew/kurtosis when variance = 0 (in MOMENT)\n");
	}
}

void indexx( int n, float *arrin, int *indx )
{
	int l,j,ir,indxt,i;
	float q;

	for (j=1;j<=n;j++) indx[j]=j;
	l=(n >> 1) + 1;
	ir=n;
	for (;;) {
		if (l > 1)
			q=arrin[(indxt=indx[--l])];
		else {
			q=arrin[(indxt=indx[ir])];
			indx[ir]=indx[1];
			if (--ir == 1) {
				indx[1]=indxt;
				return;
			}
		}
		i=l;
		j=l << 1;
		while (j <= ir) {
			if (j < ir && arrin[indx[j]] < arrin[indx[j+1]]) j++;
			if (q < arrin[indx[j]]) {
				indx[i]=indx[j];
				j += (i=j);
			}
			else j=ir+1;
		}
		indx[i]=indxt;
	}
}
