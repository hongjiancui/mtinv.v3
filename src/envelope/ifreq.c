#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "complex.h"
#include "sac.h"

int main(int ac, char **av)
{
	int npts;
	float *y, *u, *a, *b, dt;

        Complex *z, *v;
        float sc, df;
        int i, j, k, nt, n2;

	struct sac_header s;
	char FileName[128];
	FILE *fp;

        void fft1d();
        int GetPow2();

	setpar(ac,av);
	mstpar("if", "s", FileName);
	endpar();

/* read header and data from sac file */
        fp = fopen(FileName, "rb");
        fread(&s, sizeof( struct sac_header ), 1, fp);
	npts = s.npts;
        dt   = s.delta;
        y    = (float *)calloc(s.npts, sizeof(float));
	fread(&y[0], npts*sizeof(float), 1, fp);
	fclose(fp);
	fprintf(stderr, "file=%s npts=%d dt=%g\n", FileName, npts, dt);

        nt = GetPow2(npts);
        df = 1./(nt*dt);
        n2 = (int)(nt/2);
        sc = sqrt( 1.0/(float)nt * (dt/df) );

/**** u - envelope ****/
/**** v - in phase ****/
/**** z - out of phase ****/

        z = (Complex *)calloc(nt, sizeof(Complex));
        v = (Complex *)calloc(nt, sizeof(Complex));
	u = (float *)calloc(nt, sizeof(float));
	a = (float *)calloc(nt, sizeof(float));
	b = (float *)calloc(nt, sizeof(float));

/******** pad with zeros ***********/
        for(i=0; i<npts; i++) { 
                z[i].re = y[i];
                z[i].im = 0.0;
                v[i].re = y[i];
                v[i].im = 0.0;
        }

        for(i=npts; i<nt; i++) {
                z[i].re = 0.0;
                z[i].im = 0.0;
                v[i].re = 0.0;
                v[i].im = 0.0;
        }

/**** forward fourier transform ******/
        fft1d( nt, &z[0], FORWARD, sc);

/****** multiply by step function (signum function) **********/
        for( i=0; i<n2; i++) z[i].im = -1 * z[i].im;
        for( i=0; i<n2; i++) z[i].re = -1 * z[i].re;

/********* inverse fft ********/
/*** the hilbert transform is in the imaginary part ***/
        sc = 2.0 * sqrt( 1.0/(float)nt * (df/dt) );
        fft1d( nt, &z[0], INVERSE, sc);

/*** demominator is envelope function   = g(t) * conj(g(t))  *****/
	for(i=0;i<nt;i++) a[i]= sqrt(v[i].re*v[i].re + z[i].im*z[i].im);

/*** time domain differentiation 2nd order *****/
	for( i=0; i<nt-1; i++) b[i] = (z[i+1].im-z[i].im)/dt;
	b[nt-1]=0;

/****** derivative of imaginary part of g **********/
/*
	for(i=0;i<nt;i++)u[i]=z[i].im;
	mulomega( npts, dt, u, b );
*/

/***** instantaneous frequency = conj(g)*dg/dt / g*conj(g) ****/
        for( i=0; i<nt; i++) u[i] = sqrt(v[i].re*v[i].re + b[i]*b[i] )/a[i];
	wrtsac("ifreq.sac", dt, nt, &u[0], 0.);

	smooth( nt, dt, u, 50, 2 );

	wrtsac("ifreq.sac.sm", dt, nt, &u[0], 0.);

        free(v);
	free(z);
	free(u);
	free(y);
        return;
}
