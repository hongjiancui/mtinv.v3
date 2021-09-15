#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "complex.h"
#include "sac.h"

int main(int ac, char **av)
{
	int i,j,k, npts, n2, padfac=0;
	float scale, dt, df, nyquist;
	float *data, *am, *ph;
	float *power;
	Complex *samp;
	struct sac_header s;
	char FileName[128];
	FILE *fp;
	void fft1d();
	int GetPow2();

	setpar(ac,av);
	mstpar("if", "s", FileName);
	getpar("pad", "d", &padfac);
	endpar();

/* read header and data from sac file */
	fp = fopen(FileName, "rb");
	fread(&s, sizeof( struct sac_header ), 1, fp);
	dt   = s.delta;
	data = (float *)calloc(s.npts, sizeof(float));
	fread(&data[0], s.npts*sizeof(float), 1, fp);
	fclose(fp);
	fprintf(stderr, "fi=%s n=%d dt=%g\n",FileName, s.npts, s.delta);

/* set some constants */
	npts    = get_pow_2(s.npts);
	nyquist = 0.5 / dt;
	df      = 2 * nyquist / npts;
	scale   = sqrt( 1.0/(float)npts * (dt/df) );

 	fprintf(stderr, "npts=%d df=%g scale=%g\n", npts, df, scale);

/* allocate space for complex array and transfer over data points */
        samp = (Complex *)calloc(npts, sizeof(Complex));
        for(i=0; i<s.npts; i++) {
                samp[i].re = data[i];
                samp[i].im = 0.0;
        }
	for(i=s.npts; i<npts; i++) {
		samp[i].re = 0.0;
		samp[i].im = 0.0;
	}

/* forward fft and then put phase and amplitude spectra into seperate output */
	fft1d( npts, &samp[0], FORWARD, scale);

	am = (float *)calloc(npts, sizeof(float));
	ph = (float *)calloc(npts, sizeof(float));
	n2 = (int)(npts/2);

	for( i=0; i<n2; i++)
		samp[i].im = -1 * samp[i].im ;

	for( i=n2; i<npts; i++)
		samp[i].im = 1 * samp[i].im ;

	for( i=0; i<n2; i++)
		samp[i].re = -1 * samp[i].re ;

	for( i=n2; i<npts; i++)
		samp[i].re = 1 * samp[i].re ;

	scale   = sqrt( 1.0/(float)npts * (df/dt) );
	fft1d( npts, &samp[0], INVERSE, scale);

	for( i=0; i<npts; i++) 
		ph[i] = samp[i].im;

	s.b      = 0.0;
	s.delta  = dt;
	s.npts   = s.npts;

        s.nzyear = -12345.;
        s.nzjday = -12345.;
        s.nzhour = -12345.;
        s.nzmin  = -12345.;
        s.nzmsec = -12345.;
        s.t0     = -12345.;
        s.t1     = -12345.;
        s.t2     = -12345.;
        s.a      = -12345.;

	fp = fopen("hil.im", "wb");
        fwrite(&s, sizeof( struct sac_header), 1, fp);
        fwrite(&ph[0], npts*sizeof(float), 1, fp);
        fclose(fp);

	return;
}
