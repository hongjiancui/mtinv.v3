#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "sacfile.h"

char progname[128];

int main(int ac, char **av)
{
	SacFile *sf;
	double mintime, maxtime, totaltime;
	int nfiles = 0;
	float dt;
	float *data;
	int rank = 1;    /** this is the sac file that is the 1st (earliest start time) ***/
	int i, j, it;
	int ndata;

	int verbose = 0; /*** for debug ***/

/*** function prototypes ***/

	void readsacfile( SacFile *sf, int verbose );
	void set_sac_minmax( Sac_Header *s, float *data );
	void WriteMyTime2STDOUT( MyTime *t );
	int cvt2time( double minepoch, double maxepoch, double myepoch, float dt );
	void remove_mean( float *, int );
	void wrtoldsac( char *FO, Sac_Header *s, float *data );

/*** start program ***/

	strcpy( progname, av[0] );

	sf = malloc(sizeof(SacFile));
	maxtime = 0;
	mintime = 9e+25;

	for( i=1; i<ac; i++ )
	{
		sf = realloc( sf, (i+1)*sizeof(SacFile) );
		
		strcpy( sf[i].filename, av[i] );

		if(verbose)	
		{
		  fprintf( stderr, "%s: opening file %d %s\n",
                        progname, i, sf[i].filename );
		}

		readsacfile( &sf[i], verbose );

		fprintf( stdout, "%s: %s: nt=%d dt=%g b=%g e=%g dur=%g \n",
			progname,
			sf[i].filename,
			sf[i].s.npts,
			sf[i].s.delta,
			sf[i].s.b,
			sf[i].s.e,
			(sf[i].s.e - sf[i].s.b) );
			
		/*
                       WriteMyTime2STDOUT( &(sf[i].beg) );
                       WriteMyTime2STDOUT( &(sf[i].end) );
		*/

	/** rank: find the sac file with minimum(earliest) start time */

		if( sf[i].beg.epoch < mintime )
		{
			mintime = sf[i].beg.epoch;
			rank = i;
		}

		if( sf[i].end.epoch > maxtime )
		{
			maxtime = sf[i].end.epoch;
		}
		
		nfiles++;

	} /*** loop over all command line arguments ***/

	if( nfiles == 0 ) 
	{
		fprintf(stderr, "%s: Fatal error no files read\n", progname );
		exit(-1);
	}

	dt = sf[1].s.delta;
	totaltime = maxtime - mintime;
	ndata = (int) (totaltime / dt);

	fprintf( stdout, "%s: nfiles=%d rank=%d total time = %lf sec / %lf hrs ndata=%d dt=%g\n",
		progname, nfiles, rank, totaltime, totaltime/3600., ndata, dt );

	data = calloc( ndata, sizeof(float) );

	for( i=1; i<=nfiles; i++ )
	{
		j = cvt2time( mintime, maxtime, sf[i].beg.epoch, dt );

		if(verbose)
		{
		  fprintf(stderr, "%d %s mintime=%lf maxtime=%lf b=%g myepoch=%lf j=%d rank=%d\n",
                        i,
                        sf[i].filename,
                        mintime,
                        maxtime,
                        sf[i].s.b,
                        sf[i].beg.epoch,
                        j,
                        rank );
		}

		/* remove_mean( &(sf[i].data[0]), sf[i].s.npts ); */

		for( it=0; it<sf[i].s.npts; it++ )
		{
			data[j] = sf[i].data[it];
			j++;
		}

	} /*** loop over each file ***/

	sf[rank].s.npts = ndata;
	sf[rank].s.b    = sf[rank].s.b;
	sf[rank].s.e    = rint( ndata * dt );

	remove_mean( &data[0], ndata );

	set_sac_minmax( &(sf[rank].s), &data[0] );

	fprintf( stdout, "%s: writting out.sac nt=%d dt=%g b=%g e=%g dur=%g min=%g max=%g mean=%g\n",
                                progname,
                                sf[rank].s.npts,
                                sf[rank].s.delta,
                                sf[rank].s.b,
                                sf[rank].s.e,
                                (sf[rank].s.e - sf[i].s.b),
				sf[rank].s.depmin,
				sf[rank].s.depmax,
				sf[rank].s.depmen );

	wrtoldsac( "out.sac", &(sf[rank].s), &data[0] );

	free(data);
}

int cvt2time( double minepoch, double maxepoch, double myepoch, float dt )
{
        int ipt=0;
        ipt = (int)(( myepoch - minepoch )/dt);
        return ipt;
}
