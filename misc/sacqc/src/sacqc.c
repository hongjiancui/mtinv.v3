#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <sys/stat.h>
#include <math.h>

#include "sacfile.h"
#include "hist.h"

char progname[128];

int main( int ac, char **av )
{
	SacFile *sf;
	int nfiles, i, j, k;
	float *data, *diff, tmp;
	double dtmp;
	int npts;
	int retv;
	Histogram h;
	Statistics st;
	char pathname[128];
	FILE *fp;

        char kevid[32];
	char author[48];
        char *username;

/*** defaults ***/
	long evid = -1;
	int verbose = 0;
	int idiff = 0;
	int ixlog = 1;
	int iylog = 1;
	float xstart = 0;
	float xstop = 10;
	float binwidth = 0.1;
	int percentage = 1;
	int dboutput = 1;
	int SACwriteback = 1;
	int igmt5 = 1;
	int help = 0;

/***************************************************************/
/*** functional prototypes ***/
/***************************************************************/
	void Usage(void);

	int readsacfile( SacFile *sf, int verbose );
	void remove_mean( float *data, int npts );
	void rtrend( float x0, float dx, float *y, int n, int verbose );

	int setpar( int ac, char **av ), mstpar(), getpar();
	void endpar();

	void makehist( float *data, int npts, float xstart, float xstop, 
			float binwidth, Histogram *h, int verbose );

	void hist2mode( Histogram *h, int verbose );
	void hist2median( Histogram *h, int verbose );
	void hist2mean( Histogram *h, int verbose );

	void hist_check_qual( Histogram *h, int verbose );
	void hist_normalize_percentage( Histogram *h, int verbose );
	void write_bins_percentage( Histogram *h, SacFile *sf, int idiff );
	void write_bins( Histogram *h, SacFile *sf, int idiff );

	void stats( int nsamps, float *x, Statistics *st, int verbose );
	void stats2( int nsamp, float *x, Statistics *st, Histogram *h, int verbose );

	void dbwrite( Histogram *h, SacFile *sf, int idiff, long evid, char *author );
	void shorten_path( char *pathname, char *progname );
	void create_db_table();
	void create_GMT5x_script();
	void create_GMT5x_script2( Histogram *h, SacFile *sf, int idiff, Statistics *st, int verbose );

/***************************************************************/
/*** start program ***/
/***************************************************************/

	strcpy( pathname, av[0] );
	shorten_path( pathname, progname );

/***************************************************************/
/*** default values ***/
/*** noise range and realistic range of possible amplitudes ***/
/***************************************************************/

	h.maxstat  = 5;		/*** the log10 of the mode, mean or median cannot be above this value ***/
	h.minstat  = 0.5;	/*** the log10 of the mode, mean or median cannot be below this value ***/
	h.max_x    = 9;		/*** the max amp cannot be above this value ***/
	h.max_xper = 2;		/*** the max amp bin cannot contain more than 2% of the total amps ***/

	sprintf( h.reason, "-" );
	h.defmask = 'N';
	strcpy( author, "\0" );

/***************************************************************/
/*** get command line arguments ***/
/***************************************************************/

	setpar( ac, av );
	getpar( "xstart",   "f", &xstart );
        getpar( "xstop",    "f", &xstop );
        getpar( "binwidth", "f", &binwidth );
	getpar( "diff",     "b", &idiff );
	getpar( "xlog",     "b", &ixlog );
        getpar( "ylog",     "b", &iylog );
	getpar( "percent",  "b", &percentage );
	getpar( "maxstat",  "f", &(h.maxstat) );
	getpar( "max_x",    "f", &(h.max_x) );
	getpar( "max_xper", "f", &(h.max_xper) );
	getpar( "evid",     "s", kevid );
	getpar( "db",       "b", &dboutput );
	getpar( "auth",     "s", &author );
	getpar( "writeback", "b", &SACwriteback );
	getpar( "verbose",  "b", &verbose );
	getpar( "gmt5",     "b", &igmt5 );
	getpar( "help",	    "b", &help );
	endpar();

	if(help) Usage();

/* 12 and 16 bit digitizers */
/*
	fprintf( stderr, "0.5 * 2^12 = %e\n", 0.5 * pow(2,12) );
	fprintf( stderr, "0.5 * 2^16 = %e\n", 0.5 * pow(2,16) );
        fprintf( stderr, "0.5 * 2^32 = %e\n", 0.5 * pow(2,32) );
*/
	evid = atol( kevid );

	if( strcmp( author, "\0" ) == 0 )
	{
		username = getenv( "USER" );
		sprintf( author, "%s", username );
	}

/*************************************/
/*** read the sac files            ***/
/*************************************/

	sf = malloc( 3*sizeof(SacFile) );
	nfiles = 0;
	for( i = 1; i < ac; i++ )
	{
		sf = realloc( sf, (nfiles+1)*sizeof(SacFile) );

		strcpy( sf[nfiles].filename, av[i] );

		if(verbose)
		{
		  fprintf( stderr, "%s: opening file i=%d nfiles=%d %s\n",
			progname, i, nfiles, sf[nfiles].filename );
		}

		if( readsacfile( &sf[nfiles], verbose ) < 0 )
		{	
			if(verbose)
			{
			  fprintf( stderr, "%s: could not read sac file %s, skipping\n", 
				progname, sf[nfiles].filename );
			}
			continue;
		}
		else
		{
			nfiles++;
		}
	}

	if( nfiles == 0 )
	{
		fprintf( stderr, "%s: No files read\n", progname );
		exit(0);
	}

	for( i = 0; i < nfiles; i++ )
	{
		fprintf( stderr, 
	"%s: i=%d net=%s sta=%s loc=%s chan=%s npts=%d delta=%g min=%g max=%g mean=%g\n",
			progname,
			i,
			sf[i].s.knetwk,
			sf[i].s.kstnm,
			sf[i].s.khole,
			sf[i].s.kcmpnm,
			sf[i].s.npts,
			sf[i].s.delta,
			sf[i].s.depmin,
			sf[i].s.depmax,
			sf[i].s.depmen );
	}

/*************************************/
/*** compute the difference vector ***/
/*************************************/

	data = malloc( sizeof(float) );
	npts = 0;
	k = 0;
	for( i = 0; i < nfiles; i++ )
	{
		remove_mean( &(sf[i].data[0]), sf[i].s.npts );

		rtrend( 0., sf[i].s.delta, &(sf[i].data[0]), sf[i].s.npts, verbose );
		
		npts += sf[i].s.npts;
		data = realloc( data, (npts+1) * sizeof(float) );

		for( j = 0; j < sf[i].s.npts; j++ )
		{
			data[k] = sf[i].data[j];
			k++;
		}
		if(verbose)
		{
		  fprintf( stderr, "i=%d sf[i].s.npts=%d npts=%d k=%d\n", 
			i, sf[i].s.npts, npts, k );
		}
	}
	npts = k;

/*********************************************************/
/*** take the difference or report absolute Amp values ***/
/*********************************************************/

	diff = calloc( npts+1, sizeof(float) );

	if( idiff )
	{
		diff[0] = 1.0E-09;
		for( k = 1; k < npts; k++ )
		{
			diff[k] = fabs( data[k] - data[k-1] ); 
		}
	}
	else
	{
		for( k = 0; k < npts; k++ )
		{
			diff[k] = fabs( data[k] );
		}
	}

/***************************************************************/
/*** if ixlog flag set then take the logs of the diff vector ***/
/***************************************************************/

	if( ixlog )
	{
		for( k = 0; k < npts; k++ )
		{
			if( diff[k] > 0 )
			{
				tmp = log10( diff[k] );
				diff[k] = tmp;
			}
			else
			{
				diff[k] = 0.1;
			}
		}
	}

/***************************************************************/
/*** make the histograms                                     ***/
/***************************************************************/

	/* stats( npts, diff, &st, verbose ); */

	if(verbose)
	{
	  fprintf( stderr,
            "%s: stats(): n=%d min=%g med=%g max=%g ave=%g adev=%g sdev=%g svar=%g skew=%g curt=%g\n",
		progname,
                npts,
                st.min,
                st.med,
                st.max,
                st.ave,
                st.adev,
                st.sdev,
                st.svar,
                st.skew,
                st.curt );
	  fflush(stderr);
	}

	makehist( diff, npts, xstart, xstop, binwidth, &h, verbose );

	stats2( npts, diff, &st, &h, verbose );

	if(verbose)
	{
	  fprintf( stderr,
            "%s: stats2(): n=%d min=%g med=%g max=%g ave=%g adev=%g sdev=%g svar=%g skew=%g curt=%g\n",
		progname, 
                npts,
                st.min,
                st.med,
                st.max,
                st.ave,
                st.adev,
                st.sdev,
                st.svar,
                st.skew,
                st.curt );
	  fflush(stderr);
	}

/***************************************************************/
/**** make counts into percentages ***/
/***************************************************************/

	if( verbose )
	{
		fprintf( stderr, "%s: calling hist_normalize_percentage(): \n", progname );
		fflush(stderr);
	}

	hist_normalize_percentage( &h, verbose );

/***************************************************************/
/*** if iylog                                                ***/
/***************************************************************/

	if( iylog )
	{	
		for( i = 0; i < h.nbin; i++ )
		{
			if( h.b[i].counts > 0 )
			{
				dtmp = log10( h.b[i].counts );
				h.b[i].counts = dtmp;
			}
			else
			{
				h.b[i].counts = 0;
			}
		}
	}

/***************************************************************/
/*** compute mode, medan, median                             ***/
/***************************************************************/

	if( verbose )
	{
		fprintf( stderr, "%s: calling hist2mode(): \n", progname );
		fflush(stderr);
	}

	hist2mode( &h, verbose );
	hist2mean( &h, verbose );
	hist2median( &h, verbose );

/***************************************************************/
/*** check histogram and make a decision about data quality ***/
/***************************************************************/

	if( verbose )
	{
		fprintf( stderr, "%s: calling hist_check_qual(): \n", progname );
		fflush(stderr);
	}

	hist_check_qual( &h, verbose );

/***************************************************************/
/*** make histogram plot using psxy                          ***/
/***************************************************************/
	
	if( verbose ) fprintf( stderr, "%s: percentage=%d\n", progname, percentage );

	if(percentage)
	{
		if( verbose )
		{
		  fprintf( stderr,
			"%s:  calling write_bins_percentage() and create_GMT5x_script2():\n",
                        progname );
			fflush(stderr);
		}

		write_bins_percentage( &h, &(sf[0]), idiff );

		if(igmt5)
		{	
		  create_GMT5x_script2( &h, &(sf[0]), idiff, &st, verbose );
		}
	}
	else
	{
		write_bins( &h, &(sf[0]), idiff );
	}

	if(dboutput)
	{
		if( verbose )
		{
			fprintf( stderr, "%s: calling create_db_table() and dbwrite(): \n", progname );
			fflush(stderr);
		}

		create_db_table();
		dbwrite( &h, &(sf[0]), idiff, evid, author );
	}

	if( verbose )
	{
		fprintf( stderr, "%s: SACwriteback=%d\n", progname, SACwriteback );
		fflush(stderr);
	}

	if( SACwriteback )
	{
		if( h.defmask == 'Y' )
		{
			for( i = 0; i < nfiles; i++ )
			{
				sf[i].s.user2 = 999;
				sprintf( sf[i].s.kuser2, "SACQC" );

				if(verbose)
				{
				  fprintf( stderr, 
				    "%s: overwritting file=%s SAC header with kuser2=%s user2=%g\n",
					progname,
					sf[i].filename,
					sf[i].s.kuser2,
					sf[i].s.user2 );
				  fflush(stderr);
				}

				fp = fopen( sf[i].filename, "wb" );
				/* fseek( fp, (long)0, SEEK_SET ); */
				fwrite( &(sf[i].s), sizeof(Sac_Header), 1, fp );
				fwrite( &(sf[i].data[0]), (sf[i].s.npts * sizeof(float)), 1, fp );
				fclose(fp);
			}
		}

		if( h.defmask == 'N' )
		{
			/* everything ok? - do nothing */
		}
	}

	free(data);
	free(diff);
}

/***
	does not allow for more than one reason
***/

void hist_check_qual( Histogram *h, int verbose )
{
	int i;
	float max_percent;
	float rawamp_at_max_percent;

	max_percent = h->b[0].percent;
	rawamp_at_max_percent = h->b[0].center;

	for( i = 1 ; i < h->nbin; i++ )
	{
		if( h->b[i].percent > max_percent ) 
		{
			max_percent = h->b[i].percent;
			rawamp_at_max_percent = h->b[i].center;
		}
	}

	if( max_percent > 30 )
	{
		sprintf( h->reason, 
			"max_percent = %g of amp bin = %g exceeded 30%%", 
				max_percent, rawamp_at_max_percent );
		h->defmask = 'Y';
		return;
	}

	if( h->xmax_percent > h->max_xper ) 
	{
		sprintf( h->reason,
			"xmax_percent=%g exceeded xmax_xper=%g",
				h->xmax_percent, h->max_xper );
		h->defmask = 'Y';
		return;
	}

	if( h->mode > h->maxstat )
	{
		sprintf( h->reason, "mode exceeded maxstat=%g", h->maxstat );
		h->defmask = 'Y';
		return;
	}

	if( h->xmax > h->max_x )
	{
		sprintf( h->reason, "xmax exceeded max_x=%g", h->max_x );
		h->defmask = 'Y';
		return;
	}

	if( h->mode < h->minstat )
	{
		sprintf( h->reason, "mode below minstat=%g", h->minstat );
		h->defmask = 'Y';
		return;
	}
}

void dbwrite( Histogram *h, SacFile *sf, int idiff, long evid, char *author )
{
	FILE *fp;
	char filename[128];
	char type[32];

	if( idiff == 0 ) sprintf( type, "nodiff" );
	if( idiff == 1 ) sprintf( type, "diff" );

	sprintf( filename, "sacqc_insert.sql" );
	
	if( ( fp = fopen( filename, "w" ) ) == NULL )
	{
		fprintf( stderr, "%s: ERROR cannot open file %s for output\n",
			progname, filename );
		exit(-1);
	}


	fprintf( fp, "INSERT into QCID_SEQ ( lddate ) VALUES ( CURRENT_TIMESTAMP ); \n" );

	fprintf( fp, "INSERT into MT_DATA_QUALITY\n" );
	fprintf( fp, " ( qcid, evid, net, sta, loc, chan, time, endtime, npts, delta, type, \n" );
	fprintf( fp, "   qc_mode, qc_mean, qc_median, xmax, xmax_percent, defmask, \n" );
	fprintf( fp, "   mask_reason, auth, lddate ) \n" );
	fprintf( fp, "VALUES\n" );
	fprintf( fp, " ( (SELECT max(qcid) from QCID_SEQ), \n" );
	fprintf( fp,
	  " %ld, '%s', '%s', '%s', '%s', %lf, %lf, %d, %g, '%s', %g, %g, %g, %g, %g, '%c', '%s', '%s:%s', CURRENT_TIMESTAMP );\n",
		evid,
		sf->s.knetwk,
		sf->s.kstnm,
		sf->s.khole,
		sf->s.kcmpnm,
		sf->beg.epoch,
		sf->end.epoch,
		sf->s.npts,
		sf->s.delta,
		type,
		h->mode,
		h->mean,
		h->median,
		h->xmax,
		h->xmax_percent,
		h->defmask,
		h->reason,
		progname, 
		author );
	fclose(fp);
}

void create_db_table( void )
{
	FILE *fp;
	fp = fopen( "create.sql", "w" );

	fprintf( fp, "DROP TABLE MT_DATA_QUALITY;\n" );
	fprintf( fp, "DROP TABLE QCID_SEQ;\n" );
	fprintf( fp, "\n" );

	fprintf( fp, "CREATE TABLE QCID_SEQ ( qcid INTEGER PRIMARY KEY AUTOINCREMENT, lddate DATETIME default CURRENT_TIMESTAMP);\n");
	fprintf( fp, "INSERT INTO QCID_SEQ ( qcid, lddate ) VALUES ( 1, CURRENT_TIMESTAMP );\n" );
	fprintf( fp, "\n" );

	fprintf( fp, "CREATE TABLE MT_DATA_QUALITY ( \n" );
	fprintf( fp, "	qcid		NUMBER(15)	not null,\n" );
	fprintf( fp, "	evid		NUMBER(15) 	default  -1,\n");
	fprintf( fp, "	net		VARCHAR2(8) 	default  '-', \n");
	fprintf( fp, "	sta		VARCHAR2(8) 	not null,\n" );
	fprintf( fp, "	loc		VARCHAR2(8) 	default  '-', \n");
	fprintf( fp, "	chan		VARCHAR2(8) 	not null,\n" );
	fprintf( fp, "	time		FLOAT(32)	not null,\n" );
	fprintf( fp, "	endtime		FLOAT(32)	not null,\n" );
	fprintf( fp, "	npts		NUMBER(9)	not null,\n" );
	fprintf( fp, "	delta		FLOAT(24)	not null,\n" );
	fprintf( fp, "	type		VARCHAR2(16)	not null,\n" );
	fprintf( fp, "	qc_mode		FLOAT(24)	not null,\n" );
	fprintf( fp, "	qc_mean		FLOAT(24)	not null,\n" );
	fprintf( fp, "	qc_median	FLOAT(24)	not null,\n" );
	fprintf( fp, "	xmax		FLOAT(24)	not null,\n" );
	fprintf( fp, "	xmax_percent	FLOAT(24)	not null,\n" );
	fprintf( fp, "	defmask		VARCHAR2(1)	not null,\n" );
	fprintf( fp, "	mask_reason	VARCHAR2(32)	default  '-', \n");
	fprintf( fp, "	auth		VARCHAR2(20)	default  '-', \n");
	fprintf( fp, "	lddate		DATETIME     default CURRENT_TIMESTAMP ); \n" );
	fprintf( fp, "\n" );
	fclose(fp);
}

void write_bins_percentage( Histogram *h, SacFile *sf, int idiff )
{
	int i;
	FILE *fp;
	char filename[128];
	char type[32];
	char khole[8];

	if( idiff == 0 ) sprintf( type, "nodiff" );
	if( idiff == 1 ) sprintf( type, "diff" );

	sprintf( filename, "sacqc.out" );
	if( (fp=fopen(filename, "w")) == NULL )
	{
		fprintf( stderr, "%s: ERROR cannot open file %s for output\n",
			progname, filename );
		exit(-1);
	}

/*** this is the first bin ***/

	if( strcmp( sf->s.khole, "" ) == 0 )
                strcpy( khole, "-");
        else
                strcpy( khole, sf->s.khole );

	fprintf( fp,
	  "> %s %s %s %s %s mode=%g mean=%g median=%g npts=%d dt=%g xmax=%g xmax_percentage=%g\n", 
		sf->s.knetwk, 
		sf->s.kstnm, 
		khole, 
		sf->s.kcmpnm,
		type, 
		h->mode, 
		h->mean, 
		h->median, 
		sf->s.npts, 
		sf->s.delta, 
		h->xmax, 
		h->xmax_percent );

	fprintf( fp, "%g   0\n", h->b[0].start );
	fprintf( fp, "%g %lf\n", h->b[0].start,  h->b[0].percent );
	fprintf( fp, "%g %lf\n", h->b[0].stop, h->b[0].percent );

/*** following bins ***/

	for( i = 1; i < h->nbin; i++ )
	{
		if( h->b[i].percent == 0 )
		{
			fprintf( fp, "%g 0\n", h->b[i].start );
			fprintf( fp, "%g 0\n", h->b[i].stop );
		}
		else if( h->b[i-1].percent == 0 )  /*** if previous bin was zero ***/
		{
			fprintf( fp, "%g   0\n", h->b[i].start );
			fprintf( fp, "%g %lf\n", h->b[i].start, h->b[i].percent );
			fprintf( fp, "%g %lf\n", h->b[i].stop, h->b[i].percent );
		}
		else
		{
			fprintf( fp, "%g %lf\n", h->b[i].start, h->b[i].percent );
			fprintf( fp, "%g %lf\n", h->b[i].stop, h->b[i].percent );
		}
	}
	fclose(fp);
}

void write_bins( Histogram *h, SacFile *sf, int idiff )
{
	int i;
	FILE *fp;
	char filename[128];
        char type[32];
	char khole[8];

        if( idiff == 0 ) sprintf( type, "nodiff" );
        if( idiff == 1 ) sprintf( type, "diff" );

	sprintf( filename, "sacqc.out" );
        if( (fp=fopen(filename, "w")) == NULL )
        {
                fprintf( stderr, "%s: ERROR cannot open file %s for output\n",
                        progname, filename );
                exit(-1);
        }

/*** this is the first bin ***/

	if( strcmp( sf->s.khole, "" ) == 0 ) 
		strcpy( khole, "-");
	else
		strcpy( khole, sf->s.khole );

        fprintf( fp,
          "> %s %s %s %s %s mode=%g mean=%g median=%g npts=%d dt=%g xmax=%g xmax_percentage=%g\n",
                sf->s.knetwk,
                sf->s.kstnm,
                khole,
                sf->s.kcmpnm,
                type,
                h->mode,
                h->mean,
                h->median,
                sf->s.npts,
                sf->s.delta,
                h->xmax,
                h->xmax_percent );

        fprintf( fp, "%g   0\n", h->b[0].start );
        fprintf( fp, "%g %lf\n", h->b[0].start,  h->b[0].counts );
        fprintf( fp, "%g %lf\n", h->b[0].stop, h->b[0].counts );

/*** following bins ***/

        for( i = 1; i < h->nbin; i++ )
        {
                if( h->b[i].counts == 0 )
                {
                        fprintf( fp, "%g 0\n", h->b[i].start );
                        fprintf( fp, "%g 0\n", h->b[i].stop );
                }
                else if( h->b[i-1].counts == 0 )  /*** if previous bin was zero ***/
                {
                        fprintf( fp, "%g   0\n", h->b[i].start );
                        fprintf( fp, "%g %lf\n", h->b[i].start, h->b[i].counts );
                        fprintf( fp, "%g %lf\n", h->b[i].stop, h->b[i].counts );
                }
                else
                {
                        fprintf( fp, "%g %lf\n", h->b[i].start, h->b[i].counts );
                        fprintf( fp, "%g %lf\n", h->b[i].stop, h->b[i].counts );
                }
        }
        fclose(fp);
}

/*** this script is for just single file ***/

void create_GMT5x_script2( Histogram *h, SacFile *sf, int idiff, Statistics *st, int verbose )
{
	FILE *fp;

	fp = fopen( "sacqc_gmt.csh", "w" );
	fprintf( fp, "#!/bin/csh\n" );
	fprintf( fp, "set PRJ=\"-JX5i/5i\"\n" );
        fprintf( fp, "set REG=\"-R0/10/0/30\"\n" );

	fprintf( fp, "set NET=%s\n", 	sf->s.knetwk );
	fprintf( fp, "set STA=%s\n",    sf->s.kstnm );
	fprintf( fp, "set LOC=%s\n",    sf->s.khole );
	fprintf( fp, "set CHAN=%s\n",   sf->s.kcmpnm );

	fprintf( fp, "set PS=${NET}.${STA}.${LOC}.${CHAN}.sacqc.ps\n" );

        fprintf( fp, "psxy sacqc.out ${REG} ${PRJ} -W1p,red -L -P -V0 -K >! ${PS}\n" );

	fprintf( fp, "psxy ${REG} ${PRJ} -W1p,green,5_2:0p -O -V0 -K >> ${PS} << EOF\n" );
	fprintf( fp, "> \n" );
	fprintf( fp, "%g 0\n", h->mode );
	fprintf( fp, "%g 10\n", h->mode ); 
	fprintf( fp, "EOF\n" );

	fprintf( fp, "psbasemap ${REG} ${PRJ} -U\" ${NET}.${STA}.${LOC}.${CHAN} \" " );
        fprintf( fp, "  -Bxf0.1a1+l\"Log10 Raw Amplitude (counts)\" " );
        fprintf( fp, "  -Byf5a10+l\"Percent (counts)\" -BnSeW -O -V0 -K >> ${PS}\n" );

	fprintf( fp, "pstext -R0/1/0/1 -JX5i/5i -F+jMR+f10p,Times-Bold,black -D0i/0i -O >> ${PS} << EOF\n" );
	fprintf( fp, "0.95 0.95 %s.%s.%s.%s\n", sf->s.knetwk, sf->s.kstnm, sf->s.khole, sf->s.kcmpnm );
	fprintf( fp, "0.95 0.90 stats: min=%g max=%g mean=%g median=%g\n", 
			st->min, st->max, st->ave, st->med );
	fprintf( fp, "0.95 0.85 stats: skew=%g curt=%g\n", st->skew, st->curt );
	fprintf( fp, "0.95 0.80 hist: mode=%g mean=%g median=%g xmax=%g xmax_per=%g\n",
		h->mode,
                h->mean,
                h->median,
                h->xmax,
                h->xmax_percent );
	fprintf( fp, "EOF\n" );

	fprintf( fp, "psconvert -Tj -E300 -A ${PS}\n" );
	fprintf( fp, "/bin/rm -f ${PS}\n" );
	fclose(fp);

/*** run the script ***/

	chmod( "sacqc_gmt.csh", S_IRWXU|S_IRWXG|S_IRWXO );
	system( "/bin/csh ./sacqc_gmt.csh" );
}

/*** this script plots all files in directory associated with at least one event ***/

void create_GMT5x_script( )
{
	FILE *fp;

	fp = fopen( "sacqc.csh", "w" );
	fprintf( fp, "#!/bin/csh\n" );

	fprintf( fp, "set PS=nodiff.ps\n" );
	fprintf( fp, "set JPG=nodiff.jpg\n" );
	fprintf( fp, "\n" );

	fprintf( fp, "#cleanup\n" );
	fprintf( fp, "#/bin/rm -f sacqc.out sacqc_insert.sql create.sql %{PS} %{JPG}\n" );
	fprintf( fp, "\n" );

	fprintf( fp, "set SACFILES=(`/bin/ls -1 *.SAC`)\n" );
	fprintf( fp, "set maxfiles=$#SACFILES\n" );
	fprintf( fp, "\n" );

	fprintf( fp, "cat >! sacqc.par << EOF\n" );
	fprintf( fp, "nodiff\n" );
	fprintf( fp, "#diff\n" );
	fprintf( fp, "xlog\n" );
	fprintf( fp, "ylog\n" );
	fprintf( fp, "xstart=0\n" );
	fprintf( fp, "xstop=10\n" );
	fprintf( fp, "binwidth=0.1\n" );
	fprintf( fp, "percent\n" );
	fprintf( fp, "noverbose\n" );
	fprintf( fp, "db\n" );
	fprintf( fp, "#evid=9999999999\n" );
	fprintf( fp, "#auth=\"ichinose1\"\n" );
	fprintf( fp, "writeback\n" );
	fprintf( fp, "EOF\n" );
	fprintf( fp, "\n" );

	fprintf( fp, "set PRJ=\"-JX5i/5i\"\n" );
	fprintf( fp, "set REG=\"-R0/10/0/30\"\n" );

	/* fprintf( fp, " -U\"%s %s\" ", author, evid );  */

	fprintf( fp, "\n" );
	fprintf( fp, "### loop over SACFILE\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "set i=1\n" );
	fprintf( fp, "foreach SACFILE ( *.SAC )\n" );

	fprintf( fp, "set PS=${SACFILE}.ps\n" );

	fprintf( fp, "psbasemap ${REG} ${PRJ} -U\" ${SACFILE} \" " );
	fprintf( fp, "  -Bxf0.1a1+l\"Log10 Raw Amplitude (counts)\" -Byf5a10+l\"Percent Counts\" -BnSeW -P -V0 -K >! ${PS}\n" );

	fprintf( fp, "sacqc par=sacqc.par ${SACFILE}\n" );

	fprintf( fp, "sqlite3 ../sacqc.db << EOF\n" );
	fprintf( fp, ".read sacqc_insert.sql\n" );
	fprintf( fp, ".quit\n" );
	fprintf( fp, "EOF\n" );
	fprintf( fp, "\n" );

	fprintf( fp, "set line_color=( `switch_color 3 ${i}` )\n" );
	fprintf( fp, "if( ${i} < ${maxfiles} ) then\n" );
	fprintf( fp, " psxy sacqc.out ${REG} ${PRJ} -W1p,${line_color} -L -O >> ${PS}\n" );
	fprintf( fp, "else\n" );
	fprintf( fp, " psxy sacqc.out ${REG} ${PRJ} -W1p,${line_color} -L -O >> ${PS}\n" );
	fprintf( fp, "endif\n" );
	fprintf( fp, "\n" );

	fprintf( fp, "psconvert -Tj -E300 -A ${PS}\n" );

	fprintf( fp, "@ i += 1\n" );
	fprintf( fp, "end ### loop over SACFILE\n" );
	fprintf( fp, "\n" );

	fprintf( fp, "# psconvert -Tj -E300 -A ${PS}\n" );
	fprintf( fp, "\n" );

	fprintf( fp, "# /bin/rm -f sacqc.par sacqc.out sacqc_insert.sql\n" );
	fprintf( fp, "\n" );

	fclose(fp);
}

void Usage( void )
{
	fprintf( stderr, "\n" );
	fprintf( stderr, "%s : Usage: \n", progname );
	fprintf( stderr, "\t %s xstart=(float) xstop=(float) binwidth=(float) [no]diff [no]xlog [no]ylog [no]percent \n", progname );
	fprintf( stderr, "\t\t maxstat=(float) max_x=(float) max_xper=(float) evid=(long int) [no]db auth=(string) \n" );
	fprintf( stderr, "\t\t  [no]writeback [no]verbose [no]gmt5 [no]help  {list of SAC files} *.SAC \n" );
	fprintf( stderr, "\n" );

	fprintf( stderr, "Required: \n" );
	fprintf( stderr, "\t  List of SAC files. Wildcards acceptable e.g., *.SAC\n" );
	fprintf( stderr, "\n" );

	fprintf( stderr, "Optional: \n" );
	fprintf( stderr, "\t xstart=(float)   starting bin of Log10 amplitude histogram [default 0]\n" );
	fprintf( stderr, "\t xstop=(float)    maximum bin of Log10 amplitude histogram [default 10]\n" );
	fprintf( stderr, "\t binwidth=(float) bin width [default 0.1]\n" );
	fprintf( stderr, "\t [no]diff         do the difference between amplitudes rather than abs amp [default off]\n" );
	fprintf( stderr, "\t [no]xlog         amplitude log [default on]\n" );
	fprintf( stderr, "\t [no]ylog         count log [default on]\n" );
	fprintf( stderr, "\t [no]percent      plot y-axis in percent rather than counts [default on]\n" );
	fprintf( stderr, "\t maxstat=(float)  the log10 of the mode, mean or median cannot be above this value [default 5]\n" );
	fprintf( stderr, "\t minstat=(float)  the log10 of the mode, mean or median cannot be below this value [default 0.5]\n" );
	fprintf( stderr, "\t max_x=(float)    the max log10 amp cannot be above this value [default 9] \n" );
	fprintf( stderr, "\t max_xper=(float) the max amp bin cannot contain more than this %% of the total amp counts [default 2%%]\n" );
	fprintf( stderr, "\t evid=(long int)  unique event identification for database write [default -1]\n" );
	fprintf( stderr, "\t [no]db           write for sqlite3 database insert statements [default on]\n" );
	fprintf( stderr, "\t auth=(string)    author name for database write [default null]\n" );
        fprintf( stderr, "\t [no]writeback    if QC fails then write to SAC header kuser2=SACQC user2=999\n" );
        fprintf( stderr, "\t [no]verbose      verbosity or silent [default off]\n" );
        fprintf( stderr, "\t [no]gmt5         create plot using GMT version 5.x.x compatible commands [default on]\n" );
	fprintf( stderr, "\t [no]help         print this help usage [default off]\n" );
        fprintf( stderr, "\t {list of SAC files} *.SAC \n" );
	fprintf( stderr, "\n" );
	fprintf( stderr, "\t This program creates a histogram of the input SAC data. Percentage or Counts vs. Raw amplitudes\n" );
	fprintf( stderr, "\t The distribution of raw seismic data has a unique distribution.  Data with QC problems exhibit\n" );
	fprintf( stderr, "\t exotic distributions that are detectable when statistical values go out of range.\n" );
	fprintf( stderr, "\n" );
	fprintf( stderr, "\t output: histogram in GMT JPEG plot saved as {net}.{sta}.{loc}.{chan}.sacqc.jpg\n" );
	fprintf( stderr, "\t         sacqc_insert.sql create.sql  - sqlite3 scripts for loading\n" );
	fprintf( stderr, "\n" );
}
