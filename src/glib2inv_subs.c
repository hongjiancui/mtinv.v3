#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/mt.h"

char progname[128];

EventInfo *glib2inv_get_input_parameters( char *filename, 
	EventInfo *ev, int *n, int verbose )
{
	FILE *fp;
	int MAX_RECORD_LENGTH = 512;
	char rec[MAX_RECORD_LENGTH];
	int ista, nitems;
	MyTime ot;
	char timestring[24], kdum[3], comment[256], database_info[256];
	char kenvelope, grd_mo_type, kused;
	float strike, dip, rake, Mw, evdp, evla, evlo;
	float Nyquist_Frequency;

/*** timesubs.o ***/
	void parsestring( MyTime *, char * );
	void clone_mytime(  MyTime *, MyTime * );

/******************************/
/*** open glib2inv.par file ***/
/******************************/
	if( (fp = fopen( filename, "r" )) == NULL )
	{
		printf("%s: cannot open file %s\n", progname, filename );
		exit(-1);
	}

/********************/
/*** defaults     ***/
/********************/

	grd_mo_type='d';
	strike=-999; 
	dip=-999; 
	rake=-999;
	Mw=-999;
	ista = 0;
	while( fgets( rec, MAX_RECORD_LENGTH, fp ) != NULL )
	{

	/********************/
	/*** comment line ***/
	/********************/
		if( rec[0] == '#' ) continue;

	/***********************************************************/
	/*** AutoMT includes this database info line. The home   ***/
	/*** version skips it but we have it here to maintain    ***/
	/*** compatibility between operational and home versions ***/
	/***********************************************************/
		if( strncmp( rec, "DB ", 3 ) == 0 )
		{
			/*** just skip for now, G. Ichinose 2014 ***/
			/*
			sscanf( rec, "%s %[^\n]\n", kdum, database_info );
			if(verbose)printf("%s: glib2inv_subs(): database info=%s\n", 
				progname, database_info );
			*/
                        continue;
		}

	/*******************************************************/
	/*** event discription line (gets written into email)***/
	/*******************************************************/
		if( strncmp( rec, "CM ", 3 ) == 0 )
		{
			sscanf( rec, "%s %[^\n]\n", kdum, comment );
			if(verbose)printf("%s: glib2inv_subs(): comment=%s\n", progname, comment );
			continue;
		}

	/*******************************************************/
	/*** read the origin time                            ***/
	/*******************************************************/

		if( strncmp( rec, "OT ", 3 ) == 0 ) 
		{
			sscanf( rec, "%s %s", kdum, timestring );
			parsestring( &ot, timestring );
			if( verbose ) 
			{
			  printf("%s: glib2inv_subs(): origin time=\t", progname );
			  WriteMyTime2STDOUT( &ot );
			}
			continue;
		}

	/************************************************************************/
	/*** read the event information, only evlo, evla needed for inversion ***/
	/************************************************************************/

		if( strncmp(rec, "EV ", 3 ) == 0 )
		{
			sscanf( rec, "%s %f %f %f %f %f %f %f", 
				kdum, &strike, &dip, &rake, &Mw, &evlo, &evla, &evdp );

			if( verbose )
			{
			  printf("%s: glib2inv_subs(): str=%g dip=%g rak=%g Mw=%g evla=%g evlo=%g evdp=%g\n",
				progname, strike, dip, rake, Mw, evla, evlo, evdp );
			}
			continue;
		}
	
	/************************************************************************/
	/*** read the station information                                     ***/
	/************************************************************************/

		ev = (EventInfo *)realloc( ev, (ista+1)*sizeof(EventInfo) );

		/****                  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 ***/
		nitems = sscanf( rec, "%s %s %s %d %d %f %f %d %f %f %f %c %f %c %f %f", 
			ev[ista].stnm, 			/*  1 */
			ev[ista].net,			/*  2 */
			ev[ista].modfile,		/*  3 */
			&(ev[ista].npole),		/*  4 */
			&(ev[ista].npass),		/*  5 */
			&(ev[ista].lf),			/*  6 */
			&(ev[ista].hf),			/*  7 */
			&(ev[ista].nt),			/*  8 */
			&(ev[ista].dt),			/*  9 */
			&(ev[ista].tr),                 /* 10 */
			&(ev[ista].tt),			/* 11 */
			&grd_mo_type,			/* 12 */
			&(ev[ista].mul_factor),		/* 13 */
			&kused,				/* 14 */
			&(ev[ista].time_shift_all),     /* 15 */
			&(ev[ista].weight) 		/* 16 */
		);

	/***************************************/
	/*** check time number of items read ***/
	/***************************************/
		if( nitems < 16 )
		{
		  fprintf(stderr, "%s: not enough items=%d read from file %s ista=%d\n\n",
			progname, nitems, filename, ista );
		  fprintf(stderr, "%s: offending line=%s\n", progname, rec );
		  fprintf(stderr, "\n" );
		  fprintf(stderr, "%s: free-format, space seperated with 18 columns: \n", progname );
	          fprintf(stderr, "1   2   3   4     5     6    7    8   9    10 11 12 13   14   15  16     \n" );
		  fprintf(stderr, "sta net mod npole npass lf   hf   nt  dt   tr tt GM fmul used TS  weight \n" );
		  fprintf(stderr, "--- --- --- ----- ----- ---- ---- --- ---- -- -- -- ---- ---- --  ------ \n");
		  fprintf(stderr, "PAS CI  wus 3     2     0.02 0.05 512 0.15 0  0  d  1    y    0   1.0    \n" );
		  fprintf(stderr, "\n" );
		  fprintf(stderr, "Column  1: Station Code\n" );
		  fprintf(stderr, "Column  2: Netwok Code\n" );
	          fprintf(stderr, "Column  3: 1D Velocity Model Base File Name (without .mod extension)\n" );
    		  fprintf(stderr, "Column  4: Butterworth Filter number of poles\n" );
		  fprintf(stderr, "Column  5: Butterworth Filter number of passes\n" );
		  fprintf(stderr, "Column  6: Butterworth Filter lowpass corner (Hz)\n" );
		  fprintf(stderr, "Column  7: Butterworth Filter highpass corner (Hz)\n" );
		  fprintf(stderr, "Column  8: Decimated number of points\n" );
		  fprintf(stderr, "Column  9: Decimated Sampling Rate (sec/samp)\n" );
		  fprintf(stderr, "Column 10: Rise Time of Boxcar (sec)\n" );
		  fprintf(stderr, "Column 11: Duration of Boxcar (sec)\n" );
		  fprintf(stderr, "Column 12: Ground Motion (d=displacement v=velocity)\n" );
		  fprintf(stderr, "Column 13: Multuply all amplitudes of data by this value\n" );
		  fprintf(stderr, "Column 14: y=use in calculation and residual n=do not use only predict\n" );
		  fprintf(stderr, "Column 15: time shift all data by this value in (sec)\n" );
		  fprintf(stderr, "Column 16: Amplitude weight given to this station data for inversion\n" );
		  fprintf(stderr, "\n" );
		  exit(-1);
		}

		if( verbose ) 
		{
			fprintf(stdout, "%s: glib2inv_subs(): nitems read = %d for ista=%d from file %s\n",
				progname, nitems, ista, filename ); 
		}

	/***************************************/
	/*** check the npole range           ***/
	/***************************************/
		if( ev[ista].npole <= 0 || ev[ista].npole > 9 )
		{
			fprintf( stderr, "%s: glib2inv_subs(): npole=%d out of range [1,9]\n",
				progname, ev[ista].npole );
			exit(-1);
		}

	/***************************************/
	/*** check npass                     ***/
	/***************************************/
		if( ev[ista].npass < 1 || ev[ista].npass > 2 )
		{
			fprintf( stderr, "%s: glib2inv_subs(): npass=%d out of range [1,2]\n",
				progname, ev[ista].npass );
			exit(-1);
		}

	/***************************************/
	/*** check hf range below nyquist freq */
	/***************************************/
		Nyquist_Frequency = 1 / ( 2 * ev[ista].dt );
		if( ev[ista].hf > Nyquist_Frequency )
		{
			fprintf( stderr, 
			  "%s: glib2inv_subs(): hi-freq corner %g greater than Nyquist Frequency %g(Hz)\n",
				progname, ev[ista].hf, Nyquist_Frequency );

			fprintf( stderr,
			  "%s: glib2inv_subs(): lower hf for station %s.%s in file %s\n", 
				progname, ev[ista].stnm, ev[ista].net, filename );

			exit(-1);
		}

	/***************************************/
	/*** check if hf >= lf               ***/
	/***************************************/
		if( ev[ista].lf >= ev[ista].hf )
		{
			fprintf( stderr,
			  "%s: glib2inv_subs(): low-frequency corner %g greater than high-frequency corner %g\n",
				progname, ev[ista].lf, ev[ista].hf );

			fprintf( stderr,
                          "%s: glib2inv_subs(): increase hf or lower lf for station %s.%s in file %s\n",
                                progname, ev[ista].stnm, ev[ista].net, filename );
                        exit(-1);
		}

	/***************************************/
	/*** check mul factors               ***/
	/***************************************/
		if( ev[ista].mul_factor == 0 ) ev[ista].mul_factor = 1;

	/***********************************************************/
	/*** check if station is being used in the inversion of  ***/
	/*** is being carried along for prediction               ***/
	/***********************************************************/

		ev[ista].iused = 1; /** default catch all, station is being used **/
		if( kused == 'y' || kused == 'Y' || kused == '1' ) ev[ista].iused = 1;
		if( kused == 'n' || kused == 'N' || kused == '0' ) ev[ista].iused = 0;

	/*********************************************/
	/*** check to process as envelope          ***/
	/***   default is off                      ***/
	/*********************************************/
		kenvelope = 'n';
		if( kenvelope == 'y' || kenvelope == 'Y' || kenvelope == '1' ) ev[ista].ienvelope = 1;
		if( kenvelope == 'n' || kenvelope == 'N' || kenvelope == '0' ) ev[ista].ienvelope = 0;
		ev[ista].ienvelope = 0;

	/***************************************************/
	/*** ground motion type velocity or displacement ***/
	/***************************************************/

		if( grd_mo_type == 'v' )
		{
			 ev[ista].grd_mo_type = VELOCITY;
		}
		else if( grd_mo_type == 'd' )
		{
			ev[ista].grd_mo_type = DISPLACEMENT;
		}
		else
		{
		  ev[ista].grd_mo_type = DISPLACEMENT;
		  fprintf(stderr, 
		    "%s: glib2inv_subs(): Warning unknown Ground motion type=%c, assuming displacement\n",
			progname, grd_mo_type );
		}

	/******************************************/
	/*** create file names for output files ***/
	/******************************************/

		sprintf( ev[ista].data_filename, "%s.%s.%c.%02d.data", 
			ev[ista].stnm, ev[ista].net, grd_mo_type, ista );

		sprintf( ev[ista].glib_filename, "%s.%s.%s.glib",
			ev[ista].stnm, ev[ista].net, ev[ista].modfile );

		sprintf( ev[ista].ginv_filename, "%s.%s.%s.%c.%02d.ginv", 
			ev[ista].stnm, ev[ista].net, ev[ista].modfile, grd_mo_type, ista );

		if( verbose )
		{
		  printf( "%s: glib2inv_subs():           data file=%s\n", progname, ev[ista].data_filename );
		  printf( "%s: glib2inv_subs(): greens ftn lib file=%s\n", progname, ev[ista].glib_filename );
		  printf( "%s: glib2inv_subs(): greens ftn inv file=%s\n", progname, ev[ista].ginv_filename );
		}
	
	/***********************************************************************/
	/*** set some of the global data structures with the local variables ***/
	/***********************************************************************/

		clone_mytime( &ot, &(ev[ista].ot) );
		clone_mytime( &ot, &(ev[ista].ot_orig) );

		ev[ista].str  = strike;
		ev[ista].dip  = dip;
		ev[ista].rak  = rake;
		ev[ista].Mw   = Mw;
		ev[ista].my_z = evdp;
		ev[ista].evlo = evlo;
		ev[ista].evla = evla;
		ev[ista].evdp = evdp;
		strcpy( ev[ista].comment, comment );

		if( verbose )
		{
		  printf( "%s: glib2inv_subs(): ista=%03d data=%s glib=%s ginv=%s npole=%d npass=%d\n",
			progname,
			ista,
			ev[ista].data_filename,
			ev[ista].glib_filename,
			ev[ista].ginv_filename,
			ev[ista].npole, 
			ev[ista].npass );
		  printf( "%s: glib2inv_subs(): lf=%g hf=%g nt=%d dt=%g tr=%g tt=%g grd_mo_type=%d mul_factor=%g iused=%d\n",
			progname, 
			ev[ista].lf, ev[ista].hf, ev[ista].nt, 
			ev[ista].dt, ev[ista].tr, ev[ista].tt,
			ev[ista].grd_mo_type, ev[ista].mul_factor, ev[ista].iused );
		}
		ista++;
	}
	*n = ista;
	return (EventInfo *)ev;
}

void array2grn( float **garray, Greens *g )
{
        int it, nt;
        nt = g->nt;
	if( nt > 4096 ) printf("%s: array2grn: nt=%d > 4096\n", progname, nt );

	for( it=0; it<nt; it++ )
	{
		g->g.rss[it] = 0;
                g->g.rds[it] = 0;
                g->g.rdd[it] = 0;
                g->g.rep[it] = 0;
                g->g.zss[it] = 0;
                g->g.zds[it] = 0;
                g->g.zdd[it] = 0;
                g->g.zep[it] = 0;
                g->g.tss[it] = 0;
                g->g.tds[it] = 0;
	}

        for( it=0; it<nt; it++ )
        {
                g->g.rss[it] = garray[0][it];
                g->g.rds[it] = garray[1][it];
                g->g.rdd[it] = garray[2][it];
                g->g.rep[it] = garray[3][it];
                g->g.zss[it] = garray[4][it];
                g->g.zds[it] = garray[5][it];
                g->g.zdd[it] = garray[6][it];
                g->g.zep[it] = garray[7][it];
                g->g.tss[it] = garray[8][it];
                g->g.tds[it] = garray[9][it];
        }
}

void split2grn( Greens *g, float **garray )
{
        int it, nt;
        nt = g->nt;
	
        for( it=0; it<nt; it++ )
        {
                garray[0][it] = g->g.rss[it];
                garray[1][it] = g->g.rds[it];
                garray[2][it] = g->g.rdd[it];
                garray[3][it] = g->g.rep[it];
                garray[4][it] = g->g.zss[it];
                garray[5][it] = g->g.zds[it];
                garray[6][it] = g->g.zdd[it];
                garray[7][it] = g->g.zep[it];
                garray[8][it] = g->g.tss[it];
                garray[9][it] = g->g.tds[it];
        }
}

void grn2disp( Greens *g, EventInfo *ev, int verbose, int mtdegfree )
{

/*** ten greens functions ***/
        float *rss, *rds, *rdd, *rep, *zss, *zds, *zdd, *zep, *tss, *tds;

/*** directional cosine coefficients ***/
        float a1, a2, a3, a4, a5;
        Tensor M; /* float Mxx, Myy, Mzz, Mxy, Mxz, Myz; */
        float strr, dipr, rakr;
	float half=0.5;
	float Mo;

        int it;
        float pi, d2r, dt, t0, e, tt, tr, fi;
        int nt;
	void source_time_function( float *, int, float, float, float );

/*** set some constants ***/
        pi  = M_PI;
        d2r = pi/180.;

        nt = g->nt;
        dt = g->dt;
        t0 = g->t0;
        e  = t0 + (nt*dt);
        tt = ev->tt;
        tr = ev->tr;

/**** set t0 from reduction velocity redv ****/
        t0 = g->tstart;
	e  = g->tend;

/*** assign greens function to local pointers ***/

        rss = g->g.rss;
        rds = g->g.rds;
        rdd = g->g.rdd;
        rep = g->g.rep;
        zss = g->g.zss;
        zds = g->g.zds;
        zdd = g->g.zdd;
        zep = g->g.zep;
        tss = g->g.tss;
        tds = g->g.tds;

/** convolve a source time function ***/

        source_time_function( rss, nt, dt, tr, tt );
        source_time_function( rds, nt, dt, tr, tt );
        source_time_function( rdd, nt, dt, tr, tt );
        source_time_function( rep, nt, dt, tr, tt );
        source_time_function( zss, nt, dt, tr, tt );
        source_time_function( zds, nt, dt, tr, tt );
        source_time_function( zdd, nt, dt, tr, tt );
        source_time_function( zep, nt, dt, tr, tt );
        source_time_function( tss, nt, dt, tr, tt );
        source_time_function( tds, nt, dt, tr, tt );

/* w = transverse v = vertical u = radial */

        fi = g->az * d2r;
        ev->Mo = pow( 10.0, (1.5*( ev->Mw + 10.73)) );
	Mo = ev->Mo / base_moment;

        strr = ev->str * d2r;
        dipr = ev->dip * d2r;
        rakr = ev->rak * d2r;

        M.xx = -( sin(dipr) * cos(rakr) * sin(2*strr) + sin(2*dipr) * sin(rakr) * sin(strr)*sin(strr) );
        M.yy =  ( sin(dipr) * cos(rakr) * sin(2*strr) - sin(2*dipr) * sin(rakr) * cos(strr)*cos(strr) );
        M.zz =  ( sin(2*dipr) * sin( rakr ) );
        M.xy =  ( sin(dipr) * cos(rakr) * cos(2*strr) + 0.5*sin(2*dipr) * sin(rakr) * sin(2*strr) );
        M.xz = -( cos(dipr) * cos(rakr) * cos(strr) + cos(2*dipr) * sin(rakr) * sin(strr) );
        M.yz = -( cos(dipr) * cos(rakr) * sin(strr) - cos(2*dipr) * sin(rakr) * cos(strr) );

	if( verbose )
	{
          printf("%s: str=%g dip=%g rak=%g Mw=%g Mo=%g\n",
               	progname, ev->str, ev->dip, ev->rak, ev->Mw, ev->Mo );

          printf("%s: fi=%g strr=%g dipr=%g rakr=%g Mxx=%g Myy=%g Mzz=%g Mxy=%g Mxz=%g Myz=%g\n",
               	progname, fi, strr, dipr, rakr, M.xx, M.yy, M.zz, M.xy, M.xz, M.yz);
	}

/*** compute the coefficients ****/
	a1 = half * ( M.xx - M.yy ) * cos( 2 * fi ) + M.xy * sin( 2 * fi );
        a2 = M.xz * cos( fi ) + M.yz * sin( fi );
        a3 = -half*( M.xx + M.yy );
        a4 = half * ( M.xx - M.yy ) * sin( 2 * fi ) - M.xy * cos( 2 * fi );
        a5 = -M.yz * cos( fi ) + M.xz * sin( fi );

        if( verbose ) 
	{
	   printf("%s: a1=%g a2=%g a3=%g a4=%g a5=%g mtdegfree=%d\n", 
		progname, a1, a2, a3, a4, a5, mtdegfree );
	}
	
        for( it=0; it<nt; it++)
        {
		if( mtdegfree == 5 || mtdegfree == 6 )
		{
                  g->ver[it] = (a1*zss[it]+a2*zds[it]+a3*zdd[it])*Mo; /*** vertical   ***/
                  g->rad[it] = (a1*rss[it]+a2*rds[it]+a3*rdd[it])*Mo; /*** radial     ***/
                  g->tra[it] = (a4*tss[it]+a5*tds[it])*Mo;            /*** transverse ***/
		}
		else if( mtdegfree == 1 )
		{
		  g->ver[it] = zep[it] * Mo;
		  g->rad[it] = rep[it] * Mo;
		  g->tra[it] = 0;
		}
        }

	if( verbose )
	  printf("%s: done with forward calculation\n", progname );

        return;
}

void special_load( char *station_name, char *network_name, Greens *g )
{
	float **garray, *data;
	Sac_Header *s;	
	FILE *fp;
	int i, iext, it, num_ext = 10;
	char filename[128];                
		   	       /***            0    1      2     3      4     5       6      7       8     9     ***/
	char filename_extension[10][4] = { "rss", "rds", "rdd", "rex", "zss", "zds", "zdd", "zex", "tss", "tds" };
	int ig, ng=10;
	void array2grn( float **, Greens * );

/**************************/
/*** zero out the array ***/
/**************************/
	garray = (float **)malloc(ng*sizeof(float *));
	for( ig=0; ig<ng; ig++ )
		garray[ig] = (float *)malloc( 4096 * sizeof(float) );

	for( iext=0; iext<num_ext; iext++ )
	{	
		for( it=0; it < g->nt; it++ )
			garray[iext][it] = 0;
	}
	
/********************************************************/
/*** loop over all green's functions for this station ***/
/********************************************************/
	for( iext=0; iext<num_ext; iext++ )
	{

	/*************************************************************/
	/*** open the sac file for this Green's function component ***/
	/*************************************************************/
		sprintf( filename, "%s.%s.%s", 
		  station_name, network_name, filename_extension[iext] );
		fprintf( stdout, "%s: file=%s(%02d) ", progname, filename, iext );
		if( (fp=fopen(filename, "rb")) == NULL )
		{
			fprintf( stderr, "\n%s: Fatal Error, cannot open file %s\n",
				progname, filename );
			exit(-1);
		}

	/********************************************************/
	/*** load the sac header and data from file           ***/
	/*** copy the sac data into Greens structure          ***/
	/********************************************************/
		s = (Sac_Header *)malloc( sizeof(Sac_Header) );
		fread( s, sizeof(Sac_Header), 1, fp );
		data = (float *)calloc( s->npts, sizeof(float) );
		fread( &data[0], s->npts * sizeof(float), 1, fp );
		fclose(fp);

		g->nt  = s->npts;

	/********************************************************/
	/*** control overflow ***/
	/********************************************************/
		if( s->npts > 4096 ) g->nt = 4096;
		for( it = 0; it < g->nt; it++ ) 
			garray[iext][it] = data[it];

	/********************************************************/
	/*** copy the sac header into Greens structure ***/
	/********************************************************/
		strncpy( g->stnm, s->kstnm, 8 );
		strncpy( g->net, s->knetwk, 8 );

		for( i=0; i<8; i++ )
		{
			if( g->stnm[i] == ' ' ) g->stnm[i]='\0';
			if( g->net[i]  == ' ' ) g->net[i]='\0';
		}
	
		g->dt		= s->delta;
		g->t0		= s->b;
		g->evdp 	= s->evdp;
		g->stla		= s->stla;
		g->stlo		= s->stlo;
		g->evla		= s->evla;
		g->evlo		= s->evlo;
		g->rdist	= s->dist;
		g->az		= s->az;
		g->baz		= s->baz;
		g->rigidity	= 3.3E+11;
		g->redv		= -1.0;
		g->ts0		= 0.0;
		g->tstart	= s->b;
		g->tend		= s->b + s->delta * s->npts;

		fprintf( stdout, " stla=%g stlo=%g evla=%g evlo=%g evdp=%g rdist=%g ",
			g->stla, g->stlo, g->evla, g->evlo, g->evdp, g->rdist );
		fprintf( stdout, " az=%g baz=%g dt=%g nt=%d b=%g sta=%s net=%s\n",
			g->az, g->baz, g->dt, g->nt, g->t0, g->stnm, g->net );

		free(s);
		free(data);
	}
	array2grn( garray, g );
}
