#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/mt.h" /** global datatype and structure declarations **/

char progname[128];

int main( int ac, char **av )
{
	extern Greens grn_;
	int iz;
	float *z;
	Depth_Info depth_info;
	FILE *fp;
	int verbose=0, idump=0, ista;
	void Usage_Print(void);
	int it;  /*** this version computes the rigidity here ***/
	float cm2km, UnitArea, rigidity, slip4pi;
	cm2km = 1.0E+05;
	UnitArea = 1 * cm2km * cm2km;

/* see include/mt.h Mo = math.pow( 10.0, 1.5*(Mw+10.73) ) = 1.2445146117713818e+16; Reference Mw = 0.0 */

	int isrclayer;
	VelMod vmint;
	float dlaykm  = 1.0;

	void greensf_( int * );
	void wrtgrn2sac( Greens *, int );
	void getparameters( int, char **, Greens *, Depth_Info *, int *, int * );

/*** from modsubs.c ***/
        void create_mod( VelMod * );
        void print_mod0( VelMod * );
        void print_mod1( FILE *, VelMod * );
        void compute_rigidity( Greens *, int );

/*** from rayp_subs.c ***/
        void interpolate_model( VelMod *, VelMod *, float );
        void calc_1D_tt( float, float, float *, float *, float *, int *, VelMod *, int );
        void earth_flattening_transformation( VelMod * );
        float calc_takeoff_angle( float, float, float, float );

/*** misc ***/
	char pathname[128];
	char *shorten_path( char *pathname, char *filename );

/*** start program main ***/

	strcpy( pathname, av[0] );
	shorten_path( pathname, progname );

	if( ac <= 1 ) 
	{
		Usage_Print();
	}

	if(verbose)
	{
		fprintf(stdout, "%s: STDOUT: Version=%s Release Date=%s exec full path=%s\n",
		  progname, Version_Label, Version_Date, pathname );
	}
	/*
	fprintf(stderr, "%s: STDERR Version=%s Release Date=%s exec full path=%s\n",
                 progname, Version_Label, Version_Date, pathname );
	*/

/*******************************************************/
/*** get input parameters from file and command line ***/
/*******************************************************/
	getparameters( ac, av, &grn_, &depth_info, &verbose, &idump );

/**********************************************************/
/*** make an interpolated version of the original model ***/
/**********************************************************/
	interpolate_model( &(grn_.v), &vmint, dlaykm );
	if(verbose) print_mod0( &vmint );
                                                                                                                                 
/**********************************************************/
/*** Earth Flattening Transformation                    ***/
/**********************************************************/
	earth_flattening_transformation( &vmint );
	if(verbose)print_mod0( &vmint );

/*******************************************************/
/*** set the depth vector                            ***/
/*******************************************************/
	z = (float *) calloc( depth_info.nz, sizeof(float) );
	if(verbose) printf("%s: \t\t", progname );
	for( iz=0; iz < depth_info.nz; iz++ )
	{
		z[iz] = depth_info.zmin + ( iz * depth_info.zinc );
		if(verbose)
		{
			if( iz == depth_info.nz - 1 )
				printf("%g\n", z[iz] );
			else
				printf("%g ", z[iz] );
		}
	}

/*******************************************************/
/*** open the file and write the source depth        ***/
/*** references in the header of the output file     ***/
/*******************************************************/
	if( (fp = fopen(grn_.filename,"wb")) == NULL )
	{
		printf("%s: Fatal Error, cannot open file %s\n", 
			progname, grn_.filename );
		exit(-1);
	}
	fwrite(&(depth_info.nz),sizeof(int),1,fp);
	fwrite(&z[0],depth_info.nz*sizeof(float),1,fp);
	fflush(fp);

/*******************************************************/
/*** loop over depth                                 ***/
/*******************************************************/
	for( iz=0; iz< depth_info.nz; iz++ )
	{

	/*******************************************************/
	/*** compute the rigidity at the source for the      ***/
	/*** source depth's velocity and density             ***/
	/*******************************************************/
		if( verbose ) printf("%s: iz=%d z=%g\n", progname, iz, z[iz] );
		grn_.evdp = z[iz];
		compute_rigidity( &grn_, verbose );

	/********************************************************/
	/*** compute the Green's functions Ref:               ***/
	/*** Zeng and Anderson (1995) Bull. Seismol. Soc. Am. ***/
	/********************************************************/
		greensf_( &verbose );

	/***********************************************************/
	/*** dump part of the structure for debugging            ***/
	/***********************************************************/
		if(verbose)
		{
		  fprintf(stdout, "%s: filename=%s stnm=%s net=%s stla=%f stlo=%f ",
		    progname, grn_.filename, grn_.stnm, grn_.net, grn_.stla, grn_.stlo );
		  fprintf(stdout, "evla=%f evlo=%f evdp=%f rdist=%g az=%f baz=%f ",
		    grn_.evla, grn_.evlo, grn_.evdp, grn_.rdist, grn_.az, grn_.baz );
		  fprintf(stdout, "t0=%f dt=%f fmax=%f damp=%f eps=%f smin=%f\n",
		    grn_.t0, grn_.dt, grn_.fmax, grn_.damp, grn_.eps, grn_.smin );
		}
		fprintf(stdout, 
	"%s: %s evdp=%6.2f rdist=%8.2f az=%3.0f nt=%4d dt=%5.3f fmax=%5.2f t0=%g twin=%g ",
			progname,
			grn_.filename,
			grn_.evdp,
			grn_.rdist,
			grn_.az,
			grn_.nt,
			grn_.dt,
			grn_.fmax,
			grn_.t0,
			grn_.twin );

	/**************************/
	/*** apply the rigidity ***/
	/**************************/
		rigidity = grn_.rigidity;
		slip4pi =  (base_moment/(rigidity * UnitArea))/(4*M_PI);
		for( it=0; it<grn_.nt; it++ )
		{
			grn_.g.rss[it] *= slip4pi;
			grn_.g.rds[it] *= slip4pi;
			grn_.g.rdd[it] *= slip4pi;
			grn_.g.rep[it] *= slip4pi;
			grn_.g.zss[it] *= slip4pi;
			grn_.g.zds[it] *= slip4pi;
			grn_.g.zdd[it] *= slip4pi;
			grn_.g.zep[it] *= slip4pi;
			grn_.g.tss[it] *= slip4pi;
			grn_.g.tds[it] *= slip4pi;
		}

	/************************************/
	/*** introduce reduction velocity ***/
	/************************************/
		grn_.ts0	= grn_.t0;
		grn_.tstart	= grn_.ts0;
		grn_.tend	= grn_.ts0 + grn_.dt * grn_.nt;

	/*********************************************************************/
	/*** calculate travel time, ray parameter, and ray bottoming depth ***/
	/*********************************************************************/
		verbose = 0;
		calc_1D_tt( grn_.rdist, 
			grn_.evdp, 
			&(grn_.Prayparameter), 
			&(grn_.Pttime), 
			&(grn_.Praybottom), 
			&isrclayer, 
			&vmint, 
			verbose );

		grn_.Ptakeoff = calc_takeoff_angle( 
			vmint.vp[isrclayer], 
			grn_.Prayparameter, 
			grn_.rdist, 
			grn_.evdp );

		  fprintf( stdout, " p=%6.4f rb=%6.1f tt=%7.2f toa=%5.1f\n",
                    grn_.Prayparameter,
                    grn_.Praybottom,
                    grn_.Pttime,
                    grn_.Ptakeoff );

	/***********************************************************/
	/*** write out greens functions as SAC files for testing ***/
	/***********************************************************/
		ista = 0;
		if( idump ) wrtgrn2sac( &grn_, ista );

	/**************************************************************/
	/*** write out binary greens functions library for glib2inv ***/
	/**************************************************************/
		fwrite(&grn_,sizeof(Greens),1,fp);
	}

/****************/
/*** clean up ***/
/****************/
	fclose(fp);
	free(z);
	printf("%s: program finished. Bye-Bye!\n\n\n", progname );
	exit(0);
}

void getparameters( int ac, char **av, Greens *g, 
		Depth_Info *depth_info, int *verbose, int *idump )
{
	double drdist,daz,dbaz;
	float zvec[3];
	float nyquist_frequency;
	char stadb_filename[128];

	int distaz( double, double, double, double, double *, double *, double * );
	int setpar( int, char ** );
	int getpar(), mstpar();
	void endpar();
	void create_mod( VelMod * );
	void print_mod0( VelMod * );
	void getsta( char *, Greens *, int * );

	
	setpar(ac,av);
/*************************/
/**** output settings ****/
/*************************/
	getpar("verbose", "b", verbose );
	getpar("dump", "b", idump );

/*****************************/
/**** set the depth range ****/
/*****************************/
	mstpar("zrange", "vf[3]", &zvec );
	depth_info->zmin = zvec[0];
	depth_info->zinc = zvec[1];
	depth_info->zmax = zvec[2];
	depth_info->nz   = ((depth_info->zmax-depth_info->zmin)/depth_info->zinc) + 1;

	if( *verbose )
	{
	  printf("%s: Compute Greens Functions for:\n", progname );
	  printf("%s: \t\t zmin=%g zinc=%g zmax=%g nz=%d\n",
	    progname, depth_info->zmin, depth_info->zinc, 
		depth_info->zmax, depth_info->nz );
	}

/************************************************************************/
/*** get the station location information, later get this information ***/
/*** from station geometry file including station and network         ***/
/************************************************************************/
	mstpar("stnm", "s", &(g->stnm) );
	mstpar("net",  "s", &(g->net) );
	mstpar("stadb", "s", stadb_filename );

/*************************************/
/*** source latitude and longitude ***/
/*************************************/
	mstpar("evla", "f", &(g->evla) );
	mstpar("evlo", "f", &(g->evlo) );

/*** important stuff ***/
	mstpar("nt",	"d", &(g->nt) );
	mstpar("dt",	"f", &(g->dt) );
	g->twin = g->nt * g->dt;

/*** f-k stuff ***/
	mstpar("eps",	"f", &(g->eps) );
	mstpar("smin",	"f", &(g->smin) );
	g->damp = 1.0;
	getpar("damp",  "f", &(g->damp) );
	g->kmax = 10000000;
	getpar("kmax",  "d", &(g->kmax) );

/**************************/
/*** reduction velocity ***/
/**************************/
	g->redv  = -1.;
	getpar("redv",  "f", &(g->redv) );
	g->t0 = 0;
	getpar("t0", "f", &(g->t0) );

/***************************************************************/
/*** leave off the .mod in filename append only when opening ***/
/***************************************************************/
	mstpar("modeldb", "s", &(g->v.modpath) );
	mstpar("velmod", "s", &(g->v.modfile) );

/******************************************************************************/
/*** default fmax is the nyquist frequency unless otherwise set to be lower ***/
/******************************************************************************/
	nyquist_frequency = 1/(2 * g->dt);
	getpar("fmax", "f", &(g->fmax) );
	if( g->fmax > nyquist_frequency ) g->fmax = nyquist_frequency;

	g->twin = g->nt * g->dt;

	endpar();

/***********************************/
/*** get the station information ***/
/***********************************/
	g->stel = 0;
	getsta( stadb_filename, g, verbose );
	printf("%s: stnm=%s net=%s stla=%g stlo=%g stel=%g evla=%g evlo=%g\n",
		progname, g->stnm, g->net, g->stla, g->stlo, g->stel,
		g->evla, g->evlo );

	if( distaz( (double)g->evla, (double)g->evlo, 
		(double)g->stla, (double)g->stlo, &drdist, &daz, &dbaz ) == 0 )
	{
		g->rdist = (float)drdist;
		g->az    = (float)daz;
		g->baz   = (float)dbaz;
		printf("%s: distaz: r=%g az=%g baz=%g\n", 
			progname, g->rdist, g->az, g->baz );
	}
	else
	{
		printf("%s: distaz: returned a fatal error\n", progname );
		exit(-1);
	}

/********************************************************************/
/*** use the reduction velocity (km/sec) to set the t0 time shift ***/
/********************************************************************/
	if ( g->redv > 0 ) 
	{
	  g->t0     = g->rdist/g->redv;
	  g->tstart = g->ts0 = g->t0;
	  g->tend   = g->tstart + g->twin;

	  printf("%s: dist=%g reduction velocity=%g(km/sec) t0=%g ",
		progname, g->rdist, g->redv, g->t0 );
	  printf("tstart=%g tend=%g twin=%g\n",
		g->tstart, g->tend, g->twin ); 
	}

/******************************************/
/*** set the green's function file name ***/
/******************************************/
	sprintf( g->filename, "%s.%s.%s.glib", g->stnm, g->net, g->v.modfile );
	if( *verbose )
        	printf("%s: Writing glib file to %s\n", progname, g->filename );

/*********************************************/
/*** set the velocity model data structure ***/
/*********************************************/
	create_mod( &(g->v) );
	/** if( *verbose ) print_mod0( &(g->v) ); **/
	print_mod0( &(g->v) );
}

void Usage_Print()
{
	fprintf(stderr, "%s: \n\t Version=%s \n\t Release Date=%s\n",
                        progname, Version_Label, Version_Date );

	fprintf(stderr, "\n USAGE: %s par=foo.par stnm= net=\n", progname );

	fprintf(stderr, "\n REQUIRED PARAMETERS:\n" );

	fprintf(stderr, "\t stnm=ELK              Station Code\n" );
	fprintf(stderr, "\t net=US                Network Affiliation Code\n" );

	fprintf(stderr, "\t velmod=wus            velocity model name without .mod extenstion\n" );
	fprintf(stderr, "\t modeldb=/mydir/models directory location of the model files\n" );
	fprintf(stderr, "\t stadb=/mydir/stations filename and path to the station location file\n" );
	fprintf(stderr, "\t zrange=2/2/10         starting depth, depth increment, ending depth\n" );
	fprintf(stderr, "\t evla=32.567           event latitude\n" );
	fprintf(stderr, "\t evlo=-120.456         event longitude\n" );

	fprintf(stderr, "\t dt=0.15               Green's function sample per seconds\n" );
	fprintf(stderr, "\t nt=2048               Green's function num pts (power of 2)\n" );
	fprintf(stderr, "\t fmax=0.9              Green's function max frequency (max is Nyquist)\n" );
	fprintf(stderr, "\t t0=0                  Green's function starting time 0=origin time\n" );
	fprintf(stderr, "\t rdev=10               Green's function reduction velocity (km/s)\n" );
	fprintf(stderr, "\t damp=1                Green's function damping (1 is OK)\n" );
	fprintf(stderr, "\t kmax=999999           Green's function maximum wavenumber\n" );
	fprintf(stderr, "\t eps=0.0001            Green's function accuracy tolorance parameter 1\n" );
	fprintf(stderr, "\t smin=0.0001           Green's function accuracy tolorance parameter 2\n" );
	
	fprintf(stderr, "\n OPTIONAL PARAMETERS: [DEFAULTS]\n" );
	fprintf(stderr, "\t [no]verbose           be verbosy [Default is off]\n" );
	fprintf(stderr, "\t [no]dump              write out the GF synthetics in SAC files\n" );
	fprintf(stderr, "\n\n" );
}
