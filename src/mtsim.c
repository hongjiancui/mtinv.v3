#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>
#include "../include/nrutil.h"     /** numerical recipes **/
#include "../include/mt.h"         /** global datatype and structure declarations **/

char progname[128];

int main(int ac, char **av)
{

/************************/
/*** event info stuff ***/
/************************/
	EventInfo *ev;
	int ista, nsta;
	char evinfo_filename[128];
	int it;
	
/************************/
/*** Greens stuff     ***/
/************************/
        Greens **grn;
	int nz,iz;
	float *z;

/****************/
/*** Solution ***/
/****************/
	Solution *sol;
	int iz_best;
	float bestfit, mechanism_size; 

/*******************/
/*** local stuff ***/
/*******************/
	int verbose = 0;	/** verbose   1=yes 0=no default no verbose ***/
	int idump = 0;
	int igmtmap = 0;
	int mtdegfree = 5;
	int forward = 0;
	float ts0 = 0;
	int check_on_status_ok;
	int FixMyiz;
	float FixMyZ = -99;
	int Distance_Normalize = 0;
	float DistNormR0 = 1; /*** default is R0 = 1 km  in R/R0 ****/

	int ishift = 0;
	float cortol = 1.0;
	float maxtimeshift = 0; /*** about 1 cycle 25 < maxtimeshift < 5 sec ***/

	FixISOZ myfixisoz;

	int PltXcorLabel = 1;
	int idumpxy = 0;
	int iorientation = LANDSCAPE; /* orientation for gmtwf.csh */
        char orientation[12];
	int iparallel = 1;

        int mysql_db_write = 0;
        int oracle_db_write = 0;
        int iuse_snr = 0;
        float MINSNR = 3;
	int igmt5 = 1;

/**************/
/*** output ***/
/**************/

        char asc_file_name[128];
	char ps_plot_filename[128];
	char pathname[128];
        FILE *fpasc;

/*****************************/
/**** bootstrap resampling ****/
/*****************************/

	int Nmodels;     	/*** number of bootstrap resamples           ***/
	float *res;       	/*** residual vector from best fit inversion ***/
	float *best_b_vector;	/*** synthetics from best fit inversion to   ***/
				/*** be combined with bootstrap residuals to ***/
				/*** generate new data vector populations    ***/
	int nres;		/*** number of points in res vector length   ***/

	/*** subroutine to calculate the number of rows in A matrix          ***/
	int size_A_matrix( EventInfo *, Greens **, int, int );
	
/*****************************/
/*** subroutine prototypes ***/
/*****************************/
	
	int setpar(int,char **);
	int getpar();
	int mstpar();
	void endpar();
	
/*** glib2inv_subs.c ***/
	EventInfo *glib2inv_get_input_parameters( char *, EventInfo *, int *, int );

/*** mtinv_subs.c ***/
	void load_the_data( EventInfo *, int, float, int );

/*** mtinv_subs.c ***/
	float *load_greens( EventInfo *, Greens **, int, int *, int );

/*** mtsim_subs.c ***/
	void invert0(	EventInfo *ev,
			Greens **grn,
			int nsta,
			int nz,
			int iz_best,
			Solution *sol,
			float *res,
			float *best_b_vector,
			int *nres, 
			int verbose,
			int mtdegfree,
			int Distance_Normalize, 
			float DistNormR0,
			FixISOZ myfixisoz );

/*** mtsim_invert_sim_parallel.c ***/
	void invert_sim_parallel(
			EventInfo *ev,
                        Greens **grn,
                        int nsta,
                        int iz_best,
                        Solution *sol,
                        float *res,
                        float *best_b_vector,
                        int Distance_Normalize,
                        float DistNormR0,
                        FixISOZ myfixisoz,
                        int Nmodels,
			int verbose );
		
/*** mtsim_subs.c ***/
	void invert_sim_serial(
			EventInfo *ev, 
			Greens **grn,
			int nsta,
			int iz_best, 
			Solution *sol,
			float *res, 
			float *best_b_vector,
			int Distance_Normalize,
			float DistNormR0,
			FixISOZ myfixisoz, 
			int Nmodels,
			int verbose );

/*** psplot.c ***/
	void psplot(    int nsta,
                        int iz,
                        char *filenbase,
                        EventInfo *ev,
                        Solution *sol,
                        Greens **grn,
                        int units,
                        int verbose,
                        int forward,
                        int PltXcorLabel,
                        int PS_Output_Screen );

/*** check_depths.c ***/
	void check_iso_depth( FixISOZ *myfixisoz, float *z, int nz, int verbose );

/*** check_depths.c ***/
	void check_depth( float FixMyZ, int *FixMyiz, float *z, int nz, int verbose );

/*** mtinv_subs.c ***/
	void compute_synthetics( int is, int iz, EventInfo *ev, Greens **grn, Solution *sol, int mtdegfree );

/*** mtinv_subs.c ***/
	void write_email_message( FILE *fp, int nsta, int iz, Solution *sol, EventInfo *ev, Greens **grn, int ifwd );

/*** mtinv_subs.c ***/
        void plot_z_vs_fit_gmt5( int iz, float *z, Solution *sol, EventInfo *ev );
        void plot_z_vs_fit_gmt4( int iz, float *z, Solution *sol, EventInfo *ev );

/*** mtinv_subs.c ***/
	void plotmech_gmt4( int iz, Solution *sol, EventInfo *ev, float mechanism_size );
	void plotmech_gmt5( int iz, Solution *sol, EventInfo *ev, float mechanism_size );

/*** mtinv_subs.c ***/
	void write_sol_rec( FILE *fp, int iz, int nsta, EventInfo *ev, Solution *sol, Greens **grn );

/*** mtinv_subs.c ***/
	void make_map_gmt5( int, int, EventInfo *, Solution *, Greens ** );
	void make_map_gmt4( int, int, EventInfo *, Solution *, Greens ** );

/*** crosscorrelation/cross_correlation.c ***/
        void xcorr(
                float *d, float *s,
                int nt, float dt, int *ilag, float *tlag,
                float *xcorcoef,
                float cortol, int verbose, int ishift );

/*** find_best_shift.c ***/
        void find_best_shift( EventInfo *ev, float cortol, float maxtimeshift, float *time_shift_all );

/*** mtinv_subs.c ***/
        void time_shift( EventInfo *ev, int nsta, int verbose );

/*** mtinv_subs.c ***/
        void plot_results_gmt4( int iz, Solution *sol, EventInfo *ev );
	void plot_results_gmt5( int iz, Solution *sol, EventInfo *ev );

/*** mtinv_subs.c ***/
        void calc_azi_gap_and_dmin( int iz, int nsta, EventInfo *ev, Solution *sol, Greens **g );

/*** mtinv_subs.c ***/
        void write_mterror( FILE *fp, int nz, Solution *sol, EventInfo *ev, Greens **g );

/*** sacextrema/sacextrema.c ***/
        void sac_absmax( float *x, int n, float *absmax );

	void wrtnewsac( char *FO, float dt, int ns, float *ar, float b );

/*** mtsim.c ***/
	void Usage_Print( void );
        
/*** mtinv_subs.c ***/
        void gmtwfplot( EventInfo *ev, Solution *sol, Greens **grn, int nsta, int iz, int iorientation, int verbose );
       
/*** mtinv_subs.c ***/
        void dumpxy( EventInfo *ev, Solution *sol, Greens **grn, int nsta, int iz, int verbose );

/*** shorten_path.c ***/
        char *shorten_path( char *pathname, char *filename );

/********************/
/*** program name ***/
/********************/

	strcpy( pathname, av[0] );

	shorten_path( pathname, progname );

	if( verbose )
        {
          fprintf( stdout, "%s: STDOUT: Version=%s ReleaseDate=%s exec full path=%s\n",
                progname, Version_Label, Version_Date, pathname );
        }

        fprintf( stderr, "%s: STDERR: Version=%s ReleaseDate=%s exec full path=%s\n",
                progname, Version_Label, Version_Date, pathname );

/*************/
/*** USAGE ***/
/*************/
	if( ac <= 1 ) Usage_Print();

/*****************************************************/
/*** get the input parameters foreach each station ***/
/*****************************************************/
	setpar(ac,av);
	mstpar( "par", "s", &evinfo_filename );
	mstpar( "mtdegfree", "d", &mtdegfree );
	mstpar( "ts0", "f", &ts0 );
	mstpar( "fixz", "f", &FixMyZ );
	mstpar( "Nboot", "d", &Nmodels );
	getpar( "parallel", "b", &iparallel );

	getpar( "verbose", "b", &verbose );
	getpar( "norm", "b", &Distance_Normalize );
	if( Distance_Normalize )
        {
                mstpar( "R0", "f", &DistNormR0 );
        }

	getpar( "PltXcorLabel", "b", &PltXcorLabel );

        myfixisoz.z = 0;
        getpar( "FixISOZ", "f", &(myfixisoz.z) );

	getpar( "shift", "b", &ishift );
        if(ishift)
        {
                mstpar( "ctol", "f", &cortol );
                mstpar( "maxshift", "f", &maxtimeshift );
        }

        getpar( "use_snr", "b", &iuse_snr );
        if( iuse_snr )
        {
                mstpar( "minsnr", "f", &MINSNR );
        }

	getpar( "gmt5", "b", &igmt5 );

	endpar();

/*********************/
/**** START MAIN  ****/
/*********************/

	if(verbose)
        {
                fprintf(stderr, "%s: verbose ON\n", progname );
                fprintf(stdout, "%s: verbose ON\n", progname );
        }
        else
        {
                fprintf(stderr, "%s: verbose OFF\n", progname );
        }

        if(iuse_snr)
        {
                if(verbose)
                fprintf(stdout, "%s: iuse_snr ON minsnr=%g\n", progname, MINSNR );
                fprintf(stderr, "%s: iuse_snr ON minsnr=%g\n", progname, MINSNR );
        }
        else
        {
                if(verbose)
                fprintf(stdout, "%s: iuse_snr OFF\n", progname );
                fprintf(stderr, "%s: iuse_snr OFF\n", progname );
        }

        if(verbose)
        {
                if(mtdegfree==1)fprintf(stdout, "%s: mtdegfree=%d EXPLOSION\n", progname, mtdegfree );
                if(mtdegfree==5)fprintf(stdout, "%s: mtdegfree=%d DEVIATORIC\n", progname, mtdegfree );
                if(mtdegfree==6)fprintf(stdout, "%s: mtdegfree=%d FULL_MT\n", progname, mtdegfree );
        }

        if(mtdegfree==1)fprintf(stderr, "%s: mtdegfree=%d EXPLOSION\n", progname, mtdegfree );
        if(mtdegfree==5)fprintf(stderr, "%s: mtdegfree=%d DEVIATORIC\n", progname, mtdegfree );
        if(mtdegfree==6)fprintf(stderr, "%s: mtdegfree=%d FULL_MT\n", progname, mtdegfree );

        if( myfixisoz.z > 0 )
                myfixisoz.iswitch = 1;
        else
                myfixisoz.iswitch = 0;

	if( myfixisoz.iswitch )
        {
                if(verbose)
                fprintf( stdout,
                        "%s: FixISOZ z = %g iswitch= %d is ON\n",
                        progname,
                        myfixisoz.z,
                        myfixisoz.iswitch  );
                                                                                                                                                                                                         
                fprintf( stderr,
                        "%s: FixISOZ z = %g iswitch= %d is ON\n",
                        progname,
                        myfixisoz.z,
                        myfixisoz.iswitch  );
        }
        else
        {
                if(verbose)
                fprintf( stdout, "%s: FixISOZ=%d is OFF\n",
                        progname, myfixisoz.iswitch );
                fprintf( stderr, "%s: FixISOZ=%d is OFF\n",
                        progname, myfixisoz.iswitch );
        }

        if( Distance_Normalize )
        {
                if(verbose)
                fprintf( stdout, "%s: Distance_Normalize is ON R0=%g km\n", progname, DistNormR0 );
                fprintf( stderr, "%s: Distance_Normalize is ON R0=%g km\n", progname, DistNormR0 );
        }
 
	if(ishift)
        {
                if(verbose)
                {
                  fprintf( stdout,
                    "%s: Shiftng by max cross-correlation is ON cortol=%g maxtimeshift=%g\n",
                        progname, cortol, maxtimeshift );
                }
                fprintf( stderr,
                  "%s: Shiftng by max cross-correlation is ON cortol=%g maxtimeshift=%g\n",
                        progname, cortol, maxtimeshift );
        }
        else
        {
                if(verbose)
                fprintf( stdout, "%s: Shiftng by max cross-correlation is OFF\n", progname );
                fprintf( stderr, "%s: Shiftng by max cross-correlation is OFF\n", progname );
        }

/******************************************/
/*** allocate memory for parameter list ***/
/******************************************/
	if( verbose )
        {
                fprintf( stdout, "%s: allocating memory for data and event information\n",
                        progname );
        }

	ev  = (EventInfo *)malloc(sizeof(EventInfo));
	ev  = (EventInfo *)glib2inv_get_input_parameters( evinfo_filename, ev, &nsta, verbose );

/********************************************/
/*** loop over stations and load the data ***/
/********************************************/

	if( verbose )
        {
                fprintf( stdout, "%s: mtsim(): reading data nsta=%d\n", progname, nsta );
        }

        fprintf( stderr, "%s: mtsim(): reading data nsta=%d\n", progname, nsta );

        load_the_data( ev, nsta, ts0, verbose );

        for( ista = 0 ; ista < nsta ; ista++ )
        {
                  fprintf( stdout, "%s: mtsim(): STDOUT: ista=%d nt=%d %d %d dt=%g %g %g\n",
                        progname,
                        ista,
                        ev[ista].ew.s.npts,
                        ev[ista].ns.s.npts,
                        ev[ista].z.s.npts,
                        ev[ista].ew.s.delta,
                        ev[ista].ns.s.delta,
                        ev[ista].z.s.delta );
        /*** debug ***/
                for( it = 0; it < ev[ista].ew.s.npts; it++ )
                {
                        if( ev[ista].ew.data[it] > 0.1 )
                        {
                                fprintf( stdout, "--------------------------------\n" );
                                fprintf( stdout, "ista=%d it=%d ev[ista].ew.data[it]=%g\n",
                                        ista, it, ev[ista].ew.data[it] );
                        }
                }
                /* writesacfile( &ev[ista] ); */
        }

/******************************/
/*** load green's functions ***/
/******************************/

	if(verbose) 
	{
	   printf( "%s: allocating memory for Green's function\n",
		progname );
	}

	grn = (Greens **)malloc( nsta*sizeof(Greens *) );

	z = (float *)load_greens( ev, grn, nsta, &nz, verbose );

/*************************************/
/*** FixISOZ check depth set index ***/
/*************************************/
	
	if( myfixisoz.iswitch )
	{
		check_iso_depth( &myfixisoz, z, nz, verbose );
	}

/**************************************/
/*** check if fixing solution depth ***/
/*** iz_best get reset below        ***/
/**************************************/

	if( FixMyZ != -99 )
	{
		check_depth( FixMyZ, &FixMyiz, z, nz, verbose );
	}

/***********************************************************************************/
/*** turn off stations that have all three channels/components with SNR < MINSNR ***/
/*** turn off only if iused was already turned off/ do not override              ***/
/***********************a***********************************************************/
	
	if(verbose)
        {
          fprintf( stdout, "%s: iuse_snr = %d minsnr = %g\n",
                progname, iuse_snr, MINSNR );
        }

        if( iuse_snr )
        {
          for( ista=0; ista<nsta; ista++ )
          {
                if( ( ev[ista].z.P2P_snr  < MINSNR) &&
                    ( ev[ista].ns.P2P_snr < MINSNR) &&
                    ( ev[ista].ew.P2P_snr < MINSNR) &&
                    ( ev[ista].iused == 1 ) )
                {
                        ev[ista].iused = 0;
                }
                /*** else iused=0, leave alone ***/
          }
        }

/*****************************************************/
/*** Error check if any any stations are turned on ***/
/*****************************************************/

	check_on_status_ok = 0;
	for( ista=0; ista<nsta; ista++ )
	{
		if( ev[ista].iused == 1 ) check_on_status_ok=1;
	}

	if(verbose)
        fprintf( stdout, "%s: check_on_status_ok=%d\n", progname, check_on_status_ok );

	if( !check_on_status_ok )
	{
		if(verbose)
                fprintf( stdout,
                  "%s: ERROR! no stations turned on in the par file %s\n",
                        progname, evinfo_filename );
        
                fprintf( stderr,
                  "%s: ERROR! no stations turned on in the par file %s\n",
                        progname, evinfo_filename );
                exit(-1);
	}


/*********************************/
/*** set the type of inversion ***/
/*********************************/

	sol = (Solution *)malloc(nz*sizeof(Solution));

	for( iz=0; iz<nz; iz++ )
	{
		if( mtdegfree == 1 )  sol[iz].mt_type = EXPLOSION;
		if( mtdegfree == 5 )  sol[iz].mt_type = DEVIATORIC;
		if( mtdegfree == 6 )  sol[iz].mt_type = FULL_MT;
	}

	if( verbose ) 
	{
		fprintf( stdout, "%s: calling invert()\n", progname );
	}

	iz = 0;
	nres = size_A_matrix( ev, grn, nsta, iz );
	res = calloc( nres, sizeof(float) );
	best_b_vector = calloc( nres+1, sizeof(float) );

	invert0( ev, grn, nsta, nz, FixMyiz, sol, res, best_b_vector, &nres, verbose, 
		mtdegfree, Distance_Normalize, DistNormR0, myfixisoz );

	wrtnewsac( "res.sac", 1.0, nres, res, 0.0 );
	wrtnewsac( "best_b_vector.sac", 1.0, nres, best_b_vector, 0.0 );

/****************************/
/*** what is the best fit ***/
/****************************/

	if( FixMyZ != -99 )
	{
		iz_best = FixMyiz;
	} 
	else
	{
		iz_best = 0;
		bestfit = sol[0].total_fitness1;
		for( iz=0; iz<nz; iz++ )
		{
			if( sol[iz].total_fitness1 > bestfit ) 
			{
				bestfit = sol[iz].total_fitness1;
				iz_best = iz;
			}
		}
	}

/******************************************************************/
/*** write out a gmt shell script for each time shift and depth ***/
/******************************************************************/

	mechanism_size = 0.4;
	
	if( igmt5 )
        {
                if(verbose)
                {
                        fprintf( stdout, "%s: mtinv(): calling plotmech_gmt5():\n", progname );
                }
                plotmech_gmt5( iz_best, sol, ev, mechanism_size );
        }
        else
        {
                if(verbose)
                {
                        fprintf( stdout, "%s: mtinv(): calling plotmech_gmt4():\n", progname );
                }
                plotmech_gmt4( iz_best, sol, ev, mechanism_size );
        }

/********************************************/
/*** make an ascii plot for fast email    ***/
/********************************************/

	if(verbose)
          fprintf( stdout, "%s: mtsim(): calling calc_azi_gap_and_dmin():\n", progname );

	calc_azi_gap_and_dmin( iz_best, nsta, ev, sol, grn );

	sprintf( asc_file_name, "email_T%05.1fsec_Z%05.1fkm_.txt", 
			sol[iz_best].ot, grn[0][iz_best].evdp );

	if(verbose)
        {
                fprintf(stdout, "%s: mtsim.c: Writing Email Messages to %s\n",
                        progname, asc_file_name );
        }

	if( (fpasc=fopen( asc_file_name, "w" ) ) == NULL )
		fprintf(stderr, "cannot open file for writting\n");

	if(verbose)
          fprintf( stdout, "%s: mtsim(): calling write_email_message():\n", progname );

	write_email_message( fpasc, nsta, iz_best, sol, ev, grn, forward );

	fclose(fpasc);

/*************************************/
/*** write out a gmt shell script  ***/
/*************************************/

	if(verbose)
        {
                fprintf( stdout, "%s: mtsim.c: bestfit=%g iz=%d\n",
                        progname, bestfit, iz_best );
        }

	sprintf( asc_file_name, "results.%d.out", mtdegfree );

	if( (fpasc=fopen( asc_file_name, "a" ) ) == NULL )
        {
                fprintf(stdout, "cannot open file for writting\n");
        }

        if(verbose)
	{
                fprintf( stdout, "%s: mtsim(): calling write_sol_rec():\n", progname );
	}

        write_sol_rec( fpasc, iz_best, nsta, ev, sol, grn );
        fclose(fpasc);

/*************************************/
/*** write out a gmt shell script ***/
/*************************************/

	if(igmt5)
        {
                if(verbose)
                {
                  fprintf( stdout, "%s: mtinv(): calling plot_z_vs_fit_gmt5(): \n",
                        progname );
                }
                plot_z_vs_fit_gmt5( iz_best, z, sol, ev );
        }
        else
        {
                if(verbose)
                {
                  fprintf( stdout, "%s: mtinv(): calling plot_z_vs_fit_gmt4(): \n",
                        progname );
                }
                plot_z_vs_fit_gmt4( iz_best, z, sol, ev );
        }

/***************************************/
/*** compute synthetics for plotting ***/
/***************************************/

	if(verbose)
	{
	  fprintf( stdout, "%s: recalculate the synthetics for iz_best=%d\n",
		progname, iz_best );
	}

        for( ista=0; ista<nsta; ista++ )
        {
                ev[ista].syn_r.data = calloc( ev[ista].nt, sizeof(float) );
                ev[ista].syn_z.data = calloc( ev[ista].nt, sizeof(float) );
                ev[ista].syn_t.data = calloc( ev[ista].nt, sizeof(float) );

                if(verbose)
                {
                  fprintf( stdout,
                    "%s: mtsim.c: calling compute_synthetics iz_best=%d ista=%d nt=%d dt=%g\n",
                        progname, iz_best, ista, ev[ista].nt, ev[ista].dt );
                }

                compute_synthetics( ista, iz_best, ev, grn, sol, mtdegfree );
        }

/*******************************************************************************************/
/*** for each station do a cross correlation to find the lag time and correlation values ***/
/*******************************************************************************************/

	if(verbose)
	{
	  fprintf( stdout,
            "%s: mtsim.c: calling cross correlation for iz_best=%d cortol=%f maxshift=%f ishift=%d\n",
                progname, iz_best, cortol, maxtimeshift, ishift );
	}

	for( ista=0; ista<nsta; ista++ )
	{
		/*** force ishift = 0 ***/
		xcorr( &(ev[ista].z.data[0]), &(ev[ista].syn_z.data[0]), ev[ista].nt, ev[ista].dt,
                        &(ev[ista].izlag), &(ev[ista].ztlag), &(ev[ista].zxcor), cortol, verbose, 0 );
       
                xcorr( &(ev[ista].ns.data[0]), &(ev[ista].syn_r.data[0]), ev[ista].nt, ev[ista].dt,
                        &(ev[ista].irlag), &(ev[ista].rtlag), &(ev[ista].rxcor), cortol, verbose, 0 );
        
                xcorr( &(ev[ista].ew.data[0]), &(ev[ista].syn_t.data[0]), ev[ista].nt, ev[ista].dt,
                        &(ev[ista].itlag), &(ev[ista].ttlag), &(ev[ista].txcor), cortol, verbose, 0 );

	/*** based on max cross correlation shift all data by same amount ***/
                
		if(ishift)
		{
			find_best_shift( &(ev[ista]), cortol, maxtimeshift, &(ev[ista].time_shift_all) );
		} 

		if(verbose)
                {
                  fprintf( stdout, "%s: %s.%s.%s: otshift=%g stashift=%g\n",
                        progname,
                        ev[ista].ew.s.kstnm,
                        ev[ista].ew.s.knetwk,
                        ev[ista].ew.s.kcmpnm,
                        ev[ista].ts0,
                        ev[ista].time_shift_all );
                                                                                                                                                                                                         
                  fprintf( stdout, "%s: \t izlag=%3d ztlag=%5.2f zxcor=%.2f\n",
                        progname,
                        ev[ista].izlag,
                        ev[ista].ztlag,
                        ev[ista].zxcor );
                                                                                                                                                                                                         
                  fprintf( stdout, "%s: \t irlag=%3d rtlag=%5.2f rxcor=%.2f\n",
                        progname,
                        ev[ista].irlag,
                        ev[ista].rtlag,
                        ev[ista].rxcor );
                                                                                                                                                                                                         
                  fprintf( stdout, "%s: \t itlag=%3d ttlag=%5.2f txcor=%.2f\n",
                        progname,
                        ev[ista].itlag,
                        ev[ista].ttlag,
                        ev[ista].txcor );
                }

	} /*** loop over stations and cross correlate data/syn ***/

/*** do the time shift, and cross correlate again to update the plot ***/

        if(ishift)
        {
                if(verbose)
                {
                  fprintf( stdout, "%s: mtsim(): calling time_shift()\n", progname );
                }
    
                time_shift( ev, nsta, verbose );
     
                for( ista=0; ista<nsta; ista++ )
                {
                 /*** force ishift = 0 ***/
                 xcorr( &(ev[ista].z.data[0]), &(ev[ista].syn_z.data[0]), ev[ista].nt, ev[ista].dt,
                        &(ev[ista].izlag), &(ev[ista].ztlag), &(ev[ista].zxcor), cortol, verbose, 0 );
      
                 xcorr( &(ev[ista].ns.data[0]), &(ev[ista].syn_r.data[0]), ev[ista].nt, ev[ista].dt,
                        &(ev[ista].irlag), &(ev[ista].rtlag), &(ev[ista].rxcor), cortol, verbose, 0 );
       
                 xcorr( &(ev[ista].ew.data[0]), &(ev[ista].syn_t.data[0]), ev[ista].nt, ev[ista].dt,
                        &(ev[ista].itlag), &(ev[ista].ttlag), &(ev[ista].txcor), cortol, verbose, 0 );
                }
        }

/*************************************/
/*** Cgraphics PS library routines ***/
/*************************************/

	sprintf( ps_plot_filename, "plot_T%05.1fsec_Z%05.1fkm_", 
		sol[iz_best].ot, grn[0][iz_best].evdp );

	if(verbose)
        {
                fprintf( stdout, "%s: mtsim.c: mtsim(): Plotting Postscript Plot %s\n",
                        progname, ps_plot_filename );
                                                                                                                                                                                                         
                fprintf( stdout, "%s: mtsim.c: mtsim(): calling psplot()\n", progname );
        }

	psplot( nsta, iz_best, ps_plot_filename, ev, sol, grn, 0, verbose, 
		forward, PltXcorLabel, LANDSCAPE );

/**********************/
/*** make a GMT map ***/
/**********************/
/************************************** removed **********************************************/
/***
	if( igmtmap )
        {
          if( igmt5 )
          {
                if(verbose)
                {
                        fprintf(stdout,
                          "%s: mtinv.c: mtinv(): calling make_map_gmt5(), making a GMT map for plotting\n",
                                progname );
                }
                make_map_gmt5( iz_best, nsta, ev, sol, grn );
          }
          else
          {
                if(verbose)
                {
                        fprintf(stdout,
                          "%s: mtinv.c: mtinv(): calling make_map_gmt4(), making a GMT map for plotting\n",
                                progname );
                }
                make_map_gmt4( iz_best, nsta, ev, sol, grn );
          }
        }
***/
/************************************** removed **********************************************/

/**********************/
/*** plot results   ***/
/**********************/
	if( igmt5 )
        {
                if(verbose)
                {
                  fprintf(stdout, "%s: mtinv.c: mtinv(): calling plot_results_gmt5()\n",
                        progname );
                }
                plot_results_gmt5( iz_best, sol, ev );
        }
        else
        {
                if(verbose)
                {
                  fprintf(stdout, "%s: mtinv.c: mtinv(): calling plot_results_gmt4()\n",
                        progname );
                }
                plot_results_gmt4( iz_best, sol, ev );
        }

	free(sol);

/*************************************/
/*** bootstrap resample simulation ***/
/*************************************/

	sol = (Solution *)malloc(nz*sizeof(Solution));
        for( iz=0; iz<nz; iz++ )
        {
           if( mtdegfree == 1 )  sol[iz].mt_type = EXPLOSION;
           if( mtdegfree == 5 )  sol[iz].mt_type = DEVIATORIC;
           if( mtdegfree == 6 )  sol[iz].mt_type = FULL_MT;
        }

	if(iparallel)
	{

	  fprintf( stdout, "%s: calling invert_sim_parallel()\n", progname );
	  fprintf( stderr, "%s: calling invert_sim_parallel()\n", progname );

	  invert_sim_parallel(
			ev,
			grn,
			nsta, 
			iz_best,
			sol,
			res, 
			best_b_vector,
			Distance_Normalize,
			DistNormR0, 
			myfixisoz, 
			Nmodels,
			verbose );
	}
	else
	{

	  fprintf( stdout, "%s: calling invert_sim_serial()\n", progname );
	  fprintf( stderr, "%s: calling invert_sim_serial()\n", progname );

          invert_sim_serial( ev,
			grn,
			nsta,
			iz_best,
			sol,
			res,
			best_b_vector,
			Distance_Normalize,
			DistNormR0,
			myfixisoz,
			Nmodels,
			verbose );	
	}

/*****************************/
/*** unallocate the memory ***/
/*****************************/
	if(verbose) fprintf(stderr, "%s: Trying to free memory...", progname );

	free(sol);
	free(z);

	for( ista=0; ista<nsta; ista++ )
	{
		free( ev[ista].ew.data );
		free( ev[ista].ns.data );
		free( ev[ista].z.data );

		free( ev[ista].syn_t.data );
		free( ev[ista].syn_r.data );
		free( ev[ista].syn_z.data );
	}
	free(ev);
	free(grn);

	free(res);
	free(best_b_vector);

        if(verbose) fprintf(stdout, " Done.\n" );

        if(verbose)
        fprintf(stdout, "%s: mtsim.c: mtsim(): STDOUT: Program finished.  Bye-Bye!\n", progname );

        fprintf(stderr, "%s: mtsim.c: mtsim(): STDERR: Program finished.  Bye-Bye!\n", progname );

        exit(0);

} /*** end of mtsim.c ***/
