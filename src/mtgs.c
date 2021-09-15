#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>
#include "../include/nrutil.h"
#include "../include/mt.h"

char progname[128];

int main( int ac, char **av )
{
/* see include/mt.h Mo = math.pow( 10.0, 1.5*(Mw+10.73) ) = 1.2445146117713818e+16; Reference Mw = 0.0 */

/******************************************/
/*** event information                  ***/
/******************************************/
	EventInfo *ev;
	int ista, nsta;
	char evinfo_filename[128];

/******************************************/
/*** Green's funciton                   ***/
/******************************************/

	Greens **grn;
	FILE *fpg;
	int nz, iz;
	float *z;

/******************************************/
/*** solution                           ***/
/******************************************/

	Solution *sol;
	int iz_best;
	float bestfit;

/******************************************/
/*** local stuff                        ***/
/******************************************/
	int verbose = 0;
	int idump = 0;
	int mtdegfree = 5;
	int forward = 0;
	float ts0 = 0;
	int ienvelope = 0;

	int check_on_status_ok;
	int FixMyiz;
	float FixMyZ = -99;
	int Distance_Normalize = 0;
	int ishift = 0;
	float cortol = 1.0;
	FixISOZ myfixisoz;

/******************************************/
/*** function prototypes                ***/
/******************************************/
	int setpar(int,char **);
	int getpar();
	int mstpar();
	void endpar();

	EventInfo *glib2inv_get_input_parameters( char *, EventInfo *, int *, int );
	void Usage_Print();
	void load_the_data( EventInfo *, int, float, int );
	float *load_greens( EventInfo *, Greens **, int, int *, int );
	void check_iso_depth( FixISOZ *, float *, int, int );
        void check_depth( float, int *, float *, int, int );
	void gs( EventInfo *, Greens **, int, int, float *, Solution *, int, int, int, int, FixISOZ, int );

	void psplot( int, int, char *, EventInfo *, Solution *, Greens **, int, int, int );
        void compute_synthetics( int, int, EventInfo *, Greens **, Solution *, int );
        void write_email_message( FILE *, int, int, Solution *, EventInfo *, Greens **, int );
        void plot_z_vs_fit( int, float *, Solution *, EventInfo *, char * );
        void write_sol_rec( FILE *, int, int, EventInfo *, Solution *, Greens ** );
        void make_gmt_map( int, int, EventInfo *, Solution *, Greens ** );
        void xcorr( float *, float *, int, float, int *, float *, float *, float, int, int );
        void write_mterror( FILE *, int, Solution *, EventInfo *, Greens ** );
        void sac_absmax( float *, int, float * );
	void envelope( float *, int, float );

/******************************************/
/*** output                             ***/
/******************************************/

	char asc_file_name[128], ps_plot_filename[128];
	FILE *fpasc;

/******************************************/
/*** program name                       ***/
/******************************************/

	strcpy( progname, av[0] );
	
	fprintf( stderr, "%s: version=%s release date=%s\n",
		progname, Version_Label, Version_Date );

	fprintf( stdout, "%s: version=%s release date=%s\n",
		progname, Version_Label, Version_Date );

/**************/
/*** usage ***/
/**************/
	if( ac <= 1 ) Usage_Print();

	setpar(ac,av);
	mstpar( "par",       "s", &evinfo_filename );
	mstpar( "mtdegfree", "d", &mtdegfree );

	getpar( "fixz",      "f", &FixMyZ );
	getpar( "verbose",   "b", &verbose );
        getpar( "dumpsac",   "b", &idump );
	getpar( "norm",      "b", &Distance_Normalize );
	getpar( "shift",     "b", &ishift );
        getpar( "ctol",      "f", &cortol );
        getpar( "FixISOZ",   "f", &(myfixisoz.z) );
	getpar( "env",       "b", &ienvelope );
	endpar();
	
	myfixisoz.z = 0;
	if( myfixisoz.z > 0 )
		myfixisoz.iswitch = 1;
	else
		myfixisoz.iswitch = 0;

/******************************************/
/*** allocate memory for parameter list ***/
/******************************************/

	ev  = (EventInfo *)malloc(sizeof(EventInfo));
	ev  = (EventInfo *)glib2inv_get_input_parameters( evinfo_filename, ev, &nsta, verbose );

/********************************************/
/*** loop over stations and load the data ***/
/********************************************/

	if( verbose ) printf( "%s: reading data nsta=%d\n", progname, nsta );
	load_the_data( ev, nsta, ts0, verbose );

	fprintf( stdout, "done with loading data \n");

/******************************/
/*** load green's functions ***/
/******************************/

	grn = (Greens **)malloc( nsta*sizeof(Greens *) );
	z = (float *)load_greens( ev, grn, nsta, &nz, verbose );

	fprintf( stdout, "done with loading Green's functions \n" );

/*************************************/
/*** FixISOZ check depth set index ***/
/*** fix the isotropic greens      ***/
/*** function depth to shallow     ***/
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

/*****************************************************/
/*** Error check if any any stations are turned on ***/
/*****************************************************/

        check_on_status_ok = 0;
        for( ista=0; ista<nsta; ista++ )
        {
                if( ev[ista].iused == 1 ) check_on_status_ok=1;
        }
        if( !check_on_status_ok )
        {
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
           fprintf(stdout, "%s: calling gs mtdegfree=%d\n", progname, mtdegfree );

	fprintf(stdout, "nz=%d verbose=%d idump=%d mtdegfree=%d Distance_Normalize=%d\n",
		nz, verbose, idump, mtdegfree, Distance_Normalize );

	gs( ev, grn, nsta, nz, z, sol, verbose, idump, mtdegfree, Distance_Normalize, myfixisoz, ienvelope );

	fprintf(stdout, "%s: done with gs\n", progname );
	fflush(stdout);
	
	iz_best = 0;
	bestfit = sol[0].total_fitness1;
	for( iz=0; iz<nz; iz++ )
	{
		if( ev[0].evdp == z[iz] )
		{
			if(verbose)printf("best iz=%d z=evdp=%g\n", iz, z[iz] );
			iz_best = iz;
			bestfit = sol[iz].total_fitness1;
		}
		if( FixMyZ != -99 ) iz_best = FixMyiz;
	}

/******************************************************************/
/*** write out for R. B. Herrmann mteig and mtdec input         ***/
/*** overwrites previous version of of mteig.in                 ***/
/******************************************************************/
        fpasc = fopen( "mteig.in", "w" );
        fprintf( fpasc, "3\n" );
        fprintf( fpasc, "%g %g %g\n",
                sol[iz_best].moment_tensor[1][1], sol[iz_best].moment_tensor[1][2], sol[iz_best].moment_tensor[1][3] );
        fprintf( fpasc, "%g %g %g\n",
                sol[iz_best].moment_tensor[2][1], sol[iz_best].moment_tensor[2][2], sol[iz_best].moment_tensor[2][3] );
        fprintf( fpasc, "%g %g %g\n",
                sol[iz_best].moment_tensor[3][1], sol[iz_best].moment_tensor[3][2], sol[iz_best].moment_tensor[3][3] );
        fclose(fpasc);


/********************************************/
/*** write an output file with the errors ***/
/********************************************/
        sprintf( asc_file_name, "mterror.out" );
        if( (fpasc=fopen( asc_file_name, "w" )) == NULL )
                fprintf( stderr, "cannot open file %s for writting\n", asc_file_name );
        write_mterror( fpasc, nz, sol, ev, grn );
        fclose(fpasc);

/***************************************/
/*** compute synthetics for plotting ***/
/***************************************/

	for( ista=0; ista<nsta; ista++ )
        {
                ev[ista].syn_r.data = (float *)malloc( ev[ista].nt * sizeof(float) );
                ev[ista].syn_z.data = (float *)malloc( ev[ista].nt * sizeof(float) );
                ev[ista].syn_t.data = (float *)malloc( ev[ista].nt * sizeof(float) );

		if( ienvelope )
		{
		  envelope( ev[ista].z.data, ev[ista].nt, ev[ista].dt );
		  envelope( ev[ista].ns.data, ev[ista].nt, ev[ista].dt );
		  envelope( ev[ista].ew.data, ev[ista].nt, ev[ista].dt );
		}

                if( verbose )
                {
                    printf( "calling compute_synthetics iz_best=%d ista=%d nt=%d dt=%g\n",
                         iz_best, ista, ev[ista].nt, ev[ista].dt );
                }
		compute_synthetics( ista, iz_best, ev, grn, sol, mtdegfree );

		if( ienvelope )
		{
		  envelope( ev[ista].syn_r.data, ev[ista].nt, ev[ista].dt );
                  envelope( ev[ista].syn_z.data, ev[ista].nt, ev[ista].dt );
                  envelope( ev[ista].syn_t.data, ev[ista].nt, ev[ista].dt );
		}
	}

/*******************************************************************************************/
/*** for each station do a cross correlation to find the lag time and correlation values ***/
/*******************************************************************************************/
	if(verbose)
        {
                fprintf( stderr, "%s: cross correlatio for iz_best=%d cortol=%f ishift=%d\n",
                        progname, iz_best, cortol, ishift );
        }

	for( ista=0; ista<nsta; ista++ )
        {
                xcorr( &(ev[ista].z.data[0]), &(ev[ista].syn_z.data[0]), ev[ista].nt, ev[ista].dt,
                        &(ev[ista].izlag), &(ev[ista].ztlag), &(ev[ista].zxcor), cortol, verbose, ishift );

                xcorr( &(ev[ista].ns.data[0]), &(ev[ista].syn_r.data[0]), ev[ista].nt, ev[ista].dt,
                        &(ev[ista].irlag), &(ev[ista].rtlag), &(ev[ista].rxcor), cortol, verbose, ishift );

                xcorr( &(ev[ista].z.data[0]), &(ev[ista].syn_z.data[0]), ev[ista].nt, ev[ista].dt,
                        &(ev[ista].itlag), &(ev[ista].ttlag), &(ev[ista].txcor), cortol, verbose, ishift );
        }

/*************************************/
/*** Cgraphics PS library routines ***/
/*************************************/

        sprintf( ps_plot_filename, "plot_T%05.1fsec_Z%05.1fkm_",
                sol[iz_best].ot, grn[0][iz_best].evdp );

        if(verbose)
        {
                fprintf(stdout, "%s: Plotting Postscript Plot %s\n",
                        progname, ps_plot_filename );
        }
        psplot( nsta, iz_best, ps_plot_filename, ev, sol, grn, 0, verbose, forward );

/*****************************/
/*** unallocate the memory ***/
/*****************************/
        fprintf(stderr, "%s: Trying to free memory...", progname );

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

        fprintf(stderr, " Done.\n" );
        fprintf(stderr, "%s: Program finished.  Bye-Bye!\n", progname );

        exit(0);

} /*** end of mtgs.c ***/

/****************************************************************************************************/
/****************************************************************************************************/
/******************************/
/*** grid search subroutine ***/
/******************************/
/****************************************************************************************************/
/****************************************************************************************************/

void gs( EventInfo *ev, Greens **grn, int nsta, int nz, float *z, Solution *sol, int verbose, int idump, int mtdegfree,
        int Distance_Normalize, FixISOZ myfixisoz, int ienvelope )
{
/* see include/mt.h Mo = math.pow( 10.0, 1.5*(Mw+10.73) ) = 1.2445146117713818e+16; Reference Mw = 0.0 */

	int ista, iz, myiz;
	int notfound = 1;
	float **mt;
	MomentTensor Ma, Mn;

	int i, j, it, nt, irow, rows, icol, cols;
	float **a_matrix;   /*** A matrux with dimensions a[1..rows][1..cols] ***/
        float **u_matrix;   /*** temp space               u[1..rows][1..cols] ***/
        float **v_matrix;   /*** temp space               v[1..cols][1..cols] ***/
        float *w_vector;    /*** temp space               w[1..cols]          ***/
        float *b_vector;    /*** data                     b[1..rows]          ***/
        float *s_vector;    /*** synthetic                s[1..rows]          ***/
        float *x_vector;    /*** solution (mom ten)       x[1..cols]          ***/
	float *e_vector;    /*** error vector             e[1..cols]          ***/
	float **cv_matrix;  /*** covariance matrix       cv[1..cols][1..cols] ***/
	float *xtmp;

	SDRvector *sdr;
	int isdr, nsdr, MAXSDR;
	int best_iz, best_isdr;
	float best_vred, best_Mw;
	float wmin,wmax;
	FILE *fpgs;
	char outputfile[256];

	float logerrMdc, logerrMclvd, logerrMiso;
	float tmp1,tmp2,tmp3,tmpa,tmpb,tmpc;
/***************************/
/*** function prototypes ***/
/***************************/

	float **matrix( int, int, int, int );
        float *vector( int, int );

	void SDR_init( SDRvector *, int * );
        void set_moment_tensor( MomentTensor *, MomentTensor *, float *, int, int );
        void normalize_moment_tensor( MomentTensor *, MomentTensor *, int );
        void sdr_to_mt( float **, float, float, float, float, int );
        float variance_reduction( float *, float *, int, int );
        float compute_l2norm_error( float *, float *, int );
        void mt2eig( MomentTensor, Solution *, int, int );
        void eig2iso( Solution *, int, int );
        void Eig2MajorDC( Solution *, int, int );
        void Eig2MinorDC( Solution *, int, int );

	void initalize_workspace( int, int, float **, float **, float *, float *, float **, float *, float * );
	void make_amatrix3( Greens **, EventInfo *, int, int, float **, float *, int, int, FixISOZ, int );

	void svbksb( float **, float *, float **, int, int, float *, float * );
        void svdcmp( float **, int, int, float *, float ** );
        void svdvar( float **, int, float *, float ** );

	void compute_error( int, int, float **,  float *, float *, float *, float *, float *, float *, float * );

        void matmul( int, float **, int, float *, int, float * );
        void grn2disp( Greens *, EventInfo *, int, int );
	void grn2grn( Greens **, EventInfo *, int, int, int, int );
        void sac_minmax( float *, int, float *, float *, float * );

/******************************************/
/*** check the depths for each station  ***/
/******************************************/

        for( ista=0; ista<nsta; ista++ )
        {
                for( iz=0; iz<nz; iz++ )
                {
                        if( ev[ista].my_z == z[iz] && ev[ista].evdp == z[iz] )
                        {
                                notfound = 0;
                                myiz = iz;
                        }
                }
		iz = myiz;

                if( notfound )
                {
                  fprintf( stderr, "%s: ista=%d matching depth not found for station %s.%s ... looking for my_z=%g ",
                        progname, ista, ev[ista].stnm, ev[ista].net, ev[ista].my_z );
                  for( iz=0; iz<nz; iz++ )
                  {
                        fprintf( stderr, "\t iz=%d z=%18.7f\n", iz, z[iz] );
                  }
                  exit(0);
                }
                else
                {
                  if(verbose)
                  {
                     fprintf( stdout, "%s: ista=%d matching depth z=%g iz=%d found for station=%s.%s\n",
                        progname, ista, ev[ista].my_z, iz, ev[ista].stnm, ev[ista].net );
                  }
                }
        }

/*******************************************/
/*** initalize strike dip rake structure ***/
/*******************************************/
	MAXSDR = 200000;
	sdr = (SDRvector *)calloc( MAXSDR, sizeof(SDRvector) );
	SDR_init( &sdr[0], &nsdr );
	if( nsdr > MAXSDR ) 
	{	
		fprintf( stderr, "%s: reset nsdr > %d\n", progname, MAXSDR );
		exit(-1);
	}

	if(verbose) fprintf( stdout, "%s: nsdr = %d\n", progname, nsdr );
	fprintf( stderr, "%s: nsdr = %d\n", progname, nsdr );

/*************************************/
/*** allocate and initalize memory ***/
/*************************************/
	iz = myiz;

	if( verbose )
	{
	  fprintf(stdout, "%s: initalizing matrix for forward calculation nsta=%d iz=%d\n",
		progname, nsta, iz );
	}

	cols = 3;
	rows = 1;
	for( ista=0; ista < nsta; ista++ )
	{
		if( ev[ista].iused == 1 )
		{
			rows += 3 * grn[ista][iz].nt;
		}
	}

	if(verbose) 
		fprintf( stdout, "%s: rows=%d cols=%d\n", progname, rows, cols );

	if( verbose )
             fprintf(stdout, "%s: Allocating memory for iz=%d\n", progname, iz );

	if(verbose)fprintf(stdout, "%s: Allocating memory for mt\n", progname );
	mt = matrix( 1, cols+1, 1, cols+1 );
	
	if(verbose)fprintf(stdout, "%s: Allocating memory for a_matrix\n", progname );
	a_matrix  = matrix( 1, rows+1, 1, cols+1 );

	if(verbose)fprintf(stdout, "%s: Allocating memory for u_matrix\n", progname );
	u_matrix  = matrix( 1, rows+1, 1, cols+1 );
	
	if(verbose)fprintf(stdout, "%s: Allocating memory for v_matrix\n", progname );
	v_matrix  = matrix( 1, cols+1, 1, cols+1 );

	if(verbose)fprintf(stdout, "%s: Allocating memory for x_vector\n", progname );
	x_vector = vector( 1, cols+1 );

	if(verbose)fprintf(stdout, "%s: Allocating memory for w_vector\n", progname );
	w_vector = vector( 1, cols+1 );

	if(verbose)fprintf(stdout, "%s: Allocating memory for b_vector\n", progname );
	b_vector = vector( 1, rows+1 );

	if(verbose)fprintf(stdout, "%s: Allocating memory for s_vector\n", progname );
	s_vector = vector( 1, rows+1 );

	if(verbose)fprintf(stdout, "%s: initializing memory\n", progname );
	initalize_workspace( rows, cols, u_matrix, a_matrix, w_vector, x_vector, v_matrix, b_vector, s_vector );

/***********************************************/
/*** loop over all depths, sdr, and stations ***/
/***********************************************/
	best_vred = -9.99999E-13;
	best_Mw = 0;
	best_isdr = 1;
	best_iz = 0;
	iz = myiz;
	sprintf( outputfile, "mtgs.out.%g", ev[0].ts0 );
	fpgs = fopen( outputfile, "w" );

	if(verbose) {
	  fprintf(stdout,"%s: begin GS writting grid search results to file %s\n",
			progname, outputfile );
	}

	for( isdr=0; isdr<nsdr; isdr++ )
	{
		sdr[isdr].iz  = myiz;
		sdr[isdr].z   = z[iz];
		/* sdr[isdr].Mw  = ev[0].Mw; */
		
        	for( ista=0; ista<nsta; ista++ )
        	{
			ev[ista].str = sdr[isdr].s;
			ev[ista].dip = sdr[isdr].d;
			ev[ista].rak = sdr[isdr].r;
			ev[ista].Mw = 0.0;
			ev[ista].Mo = base_moment;

			/***
			  printf( "%s: ista=%d %s.%s z=%g str=%g dip=%g rak=%g Mw=%g mtdegfree=%d\n",
				progname, ista, 
				ev[ista].stnm,
				ev[ista].net,
				ev[ista].my_z,
				ev[ista].str,
				ev[ista].dip,
				ev[ista].rak,
				ev[ista].Mw,
				mtdegfree );
			****/
			grn2grn( grn, ev, ista, iz, verbose, mtdegfree );
		}

		make_amatrix3( grn, ev, nsta, iz, a_matrix, b_vector, mtdegfree, Distance_Normalize, myfixisoz, ienvelope );

	/*** do the inversion here ***/

		for( j=1; j<=cols; j++ )
		{
			for( i=1; i<=rows; i++ )
			{
				u_matrix[i][j] = a_matrix[i][j];
			}
		}
	
		svdcmp( u_matrix, rows, cols, w_vector, v_matrix );
		wmax = 0;
		for( j=1; j<=cols; j++)
		{
			if( w_vector[j] > wmax) wmax=w_vector[j];
		}
		wmin = wmax * 1.0E-5;
		for( j=1; j<=cols; j++)
		{
			if( w_vector[j] < wmin ) w_vector[j]=0.0;
		}
		svbksb( u_matrix, w_vector, v_matrix, rows, cols, b_vector, x_vector );

	/***************************************************/
	/*** SKIP ALL SOLUTIONS WITH NEGATIVE DC MOMENTS ***/
	/***************************************************/
		if( x_vector[1] < 0 )
		{
			x_vector[1] = fabs(x_vector[1]);
			sdr[isdr].r = -1 * sdr[isdr].r;
		}

		matmul( 0, a_matrix, cols, x_vector, rows, s_vector );
		sol[iz].var_red = variance_reduction( b_vector, s_vector, 1, rows+1 );
		sdr[isdr].vred = sol[iz].var_red;
		sol[iz].ot   = ev[0].ot_shift;

		sdr[isdr].Mdc   = x_vector[1] * base_moment;
                sdr[isdr].Mclvd = x_vector[2] * base_moment;
                sdr[isdr].Miso  = x_vector[3] * base_moment;

                sdr[isdr].Mtotal = fabs(sdr[isdr].Mdc) + fabs(sdr[isdr].Mclvd) + fabs(sdr[isdr].Miso);
                sdr[isdr].Mw = log10( sdr[isdr].Mtotal )/1.5 - 10.73;

		sdr[isdr].piso  = fabs( sdr[isdr].Miso  )/ sdr[isdr].Mtotal;
		sdr[isdr].pclvd = fabs( sdr[isdr].Mclvd )/ sdr[isdr].Mtotal;
		sdr[isdr].pdc   = fabs( sdr[isdr].Mdc   )/ sdr[isdr].Mtotal;

	       /***           isdr Mdc   Mclvd Miso  Mtot  Mw    S     D     Rak   Vred  Pdc   Pclvd Piso  ot    shift  ***/
		fprintf(fpgs, "%06d %5.2e %5.2e %5.2e %5.2e %5.2f %3.0f %3.0f %4.0f %7.3f %5.3f %5.3f %5.3f %4.1f %4.1f %4.1f %4.1f\n",
			isdr,
			fabs(sdr[isdr].Mdc),
			fabs(sdr[isdr].Mclvd),
			fabs(sdr[isdr].Miso),
			sdr[isdr].Mtotal,
			sdr[isdr].Mw,
			sdr[isdr].s, 
			sdr[isdr].d,
			sdr[isdr].r,
			sdr[isdr].vred,
			sdr[isdr].pdc,
			sdr[isdr].pclvd,
			sdr[isdr].piso,
			ev[0].ot_orig.fsec,
			sol[iz].ot,
			ev[0].ts0,
			z[myiz] );
		/* if( sdr[isdr].vred >= best_vred && sdr[isdr].Mdc > 0 && sdr[isdr].pclvd < 0.3 ) */
		/* if( sdr[isdr].vred >= best_vred && sdr[isdr].Mw < best_Mw ) */

		if( sdr[isdr].vred >= best_vred && sdr[isdr].pclvd < 0.3 )
		{
			best_Mw   = sdr[isdr].Mw;
			best_vred = sol[iz].var_red;
			best_isdr = isdr;
	 	}

	} /*** loop over isdr ***/

	fprintf( stdout, "%s: best_vred = %g best_isdr = %d Mdc=%5.2e Mclvd=%5.2e Miso=%5.2e Mtotal=%5.2e ",
		progname,
                sdr[best_isdr].vred,
                best_isdr,
                sdr[best_isdr].Mdc,
                sdr[best_isdr].Mclvd,
                sdr[best_isdr].Miso,
                sdr[best_isdr].Mtotal );

	fprintf( stdout, "Mw=%4.2f SDR=%g/%g/%g PDC=%3.0f PCLVD=%3.0f PISO=%3.0f\n", 
                sdr[best_isdr].Mw,
                sdr[best_isdr].s,
                sdr[best_isdr].d,
                sdr[best_isdr].r,
                sdr[best_isdr].pdc*100,
                sdr[best_isdr].pclvd*100,
                sdr[best_isdr].piso*100 );

	if(verbose)
          fprintf( stdout, "%s: freeing memory inside mtgs.c:gs()\n",
                progname );
	
	if(verbose)fprintf( stdout, "%s: freeing memory a_matrix\n", progname );
	free_matrix( a_matrix, 1, rows+1, 1, cols+1 );

	if(verbose)fprintf( stdout, "%s: freeing memory u_matrix\n", progname );
	free_matrix( u_matrix, 1, rows+1, 1, cols+1 );

	if(verbose)fprintf( stdout, "%s: freeing memory v_matrix\n", progname );
	free_matrix( v_matrix, 1, cols+1, 1, cols+1 );

	if(verbose)fprintf( stdout, "%s: freeing memory x_vector\n", progname );
	free_vector( x_vector, 1, cols+1 );

	if(verbose)fprintf( stdout, "%s: freeing memory w_vector\n", progname );
	free_vector( w_vector, 1, cols+1 );

	if(verbose)fprintf( stdout, "%s: freeing memory s_vector\n", progname );
	free_vector( s_vector, 1, rows+1 );

	if(verbose)fprintf( stdout, "%s: freeing memory b_vector\n", progname );
	free_vector( b_vector, 1, rows+1 );

	if(verbose)
	  fprintf( stdout, "%s: alllocating memory to compute best response mtgs.c:gs() rows=%d cols=%d\n",
		progname, rows, cols );

	if(verbose)fprintf(stdout, "%s: Allocating memory for u_matrix\n", progname );
	u_matrix  = matrix( 1, rows+1, 1, cols+1 );

	if(verbose)fprintf(stdout, "%s: Allocating memory for v_matrix\n", progname );	
	v_matrix  = matrix( 1, cols+1, 1, cols+1 );

	if(verbose)fprintf(stdout, "%s: Allocating memory for w_vector\n", progname );
	w_vector = vector( 1, cols+1 );

	if(verbose)fprintf( stdout, "%s: allocating memory for a_matrix\n", progname );
	a_matrix = matrix( 1, rows+1, 1, cols+1 );

	if(verbose)fprintf( stdout, "%s: allocating memory for x_vector\n", progname );
	x_vector = vector( 1, cols+1 );

	if(verbose)fprintf( stdout, "%s: allocating memory for b_vector\n", progname );
	b_vector = vector( 1, rows+1 );

	if(verbose)fprintf( stdout, "%s: allocating memory for s_vector\n", progname );
	s_vector = vector( 1, rows+1 );

	initalize_workspace( rows, cols, u_matrix, a_matrix, w_vector, x_vector, v_matrix, b_vector, s_vector );

/***********************************************/
/*** recompute the response for the best fit ***/
/***********************************************/

	isdr = best_isdr;
	iz = sdr[best_isdr].iz; 
	iz = myiz;

/***********************************************/
/*** make DC moment tensor from SDR          ***/
/***********************************************/
	ista = 0;
	ev[ista].Mw = log10(fabs(sdr[best_isdr].Mdc))/1.5 - 10.73;
	ev[ista].Mo = pow( 10.0, 1.5*( ev[ista].Mw + 10.73 ) );
	sdr_to_mt( mt, sdr[best_isdr].s, sdr[best_isdr].d, sdr[best_isdr].r, ev[ista].Mw, verbose );

/*** add the isotropic and CLVD to DC moment tensor ****/

	sol[iz].moment_tensor[1][1] = mt[1][1] + (sdr[best_isdr].Miso/base_moment) - ((1.0*(sdr[best_isdr].Mclvd))/base_moment);
	sol[iz].moment_tensor[1][2] = mt[1][2];
	sol[iz].moment_tensor[1][3] = mt[1][3];
	sol[iz].moment_tensor[2][1] = mt[2][1];
	sol[iz].moment_tensor[2][2] = mt[2][2] + (sdr[best_isdr].Miso/base_moment) - ((1.0*(sdr[best_isdr].Mclvd))/base_moment);
	sol[iz].moment_tensor[2][3] = mt[2][3];
	sol[iz].moment_tensor[3][1] = mt[3][1];
	sol[iz].moment_tensor[3][2] = mt[3][2];
	sol[iz].moment_tensor[3][3] = mt[3][3] + (sdr[best_isdr].Miso/base_moment) + ((2.0*(sdr[best_isdr].Mclvd))/base_moment);
	
	xtmp = vector( 1, 7 );
	xtmp[1] = sol[iz].moment_tensor[1][1];
	xtmp[2] = sol[iz].moment_tensor[2][2];
	xtmp[3] = sol[iz].moment_tensor[1][2];
	xtmp[4] = sol[iz].moment_tensor[1][3];
	xtmp[5] = sol[iz].moment_tensor[2][3];
	xtmp[6] = sol[iz].moment_tensor[3][3];
	set_moment_tensor( &Ma, &Mn, xtmp, mtdegfree, verbose );
	free_vector( xtmp, 1, 7 );

	normalize_moment_tensor( &Ma, &Mn, verbose );
	sol[iz].dmoment  = Ma.moment;
	sol[iz].mw       = Ma.Mw;
	sol[iz].exponent = Ma.expon;
	sol[iz].abcassa  = Ma.abcassa;
	sol[iz].mrr = Mn.rr;
	sol[iz].mtt = Mn.tt;
	sol[iz].mff = Mn.ff;
	sol[iz].mrt = Mn.rt;
	sol[iz].mrf = Mn.rf;
	sol[iz].mtf = Mn.tf;
	sol[iz].mxx = Mn.xx;
	sol[iz].mxy = Mn.xy;
	sol[iz].mxz = Mn.xz;
	sol[iz].myy = Mn.yy;
	sol[iz].myz = Mn.yz;
	sol[iz].mzz = Mn.zz;

/********************************************************/
/*** calculate the eigenvalues from the moment tensor ***/
/********************************************************/
	mt2eig( Mn, sol, iz, verbose );
	eig2iso( sol, iz, verbose );
	Eig2MajorDC( sol, iz, verbose );
	Eig2MinorDC( sol, iz, verbose ); 

/****************************************************/
/*** make the A matrix do the forward calculation ***/
/*** then compute the variance reduction          ***/
/****************************************************/

	for( ista=0; ista<nsta; ista++ )
	{
		ev[ista].str = sdr[best_isdr].s;
		ev[ista].dip = sdr[best_isdr].d;
		ev[ista].rak = sdr[best_isdr].r;
		ev[ista].Mw = 0.0;
		ev[ista].Mo = pow( 10.0, (1.5*10.73) );
		grn2grn( grn, ev, ista, iz, verbose, mtdegfree );
	}

	make_amatrix3( grn, ev, nsta, iz, a_matrix, b_vector, mtdegfree, Distance_Normalize, myfixisoz, ienvelope );

	for( j=1; j<=cols; j++ )
	{
		for( i=1; i<=rows; i++ )
		{
			u_matrix[i][j] = a_matrix[i][j];
		}
	}
	svdcmp( u_matrix, rows, cols, w_vector, v_matrix );
	wmax = 0;
	for( j=1; j<=cols; j++) 
	{
		if( w_vector[j] > wmax) wmax=w_vector[j];
	}
	wmin = wmax * 1.0E-5;
	for( j=1; j<=cols; j++)
	{
		if( w_vector[j] < wmin ) w_vector[j]=0.0;
	}
	svbksb( u_matrix, w_vector, v_matrix, rows, cols, b_vector, x_vector);

	compute_error( rows, cols, v_matrix, w_vector, x_vector, s_vector, b_vector, &logerrMdc, &logerrMclvd, &logerrMiso );

	sdr[best_isdr].Mdc   = x_vector[1] * base_moment;
	sdr[best_isdr].Mclvd = x_vector[2] * base_moment;
	sdr[best_isdr].Miso  = x_vector[3] * base_moment;

/**
	x_vector[1] = sdr[best_isdr].Mdc   / base_moment;
	x_vector[2] = sdr[best_isdr].Mclvd / base_moment;
	x_vector[3] = sdr[best_isdr].Miso  / base_moment;
***/
	matmul( 0, a_matrix, cols, x_vector, rows, s_vector );

	sol[iz].var_red = variance_reduction( b_vector, s_vector, 1, rows+1 );
	printf( "%s: iz=%d %%VRED=%g\n", progname, iz, sol[iz].var_red );

	sol[iz].l2norm_error = compute_l2norm_error( b_vector, s_vector, rows );
	printf( "%s: iz=%d %%L2NORM=%g\n", progname, iz,  sol[iz].l2norm_error );
/*
	tmp1 = pow( 10, log10(fabs(sdr[best_isdr].Mdc)) - logerrMdc );
	tmp2 = pow( 10, log10(fabs(sdr[best_isdr].Mdc)) + logerrMdc );
	tmp3 = pow( 10, log10(fabs(sdr[best_isdr].Mclvd)) - logerrMclvd );
	tmpa = pow( 10, log10(fabs(sdr[best_isdr].Mclvd)) + logerrMclvd );
	tmpb = pow( 10, log10(fabs(sdr[best_isdr].Miso)) - logerrMiso );
	tmpc = pow( 10, log10(fabs(sdr[best_isdr].Miso)) + logerrMiso );
*/
	/**           isdr Mdc  Mclvd Miso  Mtot  Mw    s     d     r    vred   pdc   pclvd piso  Mdc_err Mclvd_err Miso_err    ****/

	fprintf(fpgs, "%d %5.2e %5.2e %5.2e %5.2e %5.2f %3.0f %3.0f %4.0f %7.3f %5.3f %5.3f %5.3f %5.2e %5.2e %5.2e %g %g %g %g %4d%02d%02d %02d%02d\n",
                best_isdr,
                fabs(sdr[best_isdr].Mdc),
                fabs(sdr[best_isdr].Mclvd),
                fabs(sdr[best_isdr].Miso),
                sdr[best_isdr].Mtotal,
                sdr[best_isdr].Mw,
                sdr[best_isdr].s,
                sdr[best_isdr].d,
                sdr[best_isdr].r,
                sdr[best_isdr].vred,
                sdr[best_isdr].pdc,
                sdr[best_isdr].pclvd,
                sdr[best_isdr].piso,
                pow(10,logerrMdc),
		pow(10,logerrMclvd),
		pow(10,logerrMiso),
		ev[0].ot_orig.fsec,
		sol[iz].ot,
		ev[0].ts0,
		z[myiz],
		ev[0].ot.year, ev[0].ot.month, ev[0].ot.mday, ev[0].ot.hour, ev[0].ot.min
	);

	fflush(fpgs);
	fflush(stdout);
        fclose(fpgs);

	sol[iz].total_fitness1 = sol[iz].var_red;
	sol[iz].total_fitness2 = sol[iz].var_red;
	sol[iz].evlo = ev[0].evlo;
	sol[iz].evla = ev[0].evla;
	sol[iz].evdp = ev[0].evdp;
	sol[iz].ot   = ev[0].ot_shift;

        if(verbose) 
	  fprintf( stdout, "%s: freeing memory inside mtgs.c:gs()\n",
               	progname );

	if(verbose) fprintf( stdout, "%s: freeing memory for a_matrix\n", progname );
        free_matrix( a_matrix, 1, rows+1, 1, cols+1 );

	if(verbose) fprintf( stdout, "%s: freeing memory u_matrix\n", progname );
        free_matrix( u_matrix, 1, rows+1, 1, cols+1 );

        if(verbose) fprintf( stdout, "%s: freeing memory v_matrix\n", progname );
        free_matrix( v_matrix, 1, cols+1, 1, cols+1 );

        if(verbose) fprintf( stdout, "%s: freeing memory x_vector\n", progname );
        free_vector( x_vector, 1, cols+1 );

        if(verbose) fprintf( stdout, "%s: freeing memory w_vector\n", progname );
        free_vector( w_vector, 1, cols+1 );

	if(verbose) fprintf( stdout, "%s: freeing memory for s_vector\n", progname );
        free_vector( s_vector, 1, rows+1 );

	if(verbose) fprintf( stdout, "%s: freeing memory for b_vector\n", progname );
        free_vector( b_vector, 1, rows+1 );

       	if(verbose) fprintf(stdout, "%s: leaving mtgs.c:gs()\n", progname );
}

void Usage_Print()
{
	fprintf( stderr, "\nUSAGE: %s par= mtdegfree=(1,5,6)\n", progname );
	fprintf( stderr, "\t [no]verbose [no]dumpsac " );
        fprintf( stderr, "\t ts0=[0] fixz=[-99] [no]norm [no]shift ctol=[1] FixISOZ=[-99]\n" );

	fprintf( stderr, "\nREQUIRED PARAMETERS:\n" );
        fprintf( stderr, "\t par=(glib2inv.par) station parameter file\n" );
        fprintf( stderr, "\t mtdegfree=(1,5,6) 1=Isotropic MT, 5=Deviatoric MT, 6=Full MT\n" );

        fprintf( stderr, "\n OPTIONAL PARAMETERSL [DEFAULTS]\n" );
        fprintf( stderr, "\t ts0=[0]       Origin Time Shift Default is [0]\n" );
        fprintf( stderr, "\t fixz=[-99]    fix the depth Default is [-99] which turns off option\n" );
        fprintf( stderr, "\t [no]verbose   give verbose print to stdout. Default is off.\n" );
        fprintf( stderr, "\t [no]dumpsac   write out data and synthetics as SAC files. default is off\n" );
        fprintf( stderr, "\t [no]norm     distance normalization default is off\n" );
        fprintf( stderr, "\t [no]shift    shift the data automatically by cross correlation peak. default is off\n" );
        fprintf( stderr, "\t ctol=[0..1]  Correlation coefficient tolerance to shift the data when coef > ctol. defaut is off\n" );
        fprintf( stderr, "\t FixISOZ=     fix the depth of the rex and zex Green's function.  Default is off\n" );
        fprintf( stderr, "\n\n" );
}

void SDR_init( SDRvector *sdr, int *nsdr )
{
	float str,dip,rak;
	float dstr,ddip,drak;
	int isdr;
	int test;
	int ntot,nstr,ndip,nrak;

	dstr = 10; 
	ddip = 10; 
	drak = 10;
 
	nstr = 360/dstr;
	ndip = 90/ddip + 1;
	nrak = 180/drak + 1;
	ntot = nstr * ndip * nrak;
	fprintf(stderr, "nstr=%d ndip=%d nrak=%d ntot=%d\n", nstr, ndip, nrak, ntot );

	isdr = 0;
	for( str = 0; str < 360; str += dstr )
	{
		for( dip = 15; dip <= 90; dip += ddip )
		{
			for( rak = -90; rak <= 90; rak += drak )
			{
				sdr[isdr].s = str;
				sdr[isdr].d = dip;
				sdr[isdr].r = rak;
				sdr[isdr].id = isdr;
				isdr++;
			}
		}
	}
	
	test = 1;
	if( test )
	{
		sdr[isdr].s = 254;
		sdr[isdr].d = 82;
		sdr[isdr].r = -20;
		sdr[isdr].id = isdr;
		isdr++;
	}

	*nsdr = isdr;
}


void grn2grn( Greens **grn, EventInfo *ev, int ista, int iz, int verbose, int mtdegfree )
{

/*** ten greens functions ***/
	float *rss, *rds, *rdd, *rep, *zss, *zds, *zdd, *zep, *tss, *tds;
	float *tra, *rad, *ver;

/*** directional cosine coefficients ***/
	float a1, a2, a3, a4, a5;
	Tensor M; /* float Mxx, Myy, Mzz, Mxy, Mxz, Myz; */
	float strr, dipr, rakr;
	float half=0.5;
/* see include/mt.h Mo = math.pow( 10.0, 1.5*(Mw+10.73) ) = 1.2445146117713818e+16; Reference Mw = 0.0 */

	int nt, it;
	float pi, d2r, dt, t0, e, tt, tr, fi, r, dist, azimuth, area;

/*** set some constants ***/
        pi  = M_PI;
        d2r = pi/180.;
        nt = grn[ista][iz].nt;
        dt = grn[ista][iz].dt;
        t0 = grn[ista][iz].t0;
        e  = t0 + (nt*dt);

/*** assign pointers ***/
	rss = grn[ista][iz].g.rss;
        rds = grn[ista][iz].g.rds;
        rdd = grn[ista][iz].g.rdd;
        rep = grn[ista][iz].g.rep;
        zss = grn[ista][iz].g.zss;
        zds = grn[ista][iz].g.zds;
        zdd = grn[ista][iz].g.zdd;
        zep = grn[ista][iz].g.zep;
        tss = grn[ista][iz].g.tss;
        tds = grn[ista][iz].g.tds;
	
	fi = grn[ista][iz].az * d2r;

	/***
	ev->Mo = pow( 10.0, (1.5*( ev->Mw + 10.73)) );
	if(verbose)fprintf( stdout, "Mo=%e Mw=%g\n", ev->Mo, ev->Mw );
	****/

	strr = ev->str * d2r;
        dipr = ev->dip * d2r;
        rakr = ev->rak * d2r;

	M.xx = -( sin(dipr) * cos(rakr) * sin(2*strr) + sin(2*dipr) * sin(rakr) * sin(strr)*sin(strr) );
        M.yy =  ( sin(dipr) * cos(rakr) * sin(2*strr) - sin(2*dipr) * sin(rakr) * cos(strr)*cos(strr) );
        M.zz =  ( sin(2*dipr) * sin( rakr ) );
        M.xy =  ( sin(dipr) * cos(rakr) * cos(2*strr) + 0.5*sin(2*dipr) * sin(rakr) * sin(2*strr) );
        M.xz = -( cos(dipr) * cos(rakr) * cos(strr) + cos(2*dipr) * sin(rakr) * sin(strr) );
        M.yz = -( cos(dipr) * cos(rakr) * sin(strr) - cos(2*dipr) * sin(rakr) * cos(strr) );
	
/*** compute the coefficients ****/

        a1 = half * ( M.xx - M.yy ) * cos( 2 * fi ) + M.xy * sin( 2 * fi );
        a2 = M.xz * cos( fi ) + M.yz * sin( fi );
        a3 = -half*( M.xx + M.yy );
        a4 = half * ( M.xx - M.yy ) * sin( 2 * fi ) - M.xy * cos( 2 * fi );
        a5 = -M.yz * cos( fi ) + M.xz * sin( fi );

	for( it=0; it<nt; it++)
	{
	  grn[ista][iz].ver[it] = ( a1 * zss[it] + a2 * zds[it] + a3 * zdd[it] ); /*** vertical   ***/
          grn[ista][iz].rad[it] = ( a1 * rss[it] + a2 * rds[it] + a3 * rdd[it] ); /*** radial     ***/
          grn[ista][iz].tra[it] = ( a4 * tss[it] + a5 * tds[it] );                /*** transverse ***/
	}
}

void make_amatrix3( Greens **grn, EventInfo *ev, int nsta, int iz, float **a_matrix, float *b_vector,
        int mtdegfree, int Distance_Normalize, FixISOZ myfixisoz, int ienvelope )
{
	int it, nt, ista, irow;

	float *rss, *rds, *rdd, *rep;
	float *zss, *zds, *zdd, *zep;
	float *tss, *tds;
	float *rad, *tra, *ver;

	float R0 = 1;
	float Rnorm = 1;
	int ziso;
	float dt;

	float *etra, *erad, *ever;
	float *erdd, *erep, *ezdd, *ezep; 
	float *edatz, *edatr, *edatt;

/*****
	float etra[4096], erad[4096], ever[4096], erdd[4096];
	float erep[4096], ezdd[4096], ezep[4096];
	float edatz[4096], edatr[4096], edatt[4096];
****/

	void envelope( float *, int, float );

	irow = 1;
	for( ista = 0; ista < nsta; ista++ )
	{
		if( ev[ista].iused == 0 )  /*** skip station ***/
		{
			/*
			fprintf( stderr, "skipping station: ista=%d stnm=%s net=%s\n", 
				ista, ev[ista].stnm, ev[ista].net );
			*/
			continue;
		}

		nt = grn[ista][iz].nt;
		dt = grn[ista][iz].dt;

		if( Distance_Normalize )
		{
			Rnorm = grn[ista][iz].rdist/R0;
		}
		
		if( myfixisoz.iswitch )
			ziso = myfixisoz.indexz;
		else
			ziso = iz;

	/************************************************/
	/*** assignment of pointers                   ***/
	/************************************************/

		rss = grn[ista][iz].g.rss;
		rds = grn[ista][iz].g.rds;
		rdd = grn[ista][iz].g.rdd;
		rep = grn[ista][ziso].g.rep;
		zss = grn[ista][iz].g.zss;
		zds = grn[ista][iz].g.zds;
		zdd = grn[ista][iz].g.zdd;
		zep = grn[ista][ziso].g.zep;
		tss = grn[ista][iz].g.tss;
		tds = grn[ista][iz].g.tds;
		rad = grn[ista][iz].rad; 
		tra = grn[ista][iz].tra;
		ver = grn[ista][iz].ver;

		erad = calloc( nt, sizeof(float) );
		etra = calloc( nt, sizeof(float) );
		ever = calloc( nt, sizeof(float) );
		erdd = calloc( nt, sizeof(float) );
		erep = calloc( nt, sizeof(float) );
		ezep = calloc( nt, sizeof(float) );
		ezdd = calloc( nt, sizeof(float) );
		edatz = calloc( nt, sizeof(float) );
		edatr = calloc( nt, sizeof(float) );
		edatt = calloc( nt, sizeof(float) );

		for( it=0; it<nt; it++ )
		{
			erad[it] = rad[it];
			etra[it] = tra[it];
			ever[it] = ver[it];
			erdd[it] = rdd[it];
			erep[it] = rep[it];
			ezep[it] = zep[it];
			ezdd[it] = zdd[it];
			edatz[it] = ev[ista].z.data[it];
			edatr[it] = ev[ista].ns.data[it];
			edatt[it] = ev[ista].ew.data[it];
		}

	/************************************************/
	/*** transverse/tangential component          ***/
	/************************************************/
		if( ienvelope )
		{
			envelope( etra, nt, dt );
			envelope( edatt, nt, dt );
		}

		for( it = 0; it < nt; it++ )
		{
			a_matrix[irow][1] = etra[it] * Rnorm;   /*** Mdc ***/
			a_matrix[irow][2] = 0;                  /*** Mclvd ***/
			a_matrix[irow][3] = 0;                  /*** Miso ***/
			b_vector[irow]    = edatt[it] * Rnorm;  /*** Data ***/
			irow++;
		}
	
	/************************************************/
	/*** radial components                       ***/
	/************************************************/
		if( ienvelope )
		{
			envelope( erad, nt, dt );
			envelope( erdd, nt, dt );
			envelope( erep, nt, dt );
			envelope( edatr, nt, dt );
		}

		for( it = 0; it < nt; it++ )
		{
			a_matrix[irow][1] = erad[it] * Rnorm; 
			a_matrix[irow][2] = erdd[it] * Rnorm;
			a_matrix[irow][3] = erep[it] * Rnorm; 
			b_vector[irow]    = edatr[it] * Rnorm;
			irow++;
		}

	/************************************************/
	/*** vertical components                      ***/
	/************************************************/
		if( ienvelope )
		{
			envelope( ever, nt, dt );
			envelope( ezdd, nt, dt );
			envelope( ezep, nt, dt );
			envelope( edatz, nt, dt );
		}

		for( it = 0; it < nt; it++ )
		{
			a_matrix[irow][1] = ever[it] * Rnorm;
			a_matrix[irow][2] = ezdd[it] * Rnorm;  
			a_matrix[irow][3] = ezep[it] * Rnorm; 
			b_vector[irow]    = edatz[it] * Rnorm;
			irow++;
		}

		free( erad );
		free( etra );
		free( ever );
		free( erdd );
		free( erep );
		free( ezep );
		free( ezdd );
		free( edatz );
		free( edatr );
		free( edatt );

	} /*** loop over ista ***/

} /*** end of make_amatrix3 ***/


/********************************/
/***** initialize memory ********/
/********************************/

void initalize_workspace( int rows, int cols, float **u_matrix, float **a_matrix, float *w_vector,
	float *x_vector, float **v_matrix, float *b_vector, float *s_vector )
{
	int i, j;
	for( j=1; j<=cols; j++ )
	{
		for( i=1; i<=rows; i++ )
		{
			u_matrix[i][j] = 0;
			a_matrix[i][j] = 0;
		}
		w_vector[j] = 0;
		x_vector[j] = 0;
	}
	for( j=1; j<=cols; j++ )
	{
		for( i=1; i<=cols; i++ )
		{
			v_matrix[i][j] = 0;
		}
	}
	for( i=1; i<=rows; i++ )
	{
		b_vector[i] = 0;
		s_vector[i] = 0;
	}
}


/****************************************************/
/*** calculate the covariance matrix              ***/
/*** default is sigma = 1                         ***/
/*** get sigma from RMS preevent noise ? level    ***/
/*** error is 1.96 * sqrt( diag(CV_matrix) )      ***/
/****************************************************/

void compute_error( int rows, int cols, float **v_matrix, float *w_vector, float *x_vector, float *s_vector, float *b_vector,
        float *logerrMdc, float *logerrMclvd, float *logerrMiso )
{
/* see include/mt.h Mo = math.pow( 10.0, 1.5*(Mw+10.73) ) = 1.2445146117713818e+16; Reference Mw = 0.0 */

	int i, j;
	float data_mean, data_variance, residual_variance;
	float tmp, tmp1, tmp2, tmp3, tmpa, tmpb, tmpc, tmpx, tmpy, tmpz;

	float *e_vector;    /*** error vector             e[1..cols]          ***/
	float **cv_matrix;  /*** covariance matrix       cv[1..cols][1..cols] ***/

	float Mdc, Mclvd, Miso;
	float Mdc_err, Mclvd_err, Miso_err;
	int verbose = 1;

	float **matrix( int, int, int, int );
	float *vector( int, int );
	float mean( float *, int );
	float variance( float *, int, float );
	float root_mean_square_variance( float *, float *, int );
	void diag( int, float **, float * );
	void matmul( int, float **, int, float *, int, float * );
	void svdvar( float **, int, float *, float ** );

	if(verbose)fprintf( stdout, "%s: allocating memory for cv_matrix \n", progname );
	cv_matrix = matrix( 1, cols+1, 1, cols+1 );
	if(verbose)fprintf( stdout, "%s: allocating memory for e_vector\n", progname );
	e_vector = vector( 1, cols+1 );

	svdvar( v_matrix, cols, w_vector, cv_matrix );

	residual_variance = root_mean_square_variance( b_vector, s_vector, rows );
	data_mean = mean( b_vector, rows );
	data_variance = variance( b_vector, rows, data_mean );
	
	for( j=1; j<=cols; j++ )
	{
		for( i=1; i<=cols; i++ )
		{
			tmp = w_vector[j] * cv_matrix[i][j];
			cv_matrix[i][j] = tmp;
		}
	}

	diag( cols, cv_matrix, e_vector );

	printf( "\nresidual_variance=%g datamean=%g datavariance=%g\n",
		residual_variance, data_mean, data_variance );

	printf( "            Mdc           Mclvd            Miso\n" );

	for( j=1; j<=cols; j++ )
	{
		for( i=1; i<=cols; i++ )
		{
			printf("%15.5e ", sqrt(fabs(cv_matrix[i][j])) * base_moment );
		}
		printf("\n");
	}

	printf( "Diag covm scaled by data variance=\n" );
	for( j=1; j<=cols; j++ )
		printf("%15.5e ", data_variance * e_vector[j] * base_moment );
	printf("\n");

	printf("Diag covm scaled by residual_variance=\n" );
	residual_variance = 1.96;
	for( j=1; j<=cols; j++ )
		printf("%15.5e ", residual_variance * sqrt(e_vector[j]) * base_moment );
	printf("\n");

	Mdc_err   = residual_variance * sqrt(e_vector[1]) * base_moment;
	Mclvd_err = residual_variance * sqrt(e_vector[2]) * base_moment;
	Miso_err  = residual_variance * sqrt(e_vector[3]) * base_moment;
/*
	Mdc_err   = data_variance * e_vector[1] * base_moment;
	Mclvd_err = data_variance * e_vector[2] * base_moment;
	Miso_err  = data_variance * e_vector[3] * base_moment;
*/
	Mdc   = fabs(x_vector[1]) * base_moment;
	Mclvd = fabs(x_vector[2]) * base_moment;
	Miso  = fabs(x_vector[3]) * base_moment;

	printf("Solution vector:\n");
	printf("%15.5e %15.5e %15.5e\n", Mdc, Mclvd, Miso );

	tmpa = log10(Mdc);
	tmpb = log10(Mclvd);
	tmpc = log10(Miso);

	tmp1 = log10(Mdc_err)/log10(Mdc);
	tmp2 = log10(Mclvd_err)/log10(Mclvd);
	tmp3 = log10(Miso_err)/log10(Miso);

	*logerrMdc   = log10(Mdc_err);
	*logerrMclvd = log10(Mclvd_err);
	*logerrMiso  = log10(Miso_err);

	printf("%15.5e %15.5e %15.5e\n", tmp1, tmp2, tmp3 );

	tmpx = pow( 10, tmpa - tmp1 );
	tmpy = pow( 10, tmpb - tmp2 );
	tmpz = pow( 10, tmpc - tmp3 );
	printf("%15.5e %15.5e %15.5e\n", tmpx, tmpy, tmpz );

	tmpx = pow( 10, tmpa + tmp1 );
	tmpy = pow( 10, tmpb + tmp2 );
	tmpz = pow( 10, tmpc + tmp3 );
	printf("%15.5e %15.5e %15.5e\n", tmpx, tmpy, tmpz );
	printf("\n");

	if(verbose) fprintf( stdout, "%s: freeing memory cv_matrix\n", progname );
	free_matrix( cv_matrix, 1, cols+1, 1, cols+1 );

	if(verbose) fprintf( stdout, "%s: freeing memory e_vector\n", progname );
	free_vector( e_vector, 1, cols+1 );

}
