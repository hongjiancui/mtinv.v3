#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
 
#include <pthread.h>    /* POSIX Threads */

#include "../include/nrutil.h"
#include "../include/mt.h"
                                                                                                                                                               
char progname[128];

typedef struct {
        int ista;
        int nz;
        float *z;
} DepthVector;

typedef struct thread_data
{
        Greens *grn;
	EventInfo *ev;
	DepthVector *z;
	int ista;
	int nsta;
	int idumpsac;
	int idumpgrn;
	int verbose;
} ThreadData;

void glib2inv_parallel( Greens **grn, EventInfo *ev, DepthVector *z, int nsta, int idumpsac, int idumpgrn, int verbose )
{

/**********************/
/*** pthreads stuff ***/
/**********************/

        pthread_t *thread;
        ThreadData *td;
        pthread_attr_t attr;

	void *glib2inv_staproc_pthread( void *ptr );
	int iz, nz, ista;

/*********************************************************************/
/*** start the subroutine                                          ***/
/*********************************************************************/

	if(verbose)
	{
	  fprintf( stdout,
	    "%s: glib2inv_parallel.c: glib2inv_parallel(): nsta=%d idumpsac=%d idumpgrn=%d verbose=%d\n",
		progname, nsta, idumpsac, idumpgrn, verbose );
	}

	thread = (pthread_t *)malloc( nsta * sizeof(pthread_t) );

	td = malloc( nsta * sizeof(ThreadData) );

	for( ista = 0; ista < nsta; ista++ )
	{
		td[ista].grn          = &grn[ista][0];
		td[ista].ev           = ev;
		td[ista].z            = &z[ista];
		td[ista].ista         = ista;
		td[ista].nsta         = nsta;
		td[ista].idumpsac     = idumpsac;
		td[ista].idumpgrn     = idumpgrn;
		td[ista].verbose      = verbose;
	}

/***
	for( ista = 0; ista < nsta; ista++ )
	{
		for( iz = 0; iz < z[ista].nz; iz++ )
		{
			fprintf( stdout, 
"%s: pthreads: ista=%d iz=%d %s sta=%s net=%s stla=%g stlo=%g stel=%g evla=%g evlo=%g evdp=%g rdist=%g az=%g baz=%g t0=%g dt=%g nt=%d\n", 
				progname,
				ista, 
				iz, 
				td[ista].grn[iz].filename,
				td[ista].grn[iz].stnm,
				td[ista].grn[iz].net,
				td[ista].grn[iz].stla,
				td[ista].grn[iz].stlo,
				td[ista].grn[iz].stel,
				td[ista].grn[iz].evla,
				td[ista].grn[iz].evlo,
				td[ista].grn[iz].evdp,
				td[ista].grn[iz].rdist,
				td[ista].grn[iz].az,
				td[ista].grn[iz].baz,
				td[ista].grn[iz].t0,
				td[ista].grn[iz].dt,
				td[ista].grn[iz].nt );
		}
	}
***/

/*********************************************************************/
/*** start the pthreads                                            ***/
/*********************************************************************/

	pthread_attr_init( &attr );

	pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_JOINABLE );

	/* pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_DETACHED ); */

	for( ista = 0; ista < nsta; ista++ )
	{
		fprintf( stderr, "%s: glib2inv_parallel(): STDERR: calling pthread_create(): ista=%d\n",
			progname, ista );

		pthread_create( &thread[ista], &attr, glib2inv_staproc_pthread, &td[ista] );
	}

	for( ista = 0; ista < nsta; ista++ )
	{
		pthread_join( thread[ista], NULL );
	}

	pthread_attr_destroy( &attr );

} /**** end of glib2inv_parallel.c ****/



/*******************************************************************************/
/*******************************************************************************/
/***** void *glib2inv_staproc_pthread( void *ptr )                          ****/
/*******************************************************************************/
/*******************************************************************************/

void *glib2inv_staproc_pthread( void *ptr )
{
	int DIFFoperator = 3;
	int mtdegfree = 5;

/********************/
/*** Greens stuff ***/
/********************/

	int ig, ng=10, MAX_ARRAY_SIZE=4096;
	float **garray;

	FILE *fpout;

	int old_nt;
	char taper_type[3];
	float taper_frac_width;

/***************************/
/*** function prototypes ***/
/***************************/

        /*** tdif/Differentiates.c ***/
        void differentiate( float *x, int n, float dt, int op, int verbose );
                                                                                                                                                               
        /*** source/source_subs.c ***/
        void source_time_function( float *data, int nt, float dt, float tr, float tt );
                                                                                                                                                               
        /*** filter/filtersubs.c ***/
        void iir_filter( float *data,
                        int nsamps,
                        char *filter_type,
                        float trbndw,
                        float a,
                        int iord,
                        char *operation_type,
                        float flo,
                        float fhi,
                        float ts,
                        int passes );

        /*** Ichinose Feb2010 ***/
        /*** substitute interpolate_fft with interpolate_wiggins ***/
        /*** interpolate/interpolate_subs.c ***/
        void interpolate_fft(   float *data,
                                int old_npts, float old_delta,
                                int *new_npts, float new_delta );
                                                                                                                                                               
        /*** wiggins/wiggins_sub.c ***/
        void interpolate_wiggins2( float *data, int npts, float delta,
                                float b, int new_nt, float new_dt, int verbose );
                                                                                                                                                               
        /*** glib2inv_subs.c ***/
        void grn2disp( Greens *g, EventInfo *ev, int verbose, int mtdegfree );
                                                                                                                                                               
        /*** wrtgrn2sac.c ***/
        void wrtgrn2sac( Greens *g, int ista );
                                                                                                                                                               
        /*** glib2inv_subs.c ***/
        void split2grn( Greens *g, float **garray );
                                                                                                                                                               
        /*** glib2inv_subs.c ***/
        void array2grn( float **garray, Greens *g );
                                                                                                                                                               
        /*** envelope/envelope_sub.c ***/
        void envelope( float *y, int npts, float dt );
 
/******************************/
/*** thread Data interface ****/
/******************************/

	ThreadData *td;
	td = (ThreadData *) ptr;

	Greens *grn;
	EventInfo *ev;
	DepthVector *z;
	int nsta;
	int verbose;
	int iz;
	int ista;
	int idumpsac;
	int idumpgrn;

/***************************************************************************************************/
/**** begin program                                                                              ***/
/***************************************************************************************************/

	ista  = td->ista;
	nsta  = td->nsta;
	grn   = (Greens *)    td->grn;
	ev    = (EventInfo *) td->ev;
	z     = (DepthVector *) td->z;
	idumpsac = td->idumpsac;
	idumpgrn = td->idumpgrn;
	
	garray = (float **)malloc( ng * sizeof(float *) );
	for( ig=0; ig<ng; ig++ )
	{
		garray[ig] = (float *)calloc( MAX_ARRAY_SIZE, sizeof(float) );
	}

	for( iz = 0; iz < z->nz; iz++ )
	{

        /******************************************************************/
        /*** demultiplex from structure to array of 10 greens functions ***/
        /*** and loop over the 10 fundamental faulting orientations     ***/
        /******************************************************************/
                                                                                                                                                               
		split2grn( &grn[iz], garray );
                                                                                                                                                               
		for( ig = 0 ; ig < ng; ig++ )
		{
			if( ev[ista].grd_mo_type == VELOCITY )
			{
				DIFFoperator = 3;
				differentiate(
					garray[ig],
					grn[iz].nt,
					grn[iz].dt,
					DIFFoperator,
					verbose );
			}


		/***********************************/
		/*** convolve trapazoid function ***/
		/***********************************/
                                                                                                                                                               
			source_time_function(
				garray[ig],
				grn[iz].nt,
				grn[iz].dt,
				ev[ista].tr,
				ev[ista].tt );
                                                                                                                                                       
		/****************************************/
		/*** bandpass filter greens functions ***/
		/****************************************/
                                                                                                                                                               
			iir_filter(
				garray[ig],
				grn[iz].nt,
				"BU",
				ev[ista].trbndw,
  				ev[ista].a,
  				ev[ista].npole,
  				"BP",
  				ev[ista].lf,
  				ev[ista].hf,
  				grn[iz].dt,
  				ev[ista].npass );
                                                                                                                                                               
		/*****************************************************/
		/*** interpolate greens functions to new nt and dt ***/
		/*****************************************************/
                                                                                                                                                               
			interpolate_fft(
				garray[ig],
				grn[iz].nt,
				grn[iz].dt,
				&old_nt,
				ev[ista].dt );
                                                                                                                                                               
		/*******************************************/
		/*** taper with the new parameters       ***/
		/*** do not need to taper the synthetics ***/
		/*******************************************/
		/***
			strcpy( taper_type, "h\0" );
			taper_frac_width = 0.40;
			if( ev[ista].dt < 0.85 ) taper_frac_width = 0.25;
			if( ev[ista].dt < 0.50 ) taper_frac_width = 0.10;
			taper_frac_width = 0.05;
			taper( garray[ig], ev[ista].nt,
			taper_type, taper_frac_width );
		***/
                                                                                                                                                       
		/*******************************************************************/
		/*** compute envelope function of the greens function synthetics ***/
		/*******************************************************************/
                                                                                                                                                               
			if( ev[ista].ienvelope == 1 )
			{
  				if( verbose )
  				{
					fprintf( stdout, "%s: ista = %d computing envelope \n",
						progname, ista );
  				}
  				envelope( garray[ig], ev[ista].nt, ev[ista].dt );
			}

		} /*** loop over Green function - ig ***/
                                                                                                                                                              
	/***********************************************************/
	/*** reset nt and delta from interpolation or decimation ***/
	/***********************************************************/
                                                                                                                                                               
		grn[iz].nt = ev[ista].nt;
		grn[iz].dt = ev[ista].dt;
                                                                                                                                                               
	/********************************************/
	/*** convert back from array to structure ***/
	/********************************************/
                                                                                                                                                               
		array2grn( garray, &grn[iz] );
                                                                                                                                                               
	} /*** loop over depth - iz ***/

	/* free(garray); */

/**************************************************/
/*** write out greens functions as seperate sac ***/
/***    files for inspection                    ***/
/**************************************************/

	if( idumpgrn )
	{
		for( iz = 0; iz < z->nz; iz++ )
		{
			wrtgrn2sac( &grn[iz], ista );
		}
	}

/****************************************************/
/*** if event tag present in input PAR file then  ***/
/*** compute displacement synthetics              ***/
/****************************************************/

	if( idumpsac )
	{
	 	for( iz = 0; iz < z->nz; iz++ )
		{
			if(     ( ev[ista].my_z == z->z[iz] ) &&
				( ev[ista].str  != -999  ) &&
				( ev[ista].dip  != -999  ) &&
				( ev[ista].rak  != -999  ) &&
				( ev[ista].Mw   != -999  )  )
			{
				grn[iz].ver = calloc( grn[iz].nt, sizeof(float) );
                        	grn[iz].rad = calloc( grn[iz].nt, sizeof(float) );
                        	grn[iz].tra = calloc( grn[iz].nt, sizeof(float) );

                        	grn2disp( &(grn[iz]), &ev[ista], verbose, mtdegfree );

                        /*** write out the synthetics here ***/

			}
		}
	}

/****************************************************/
/*** write out the ginv file ***/
/****************************************************/

	if( (fpout = fopen( ev[ista].ginv_filename, "wb" )) == NULL )
	{
		fprintf( stderr,
			"%s: glib2inv.c: glib2inv(): STDERR: Fatal Error, cannot open file %s\n",
				progname, ev[ista].ginv_filename );
		exit(-1);
	}

	fprintf( stderr,
	  "%s: glib2inv_parallel.c: glib2inv_staproc_pthreads(): STDERR: ista=%d %s.%s nz=%d nt=%d dt=%g writing file %s\n",
		progname,
		ista,
		grn[0].stnm,
		grn[0].net,
		z->nz,
		grn[0].nt,
		grn[0].dt,
		ev[ista].ginv_filename );

	fwrite( &(z->nz), sizeof(int), 1, fpout );
	fwrite( &(z->z[0]), z->nz * sizeof(float), 1, fpout );

	for( iz = 0; iz < z->nz; iz++ )
	{
		fwrite( &grn[iz], sizeof(Greens), 1, fpout );
	}

	fclose(fpout);

	pthread_exit((void *)0);

} /*** end glib2inv_staproc_pthreads() ***/
