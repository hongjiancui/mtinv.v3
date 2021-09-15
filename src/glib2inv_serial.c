#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
                                                                                                                                                               
#include "../include/nrutil.h"
#include "../include/mt.h"
                                                                                                                                                               
char progname[128];

typedef struct {
        int ista;
        int nz;
        float *z;
} DepthVector;

void glib2inv_serial( Greens **grn, EventInfo *ev, DepthVector *z, int nsta, int idumpsac, int idumpgrn, int verbose )
{

	int iz, ista;
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
 
/***************************************************************************************************/
/**** begin program                                                                              ***/
/***************************************************************************************************/

/******************************************************************************/
/*** loop over stations process Green functions                             ***/
/******************************************************************************/
                                                                                                                                                               
/*******************************************************/
/*** allocate memory for greens function demultiplex ***/
/*******************************************************/

	garray = (float **)malloc( ng * sizeof(float *) );
	for( ig=0; ig<ng; ig++ )
	{
		garray[ig] = (float *)calloc( MAX_ARRAY_SIZE, sizeof(float) );
	}

	for( ista = 0; ista < nsta; ista++ )
	{
		for( iz = 0; iz < z[ista].nz; iz++ )
		{

        /******************************************************************/
        /*** demultiplex from structure to array of 10 greens functions ***/
        /*** and loop over the 10 fundamental faulting orientations     ***/
        /******************************************************************/
                                                                                                                                                               
                        split2grn( &grn[ista][iz], garray );
                                                                                                                                                               
                        for( ig = 0 ; ig < ng; ig++ )
                        {
                                                                                                                                                               
                        /**********************************************************/
                        /*** ground motion type -> differentiate for velocity   ***/
                        /**********************************************************/
                                                                                                                                                               
                                if( ev[ista].grd_mo_type == VELOCITY )
                                {
                                        DIFFoperator = 3;
                                        differentiate(
                                                garray[ig],
                                                grn[ista][iz].nt,
                                                grn[ista][iz].dt,
                                                DIFFoperator,
                                                verbose );
                                }
                                                                                                                                                               
                        /***********************************/
                        /*** convolve trapazoid function ***/
                        /***********************************/
                                                                                                                                                               
                                source_time_function(
                                        garray[ig],
                                        grn[ista][iz].nt,
                                        grn[ista][iz].dt,
                                        ev[ista].tr,
                                        ev[ista].tt );
                                                                                                                                                               
                        /****************************************/
                        /*** bandpass filter greens functions ***/
                        /****************************************/
                                                                                                                                                               
                                iir_filter(
                                        garray[ig],
                                        grn[ista][iz].nt,
                                        "BU",
                                        ev[ista].trbndw,
                                        ev[ista].a,
                                        ev[ista].npole,
                                        "BP",
                                        ev[ista].lf,
                                        ev[ista].hf,
                                        grn[ista][iz].dt,
                                        ev[ista].npass );
                                                                                                                                                               
                        /*****************************************************/
                        /*** interpolate greens functions to new nt and dt ***/
                        /*****************************************************/
                                                                                                                                                               
                                interpolate_fft(
                                        garray[ig],
                                        grn[ista][iz].nt,
                                        grn[ista][iz].dt,
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
                                                                                                                                                               
                        grn[ista][iz].nt = ev[ista].nt;
                        grn[ista][iz].dt = ev[ista].dt;
                                                                                                                                                               
                /********************************************/
                /*** convert back from array to structure ***/
                /********************************************/
                                                                                                                                                               
                        array2grn( garray, &grn[ista][iz] );
                                                                                                                                                               
                } /*** loop over depth - iz ***/

        } /*** loop over stations - ista  ***/
 
	free(garray);

/**************************************************/
/*** write out greens functions as seperate sac ***/
/***    files for inspection                    ***/
/**************************************************/
                                                                                                                                                                                 
        if( idumpgrn )
        {
                for( ista = 0; ista < nsta; ista++ )
                {
                        for( iz = 0; iz < z[ista].nz; iz++ )
                        {
                                wrtgrn2sac( &grn[ista][iz], ista );
                        }
                }
        }
                                                                                                                                                                                 
/****************************************************/
/*** if event tag present in input PAR file then  ***/
/*** compute displacement synthetics              ***/
/****************************************************/
                                                                                                                                                                                 
        if( idumpsac )
        {
                                                                                                                                                                                 
          for( ista = 0; ista < nsta; ista++ )
          {
            for( iz = 0; iz < z[ista].nz; iz++ )
            {
                if(     ( ev[ista].my_z == z[ista].z[iz] ) &&
                        ( ev[ista].str  != -999  ) &&
                        ( ev[ista].dip  != -999  ) &&
                        ( ev[ista].rak  != -999  ) &&
                        ( ev[ista].Mw   != -999  )  )
                {
                        grn[ista][iz].ver = calloc( grn[ista][iz].nt, sizeof(float) );
                        grn[ista][iz].rad = calloc( grn[ista][iz].nt, sizeof(float) );
                        grn[ista][iz].tra = calloc( grn[ista][iz].nt, sizeof(float) );
                                                                                                                                                                                 
                        grn2disp( &(grn[ista][iz]), &ev[ista], verbose, mtdegfree );
                                                                                                                                                                                 
                        /*** write out the synthetics here ***/
                }
                                                                                                                                                                                 
            } /*** iz loop ***/
                                                                                                                                                                                 
          } /*** ista loop ***/
                                                                                                                                                                                 
        } /*** if idumpsac ***/

/******************************************************************************/
/*** loop over stations and write out green functions                       ***/
/******************************************************************************/
                                                                                                                                                               
        for( ista = 0; ista < nsta; ista++ )
        {
                if( (fpout = fopen( ev[ista].ginv_filename, "wb" )) == NULL )
                {
                        fprintf( stderr,
                          "%s: glib2inv.c: glib2inv(): STDERR: Fatal Error, cannot open file %s\n",
                                progname, ev[ista].ginv_filename );
                        exit(-1);
                }
                                                                                                                                                               
                fprintf( stderr,
                  "%s: glib2inv.c: glib2inv(): STDERR: %s.%s nz=%d nt=%d dt=%g writing file %s\n",
                        progname,
                        grn[ista][0].stnm,
                        grn[ista][0].net,
                        z[ista].nz,
                        grn[ista][0].nt,
                        grn[ista][0].dt,
                        ev[ista].ginv_filename );
                                                                                                                                                               
                fwrite( &(z[ista].nz), sizeof(int), 1, fpout );
                fwrite( &(z[ista].z[0]), z[ista].nz * sizeof(float), 1, fpout );
                                                                                                                                                               
                for( iz = 0; iz < z[ista].nz; iz++ )
                {
                        fwrite( &grn[ista][iz], sizeof(Greens), 1, fpout );
                }
                                                                                                                                                               
                fclose(fpout);
        }
                                                                                                                                                               
} /*** end of glib2inv_serial.c ***/
