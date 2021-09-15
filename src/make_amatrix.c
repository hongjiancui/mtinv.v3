#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../include/nrutil.h"
#include "../include/mt.h"

static float half=0.5;
static float sixth=0.166666667;
static float third=0.333333333;


/*********************************************************************************/
/*** this version of make_amatrix2 uses the formaution of Jost&Herrmann (1989) ***/
/*** only deviatoric moment tensors not full mt                                ***/
/*********************************************************************************/

void make_amatrix2( Greens **grn, EventInfo *ev, int nsta, int iz, float **a_matrix, 
	float *b_vector, int mtdegfree, int Distance_Normalize, float DistNormR0, FixISOZ myfixisoz )
{
	int it, nt, ista, irow;
	float fi, d2r;
	float *rss, *rds, *rdd, *rep, *zss, *zds, *zdd, *zep, *tss, *tds;
	float Rnorm = 1;
	float wt = 1;
	int verbose = 0;
	int ziso;

/************************************************/
/*** set some constants                       ***/
/************************************************/

	d2r = M_PI/180;
	if(verbose)fprintf(stderr, "make_amatrix2(iz=%d)\n", iz );

/************************************************/
/*** station loop                             ***/
/************************************************/

	irow = 1;
	for( ista = 0; ista < nsta; ista++ )
	{

	/*****************************************************/
	/*** skip comment out stations in the calculations ***/
	/*****************************************************/
		if( ev[ista].iused < 1 )
		{
			if(verbose)
			{
			  fprintf(stderr, 
			    "make_amatrix2(iz=%d): skiping %d %s %s\n",
				iz, ista, ev[ista].stnm, ev[ista].net );
			}
			continue;
		}

	/************************************************/
	/*** local variables                          ***/
	/************************************************/

		nt = grn[ista][iz].nt;
		fi = grn[ista][iz].az;
		wt = ev[ista].weight;

		if( Distance_Normalize )
		{
			Rnorm = grn[ista][iz].rdist/DistNormR0;
			if(verbose)
			{
			  fprintf(stdout, 
			    "make_amatrix2: sta=%s R0=%g dist=%g Rnorm=%g\n",
				ev[ista].stnm, DistNormR0, grn[ista][iz].rdist, Rnorm );
			}
		}

		if( myfixisoz.iswitch )
			ziso = myfixisoz.indexz;
		else
			ziso = iz;

	/************************************************/
	/*** create pointers to the Green's functions ***/
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

	/************************************************/
	/*** transverse/tangential component          ***/
	/************************************************/

		for( it = 0; it < nt; it++ )
		{
		  a_matrix[irow][1] = (  half * sin(2.0*fi) * tss[it] ) * Rnorm * wt;
		  a_matrix[irow][2] = ( -half * sin(2.0*fi) * tss[it] ) * Rnorm * wt;
		  a_matrix[irow][3] = (        -cos(2.0*fi) * tss[it] ) * Rnorm * wt;
		  a_matrix[irow][4] = (             sin(fi) * tds[it] ) * Rnorm * wt;
		  a_matrix[irow][5] = (            -cos(fi) * tds[it] ) * Rnorm * wt;
		  a_matrix[irow][6] = 0;
		  b_vector[irow]    = ev[ista].ew.data[it] * Rnorm * wt;
		  irow++;
		}

	/************************************************/
	/*** radial components ***/
	/************************************************/

		for( it = 0; it < nt; it++ )
		{
		  a_matrix[irow][1] = ( +half*cos(2*fi) * rss[it] -half*rdd[it] ) * Rnorm * wt;
		  a_matrix[irow][2] = ( -half*cos(2*fi) * rss[it] -half*rdd[it] ) * Rnorm * wt;
		  a_matrix[irow][3] = (       sin(2*fi) * rss[it] ) * Rnorm * wt;
		  a_matrix[irow][4] = (         cos(fi) * rds[it] ) * Rnorm * wt;
		  a_matrix[irow][5] = (         sin(fi) * rds[it] ) * Rnorm * wt;
		  a_matrix[irow][6] = 0;
		  b_vector[irow]    = ev[ista].ns.data[it] * Rnorm * wt;
		  irow++;
		}

	/************************************************/
	/*** vertical component ***/
	/************************************************/
		for( it = 0; it < nt; it++ )
		{
		  a_matrix[irow][1] = ( +half*cos(2*fi) * zss[it] -half*zdd[it] ) * Rnorm * wt;
		  a_matrix[irow][2] = ( -half*cos(2*fi) * zss[it] -half*zdd[it] ) * Rnorm * wt;
		  a_matrix[irow][3] = (       sin(2*fi) * zss[it] ) * Rnorm * wt;
		  a_matrix[irow][4] = (         cos(fi) * zds[it] ) * Rnorm * wt;
		  a_matrix[irow][5] = (         sin(fi) * zds[it] ) * Rnorm * wt;
		  a_matrix[irow][6] = 0;
		  b_vector[irow]    = ev[ista].z.data[it] * Rnorm * wt;
		  irow++;
		}
	}
}

/***********************************************************************************/
/*** this version of make_amatrix uses the formaution of Herrmann and Hutchenson ***/
/*** (1993) and can be used for both deviatoric and full moment tensors          ***/
/***********************************************************************************/

void make_amatrix( Greens **grn, EventInfo *ev, int nsta, int iz, float **a_matrix, float *b_vector,
	int mtdegfree, int Distance_Normalize, float DistNormR0, FixISOZ myfixisoz )
{
	int it, nt, ista, irow;
	float fi, d2r;
	float *rss, *rds, *rdd, *rep, *zss, *zds, *zdd, *zep, *tss, *tds;
        float Rnorm = 1;
	float wt = 1.0;
        int verbose = 0;
        int ziso;

	d2r = M_PI/180;
	if(verbose)
	{
		if( mtdegfree == FORCE_EXPLOSION )
			printf("inversion = FORCE_EXPLOSION\n");
		if( mtdegfree == DEVIATORIC_MOMENT_TENSOR ) 
			printf("inversion = DEVIATORIC_MOMENT_TENSOR\n");
		if( mtdegfree == FULL_MOMENT_TENSOR )
			printf("inversion = FULL_MOMENT_TENSOR\n");
	}

/************************************************/
/*** station loop                             ***/
/************************************************/
                                                                                                                                    
        irow = 1;
        for( ista = 0; ista < nsta; ista++ )
        {
                                                                                                                                    
        /*****************************************************/
        /*** skip comment out stations in the calculations ***/
        /*****************************************************/
                if( ev[ista].iused < 1 )
                {
                        if(verbose)
                        {
                          fprintf(stderr,
                            "make_amatrix(iz=%d): skiping %d %s %s\n",
                                iz, ista, ev[ista].stnm, ev[ista].net );
                        }
                        continue;
                }
                                                                                                                                    
        /************************************************/
        /*** local variables                          ***/
        /************************************************/

                nt = grn[ista][iz].nt;
                fi = grn[ista][iz].az * d2r;
		wt = ev[ista].weight;

                if( Distance_Normalize )
                {
                        Rnorm = grn[ista][iz].rdist / DistNormR0;

                        if(verbose)
                        {
                          fprintf(stdout,
                            "make_amatrix: sta=%s R0=%g dist=%g Rnorm=%g\n",
                                ev[ista].stnm, DistNormR0, grn[ista][iz].rdist, Rnorm );
                        }
                }
                                                                                                                                    
                if( myfixisoz.iswitch )
                        ziso = myfixisoz.indexz;
                else
                        ziso = iz;
                                                                                                                                    
        /************************************************/
        /*** create pointers to the Green's functions ***/
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
                                                                                                                                    
        /************************************************/
        /*** transverse/tangential component          ***/
        /************************************************/
                                                                                                                                    
                for( it = 0; it < nt; it++ )
                {
		  if( mtdegfree == FORCE_EXPLOSION )
		  {
			a_matrix[irow][1] = 0;
			a_matrix[irow][2] = 0;
			a_matrix[irow][3] = 0;
			a_matrix[irow][4] = 0;
			a_matrix[irow][5] = 0;
			a_matrix[irow][6] = 0;
		  }
		  else
		  {
                     a_matrix[irow][1] = (  half * sin(2*fi) * tss[it] ) * Rnorm * wt;
                     a_matrix[irow][2] = ( -half * sin(2*fi) * tss[it] ) * Rnorm * wt;
                     a_matrix[irow][3] = (        -cos(2*fi) * tss[it] ) * Rnorm * wt;
                     a_matrix[irow][4] = (           sin(fi) * tds[it] ) * Rnorm * wt;
                     a_matrix[irow][5] = (          -cos(fi) * tds[it] ) * Rnorm * wt;
                     a_matrix[irow][6] = 0;
		  }
                  b_vector[irow]    = ev[ista].ew.data[it] * Rnorm * wt;
                  irow++;
		}

        /************************************************/
        /*** radial components ***/
        /************************************************/
                                                                                                                                    
                for( it = 0; it < nt; it++ )
                {
		  if( mtdegfree == FORCE_EXPLOSION )
		  {
			a_matrix[irow][1] = ( third * rep[it] ) * Rnorm * wt;
			a_matrix[irow][2] = ( third * rep[it] ) * Rnorm * wt;
			a_matrix[irow][3] = 0;
			a_matrix[irow][4] = 0;
			a_matrix[irow][5] = 0;
			a_matrix[irow][6] = ( third * rep[it] ) * Rnorm * wt;
		  }
		  else if( mtdegfree == DEVIATORIC_MOMENT_TENSOR )
		  {
                    a_matrix[irow][1] = ( +half*cos(2*fi) * rss[it] -half*rdd[it] ) * Rnorm * wt;
                    a_matrix[irow][2] = ( -half*cos(2*fi) * rss[it] -half*rdd[it] ) * Rnorm * wt;
                    a_matrix[irow][3] = (       sin(2*fi) * rss[it] ) * Rnorm * wt;
                    a_matrix[irow][4] = (         cos(fi) * rds[it] ) * Rnorm * wt;
                    a_matrix[irow][5] = (         sin(fi) * rds[it] ) * Rnorm * wt;
                    a_matrix[irow][6] = 0;
		  }
		  else if( mtdegfree == FULL_MOMENT_TENSOR )
		  {
		    a_matrix[irow][1] = ( +half*cos(2*fi) * rss[it] -sixth*rdd[it] + third*rep[it] ) * Rnorm * wt;
                    a_matrix[irow][2] = ( -half*cos(2*fi) * rss[it] -sixth*rdd[it] + third*rep[it] ) * Rnorm * wt;
                    a_matrix[irow][3] = (       sin(2*fi) * rss[it] ) * Rnorm * wt;
                    a_matrix[irow][4] = (         cos(fi) * rds[it] ) * Rnorm * wt;
                    a_matrix[irow][5] = (         sin(fi) * rds[it] ) * Rnorm * wt;
                    a_matrix[irow][6] = ( third * rdd[it] + third * rep[it] ) * Rnorm * wt;
		  }
                  b_vector[irow]    = ev[ista].ns.data[it] * Rnorm * wt;
                  irow++;
                }

        /************************************************/
        /*** vertical component ***/
        /************************************************/
                for( it = 0; it < nt; it++ )
                {
		  if( mtdegfree == FORCE_EXPLOSION )
		  {
                        a_matrix[irow][1] = ( third * zep[it] ) * Rnorm * wt;
                        a_matrix[irow][2] = ( third * zep[it] ) * Rnorm * wt;
                        a_matrix[irow][3] = 0;
                        a_matrix[irow][4] = 0;
                        a_matrix[irow][5] = 0;
                        a_matrix[irow][6] = ( third * zep[it] ) * Rnorm * wt;
                  }
                  else if( mtdegfree == DEVIATORIC_MOMENT_TENSOR )
                  {
                    a_matrix[irow][1] = ( +half*cos(2*fi) * zss[it] -half*zdd[it] ) * Rnorm * wt;
                    a_matrix[irow][2] = ( -half*cos(2*fi) * zss[it] -half*zdd[it] ) * Rnorm * wt;
                    a_matrix[irow][3] = (       sin(2*fi) * zss[it] ) * Rnorm * wt;
                    a_matrix[irow][4] = (         cos(fi) * zds[it] ) * Rnorm * wt;
                    a_matrix[irow][5] = (         sin(fi) * zds[it] ) * Rnorm * wt;
                    a_matrix[irow][6] = 0;
		  }
		  else if( mtdegfree == FULL_MOMENT_TENSOR )
                  {
                    a_matrix[irow][1] = ( +half*cos(2*fi) * zss[it] -sixth*zdd[it] + third*zep[it] ) * Rnorm * wt;
                    a_matrix[irow][2] = ( -half*cos(2*fi) * zss[it] -sixth*zdd[it] + third*zep[it] ) * Rnorm * wt;
                    a_matrix[irow][3] = (       sin(2*fi) * zss[it] ) * Rnorm * wt;
                    a_matrix[irow][4] = (         cos(fi) * zds[it] ) * Rnorm * wt;
                    a_matrix[irow][5] = (         sin(fi) * zds[it] ) * Rnorm * wt;
                    a_matrix[irow][6] = ( third * zdd[it] + third * zep[it] ) * Rnorm * wt;
                  }
                  b_vector[irow]    = ev[ista].z.data[it] * Rnorm * wt;
                  irow++;
                }
	}
}
