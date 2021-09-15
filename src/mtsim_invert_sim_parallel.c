#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>
#include <pthread.h>

#include "../include/nrutil.h"     /** numerical recipes **/
#include "../include/mt.h"         /** global datatype and structure declarations **/

char progname[128];

typedef struct thread_data
{
	int k;
	EventInfo *ev;
	Greens **grn;
	Solution *psol;
	float *res;
	int rows;
	int cols;
	float *best_b_vector;
	int Nmodels;
	int verbose;

	int nsta;
	int Distance_Normalize;
        float DistNormR0;
        int iz_best;
	
	int indexz;
	int iswitch;
	float z;

} ThreadData;

/***********************************************/
/*** INVERT_SIM_PARALLEL()                   ***/
/***********************************************/

#define NTHREADS 8

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
		int verbose )
{

/**************************************************/
/*** pthreads stuff ***/
/**************************************************/

	pthread_t thread[NTHREADS];
	ThreadData *td;
	pthread_attr_t attr;
	int ith, NUM_SIM_PER_THREAD;

/**************************************************/
/*** variable declarations ***/
/**************************************************/
	int mtdegfree = 6;
	FILE *fp;
	Solution *psol;
	int i, j, k;
	int rows, cols;

/****************************/
/**** function prototypes ***/
/****************************/

/*** mtsim_invert_sim_parallel.c ***/

	void *invert_sim_thread( void *ptr );

/*** math/math.c ***/
	float **matrix(  int nrl, int nrh, int ncl, int nch );
	void free_matrix(float **m, int nrl, int nrh, int ncl, int nch );

/*** math/math.c ***/
	float *vector(int nl, int nh );
	void free_vector(float *v, int nl, int nh );

/***  mtsim_subs.c ***/
	int size_A_matrix( EventInfo *ev, Greens **grn, int nsta, int iz );

/*** mtsim_subs.c ***/
	void inversion_init_mem( int rows, int cols,
                        float **a_matrix,
                        float **u_matrix,
                        float *w_vector,
                        float *e_vector,
                        float *x_vector,
                        float **cv_matrix,
                        float **v_matrix,
                        float *b_vector,
                        float *s_vector );

/*** make_amatrix.c ***/
	void make_amatrix( Greens **grn,
                           EventInfo *ev,
                           int nsta,
                           int iz_best,
                           float **a_matrix,
                           float *b_vector,
                           int mtdegfree,
                           int Distance_Normalize,
                           float DistNormR0,
                           FixISOZ myfixisoz );

/*** mtsim_subs.c ***/
	void matrix_copy( int rows, int cols, float **out, float **in );

/*** mtsim_subs.c ***/
	void vector_copy( int rows, float *out, float *in );

/*** mtsim_subs.c ***/
	void singular_value_decomposition( 
		int rows,
		int cols,
		float **u_matrix,
                float *w_vector,
		float **v_matrix,
		float *b_vector,
		float *x_vector );

/*** math/math.c ***/
	void matmul( int, float **, int, float *, int, float * );
	
/*** mtsim_subs.c ***/
	void fitness( Solution *sol, int iz, float *b_vector, float *s_vector, int rows, int mtdegfree, int verbose );

/*** mtsim_subs.c ***/
	void set_mt_sol( EventInfo *ev, Solution *sol, float *x_vector, int iz, int mtdegfree, int verbose );

/*** mtsim_subs.c ***/
	void write_mt_sim( FILE *fp, int nsol, Solution *sol, int verbose );

/*** mtsim_subs.c ***/
	void BS_replace( float *res, int rrows, float *best_b_vector, int brows, float *b_vector, int iseed, int verbose );

	fprintf( stdout, "%s: FixISOZ: iswitch=%d indexz=%d z=%g\n", 
			progname, myfixisoz.iswitch, myfixisoz.indexz, myfixisoz.z );
	fprintf( stderr, "%s: FixISOZ: iswitch=%d indexz=%d z=%g\n",
                        progname, myfixisoz.iswitch, myfixisoz.indexz, myfixisoz.z );

/*************************************************/
/*** get the total data length rows the matrix ***/
/*************************************************/
	
	cols = 6;
	rows = size_A_matrix( ev, grn, nsta, iz_best );

	if(verbose)
	{
	  fprintf( stdout, "%s: invert_sim_parallel(): A matrix: rows=%d cols=%d\n",
                progname, rows, cols);
	  fflush(stdout);
	}

/*************************/
/*** Allocating memory ***/
/*************************/

        if(verbose)
        {
          fprintf(stdout, "%s: allocating memory for iz=%d\n", progname, iz_best );
          fprintf(stdout, "%s: Allocating memory for a_matrix rows=%d cols=%d\n",
                progname, rows, cols );
	  fflush(stdout);
        }

	psol = malloc( Nmodels * sizeof(Solution) );

	td = malloc( Nmodels * sizeof(ThreadData) );

	NUM_SIM_PER_THREAD = rint( Nmodels/NTHREADS );

	fprintf( stdout, "%s: NTHREADS=%d NUM_SIM_PER_THREAD=%d\n", progname, NTHREADS, NUM_SIM_PER_THREAD );
	fprintf( stderr, "%s: NTHREADS=%d NUM_SIM_PER_THREAD=%d\n", progname, NTHREADS, NUM_SIM_PER_THREAD );

	for( k = 0; k < Nmodels; k++ )
	{
		td[k].k        = k;
		td[k].ev       = ev;
		td[k].psol     = psol;
		td[k].grn      = grn;
		td[k].res      = res;
		td[k].rows     = rows;
		td[k].cols     = cols;
		td[k].best_b_vector = best_b_vector;
		td[k].Nmodels  = Nmodels;
		td[k].verbose  = verbose;
	
		td[k].nsta = nsta;
		td[k].Distance_Normalize =  Distance_Normalize;
		td[k].iz_best     = iz_best;
		td[k].DistNormR0  = DistNormR0;

		td[k].iswitch   = (int)myfixisoz.iswitch;
		td[k].indexz    = (int)myfixisoz.indexz;
		td[k].z         = (float)myfixisoz.z;
	
	}

/*********************************************************************/
/*** start the pthreads                                            ***/
/*********************************************************************/

	pthread_attr_init( &attr );

	pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_JOINABLE );

	for( i = 0; i < NUM_SIM_PER_THREAD; i++ )
	{
		for( ith = 0; ith < NTHREADS; ith++ )
		{
			k = ith + i * NTHREADS;

			if( k > Nmodels-1 )
                	{
				break;
			}
			else
			{
				if(verbose)
				{
					fprintf( stderr, "%s: launching k=%05d i=%05d ith=%05d\n", progname, k, i, ith );
					fflush( stderr );
				}

				pthread_create( &thread[ith], &attr, invert_sim_thread, &td[k] );
			}

		} /*** ith - thread id pcreate ***/

		for( ith = 0; ith < NTHREADS; ith++ )
		{
			k = ith + i * NTHREADS;

			if( k > Nmodels-1 )
			{
				break;
			}
			else
			{
				if(verbose)
				{
					fprintf( stderr, "%s: joining thread k=%05d i=%05d ith=%05d\n",
						progname, k, i, ith );
					fflush(stderr);
				}
			}
			pthread_join( thread[ith], NULL );
		}
	}

	pthread_attr_destroy( &attr );

/*** write out results ***/

	if( (fp = fopen( "mtsim.out", "w" )) == NULL )
        {
                fprintf( stdout, "%s: cannot open file mtsim.out\n", progname );
                exit(-1);
        }
	write_mt_sim( fp, Nmodels, psol, verbose );
	fclose(fp);
	free(psol);

/****************************/
/*** free memory clean up ***/
/****************************/

	if(verbose)
        {
          fprintf( stdout, "%s: mtsim_invert_sim_parallel.c: freeing memory inside invert\n\trows=%d\n\tcols=%d\n\n",
                progname, rows, cols );
	  fflush(stdout);
	}

/***
	if(verbose) { fprintf(stdout, "%s: freemem a_matrix\n", progname); fflush(stdout); }
	free_matrix( a_matrix, 0, rows+1, 0, cols+1 ); 

	if(verbose) { fprintf(stdout, "%s: freemem u_matrix\n", progname); fflush(stdout); }
	free_matrix( u_matrix, 0, rows+1, 0, cols+1 );

	if(verbose) { fprintf(stdout, "%s: freemem v_matrix\n", progname); fflush(stdout); }
	free_matrix( v_matrix, 0, cols+1, 0, cols+1 );

	if(verbose) { fprintf(stdout, "%s: freemem cv_matrix\n",progname); fflush(stdout); }
	free_matrix( cv_matrix, 0, cols+1, 0, cols+1 );

	if(verbose) { fprintf(stdout, "%s: freemem x_vector\n", progname); fflush(stdout); }
	free_vector( x_vector, 0, cols+1 );

	if(verbose) { fprintf(stdout, "%s: freemem s_vector\n", progname); fflush(stdout); }
	free_vector( s_vector, 0, rows+1 );

	if(verbose) { fprintf(stdout, "%s: freemem w_vector\n", progname); fflush(stdout); }
	free_vector( w_vector, 0, cols+1 );

	if(verbose) { fprintf(stdout, "%s: freemem b_vector\n", progname); fflush(stdout); }
	free_vector( b_vector, 0, rows+1 );

	if(verbose) { fprintf(stdout, "%s: freemem e_vector\n", progname); fflush(stdout); }
	free_vector( e_vector, 0, cols+1 );
**/
	if(verbose)
	{
		fprintf( stdout, "%s: Exiting from invert()\n", progname );
		fflush(stdout);
	}

} /*** END OF INVERT_SIM_PARALLEL() ***/
	

/*************************************************************************/
/*** INVERT_SIM_THREAD()                                               ***/
/*************************************************************************/

void *invert_sim_thread( void *ptr )
{
	int mtdegfree = 6;

	float **a_matrix;   /*** A matrux with dimensions a[1..rows][1..cols] ***/
        float **u_matrix;   /*** temp space               u[1..rows][1..cols] ***/
        float **v_matrix;   /*** temp space               v[1..cols][1..cols] ***/
        float *w_vector;    /*** temp space               w[1..cols]          ***/
        float *b_vector;    /*** data                     b[1..rows]          ***/
        float *s_vector;    /*** synthetic                s[1..rows]          ***/
        float *x_vector;    /*** solution (mom ten)       x[1..cols]          ***/
        float *e_vector;    /*** error vector             e[1..cols]          ***/
        float **cv_matrix;  /*** covariance matrix       cv[1..cols][1..cols] ***/

/***  mtsim_subs.c ***/
        int size_A_matrix( EventInfo *ev, Greens **grn, int nsta, int iz );
                                                                                                                                                                
/*** mtsim_subs.c ***/
        void inversion_init_mem( int rows, int cols,
                        float **a_matrix,
                        float **u_matrix,
                        float *w_vector,
                        float *e_vector,
                        float *x_vector,
                        float **cv_matrix,
                        float **v_matrix,
                        float *b_vector,
                        float *s_vector );
                                                                                                                                                                
/*** make_amatrix.c ***/
        void make_amatrix( Greens **grn,
                           EventInfo *ev,
                           int nsta,
                           int iz_best,
                           float **a_matrix,
                           float *b_vector,
                           int mtdegfree,
                           int Distance_Normalize,
                           float DistNormR0,
                           FixISOZ myfixisoz );

/*** mtsim_subs.c ***/
        void matrix_copy( int rows, int cols, float **out, float **in );
                                                                                                                                                                
/*** mtsim_subs.c ***/
        void vector_copy( int rows, float *out, float *in );
                                                                                                                                                                
/*** mtsim_subs.c ***/
        void BS_replace( float *res, int rrows, float *best_b_vector, int brows, float *b_vector, int iseed, int verbose );

/*** mtsim_subs.c ***/
        void matrix_copy( int rows, int cols, float **out, float **in );

/*** mtsim_subs.c ***/
        void singular_value_decomposition(
                int rows,
                int cols,
                float **u_matrix,
                float *w_vector,
                float **v_matrix,
                float *b_vector,
                float *x_vector );

/*** math/math.c ***/
        void matmul( int, float **, int, float *, int, float * );
                                                                                                                                                                
/*** mtsim_subs.c ***/
        void fitness( Solution *sol, int iz, float *b_vector, float *s_vector, int rows, int mtdegfree, int verbose );
                                                                                                                                                                
/*** mtsim_subs.c ***/
        void set_mt_sol( EventInfo *ev, Solution *sol, float *x_vector, int iz, int mtdegfree, int verbose );

/*** pthreads stuff ***/

	ThreadData *td;
	td = (ThreadData *) ptr;

	int k;
	EventInfo *ev;
	Greens **grn;
	Solution *psol;
	float *res;
	int rows;
	int cols;
	float *best_b_vector;
	int Nmodels;
	int verbose;
	int Distance_Normalize;
	FixISOZ myfixisoz;
	float DistNormR0;
	int iz_best;
	int nsta;

	k        = (int) td->k;
	ev       = (EventInfo *) td->ev;
	grn      = (Greens **) td->grn;
	psol     = (Solution *) td->psol;
	res      = (float *) td->res;
	rows     = (int) td->rows;
	cols     = (int) td->cols;
	best_b_vector = (float *) td->best_b_vector;
	Nmodels  = (int) td->Nmodels;
        verbose  = (int) td->verbose;
	nsta     = (int) td->nsta;
	iz_best  = (int) td->iz_best;
	Distance_Normalize = (int) td->Distance_Normalize;
	DistNormR0 = (float) DistNormR0;

/* myfixisoz  = (FixISOZ) td->myfixisoz; */
	myfixisoz.indexz = (int)td->indexz;
	myfixisoz.iswitch = (int)td->iswitch;
	myfixisoz.z = (float)td->z;

/********************/
/*** begin loop k ***/
/********************/

	a_matrix  = matrix( 0, rows+1, 0, cols+1 );
        u_matrix  = matrix( 0, rows+1, 0, cols+1 );
        v_matrix  = matrix( 0, cols+1, 0, cols+1 );
        cv_matrix = matrix( 0, cols+1, 0, cols+1 );
        w_vector = vector( 0, cols+1 );
        e_vector = vector( 0, cols+1 );
        x_vector = vector( 0, cols+1 );
        b_vector = vector( 0, rows+1 );
        s_vector = vector( 0, rows+1 );

/*** initalize the memory with zeros ***/

        inversion_init_mem( rows, cols, a_matrix, u_matrix, w_vector, e_vector,
                        x_vector, cv_matrix, v_matrix, b_vector, s_vector );

        make_amatrix( grn, ev, nsta, iz_best, a_matrix, b_vector, mtdegfree,
                        Distance_Normalize, DistNormR0, myfixisoz );

        psol[k].mt_type = FULL_MT;
        psol[k].evlo    = ev[0].evlo;
        psol[k].evla    = ev[0].evla;
        psol[k].evdp    = ev[0].evdp;
        psol[k].ot      = ev[0].ot_shift;

        BS_replace( res, rows, best_b_vector, rows, b_vector, k, verbose );

        matrix_copy( rows, cols, u_matrix, a_matrix );

        singular_value_decomposition(
                        rows,
                        cols,
                        u_matrix,
                        w_vector,
                        v_matrix,
                        b_vector,
                        x_vector );
                                                                                                                                                        
        matmul( 0, a_matrix, cols, x_vector, rows, s_vector );
                                                                                                                                                        
/**** debug by writting out vectors ***/
  /*** sprintf( sacfilename, "dat.%05d.sac", k ); ***/
  /*** wrtnewsac( sacfilename, 1.0, rows, b_vector, 0.0 ); ***/
                                                                                                                                                                
        fitness( psol, k, b_vector, s_vector, rows, mtdegfree, verbose );
                                                                                                                                                        
        set_mt_sol( ev, psol, x_vector, k, mtdegfree, verbose );

/*** free memory ***/

	if(verbose) { fprintf(stdout, "%s: freemem a_matrix\n", progname); fflush(stdout); }
	free_matrix( a_matrix, 0, rows+1, 0, cols+1 );

        if(verbose) { fprintf(stdout, "%s: freemem u_matrix\n", progname); fflush(stdout); }
        free_matrix( u_matrix, 0, rows+1, 0, cols+1 );

        free_matrix( v_matrix, 0, cols+1, 0, cols+1 );
        free_matrix( cv_matrix, 0, cols+1, 0, cols+1 );
        free_vector( x_vector, 0, cols+1 );
        free_vector( s_vector, 0, rows+1 );
        free_vector( w_vector, 0, cols+1 );
        free_vector( b_vector, 0, rows+1 );
        free_vector( e_vector, 0, cols+1 );

	if(verbose) { fprintf(stdout, "%s: freemem done...\n", progname ); fflush(stdout); }

	pthread_exit((void *)0);

} /*** End of INVERT_SIM_THREAD() ***/
