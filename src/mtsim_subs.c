#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>
#include "../include/nrutil.h"     /** numerical recipes **/
#include "../include/mt.h"         /** global datatype and structure declarations **/

char progname[128];

/*******************************************************/
/*** this subroutine does the inversion by looping   ***/
/*** over all the greens functions in the library    ***/
/*******************************************************/

void invert0( EventInfo *ev, 
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
	FixISOZ myfixisoz ) 
{
	int ista, iz;
	char sac_file_name[128];
	int idumpsac = 0;
	FILE *fp;

/*********************************************************/
/*** cols = 6 colums of symetic moment tensor elements ***/
/*********************************************************/

	int i, j, rows, cols;

	float **a_matrix;   /*** A matrux with dimensions a[1..rows][1..cols] ***/
	float **u_matrix;   /*** temp space               u[1..rows][1..cols] ***/
	float **v_matrix;   /*** temp space               v[1..cols][1..cols] ***/
	float *w_vector;    /*** temp space               w[1..cols]          ***/
	float *b_vector;    /*** data                     b[1..rows]          ***/
	float *s_vector;    /*** synthetic                s[1..rows]          ***/
	float *x_vector;    /*** solution (mom ten)       x[1..cols]          ***/
	float *e_vector;    /*** error vector             e[1..cols]          ***/
	float **cv_matrix;  /*** covariance matrix       cv[1..cols][1..cols] ***/

/*****************************************/
/*** subroutines: function prototypes  ***/
/*****************************************/

	float **matrix( int, int, int, int );

	float *vector( int, int );

	int size_A_matrix( EventInfo *, Greens **, int, int );

	void inversion_init_mem( int, int, float **, float **, float *, float *, float *, 
			float **, float **, float *, float * );

	void make_amatrix(   Greens **grn,
                             EventInfo *ev,
                             int nsta,
                             int iz_best,
                             float **a_matrix,
                             float *b_vector,
                             int mtdegfree,
                             int Distance_Normalize,
                             float DistNormR0,
                             FixISOZ myfixisoz );

	void matrix_copy( int, int, float **, float ** );

	void vector_copy( int, float *, float * );

	void singular_value_decomposition(int rows, int cols, float **u_matrix, 
		float *w_vector, float **v_matrix, float *b_vector, float *x_vector );

	void matmul( int, float **, int, float *, int, float * );

	void compute_residual(float *b_vector,float *s_vector,int rows,float *res,float fac);

	void wrtnewsac( char *FO, float dt, int ns, float *ar, float b );

	void svdvar( float **, int, float *, float ** );
	
	void fitness( Solution *, int, float *, float *, int, int, int );

	void set_mt_sol( EventInfo *, Solution *, float *, int, int, int );

	/*** mtinv_subs.c ***/
        void write_gmt_xy_values( Solution *sol,
                                  EventInfo *ev,
                                  Greens **grn,
                                  int iz,
                                  int nsta,
                                  int verbose );

	void dumpSAC( EventInfo *, Greens **, int, int, int, float *, int );

	void write_mt_sim( FILE *, int, Solution *, int );

/*********************************************************************************************/
/*** LOOP OVER DEPTH *************************************************************************/
/*********************************************************************************************/
	cols = 6;

	for( iz = 0 ; iz < nz ; iz++ )
	{

	/*************************************************/
	/*** get the total data length rows the matrix ***/
	/*************************************************/
		rows = size_A_matrix( ev, grn, nsta, iz );

		if( verbose ) 
		{
		  printf("%s: A matrix: iz=%d rows=%d cols=%d\n", 
			progname, iz, rows, cols);
		}

	/*************************/
	/*** Allocating memory ***/
	/*************************/
		if(verbose)
		{
		 printf("%s: allocating memory for iz=%d\n", progname, iz );
		 fprintf(stdout, "%s: Allocating memory for a_matrix rows=%d cols=%d\n", 
			progname, rows, cols );
		}
		a_matrix  = matrix( 0, rows+1, 0, cols+1 );
		u_matrix  = matrix( 0, rows+1, 0, cols+1 );
		v_matrix  = matrix( 0, cols+1, 0, cols+1 );
		cv_matrix = matrix( 0, cols+1, 0, cols+1 );
		w_vector = vector( 0, cols+1 );
		e_vector = vector( 0, cols+1 );
		x_vector = vector( 0, cols+1 );
		b_vector = vector( 0, rows+1 );
		s_vector = vector( 0, rows+1 );

	/********************************/
	/***** initialize memory ********/
	/********************************/
	
		inversion_init_mem( rows, cols, a_matrix, u_matrix, w_vector, e_vector,
			x_vector, cv_matrix, v_matrix, b_vector, s_vector );	

	/*********************************************************/
	/*** set up the A matrix and data vector for inversion ***/
	/*********************************************************/

		make_amatrix( grn, ev, nsta, iz, a_matrix, b_vector, mtdegfree, 
			Distance_Normalize, DistNormR0, myfixisoz );

	/***************************************/
	/*** copy A matrix into the U matrix ***/
	/***************************************/

		matrix_copy( rows, cols, u_matrix, a_matrix );

	/******************************/
	/***** do the inversion *******/
	/******************************/

		if( verbose )
		{
		  fprintf( stdout, "%s: inverting... calling svdcmp() rows=%d cols=%d\n",
                        progname, rows, cols );
		}

		singular_value_decomposition( rows, cols, u_matrix, w_vector, 
			v_matrix, b_vector, x_vector );

		matmul( 0, a_matrix, cols, x_vector, rows, s_vector );

		fitness( sol, iz, b_vector, s_vector, rows, mtdegfree, verbose );

		set_mt_sol( ev, sol, x_vector, iz, mtdegfree, verbose );

		if( iz == iz_best )
		{
		  *nres = rows;
		  vector_copy( rows, best_b_vector, s_vector );
		  compute_residual( b_vector, s_vector, rows, res, 1.0 );
		  wrtnewsac( "dat.sac", 1.0, rows, b_vector, 0.0);

		  fp = fopen("mtsim.out.orig","w");
		  write_mt_sim( fp, 1, &sol[iz], verbose );
		  fclose(fp);
		}

	/************************************************************/
	/*** calculate the covariance matrix default is sigma = 1 ***/
	/*** get sigma from RMS preevent noise ? level            ***/
	/*** error is 1.96 * sqrt( diag(CV_matrix) )              ***/
	/************************************************************/

		/* svdvar( v_matrix, cols, w_vector, cv_matrix ); */

		if(idumpsac)
		  dumpSAC( ev, grn, nsta, iz, rows, s_vector, verbose );

		sol[iz].evlo = ev[0].evlo;
		sol[iz].evla = ev[0].evla;
		sol[iz].evdp = ev[0].evdp;
		sol[iz].ot   = ev[0].ot_shift;

	/******************************************************************/
	/*** write the output to GMT files for plotting in plotmech.csh ***/
	/******************************************************************/

		write_gmt_xy_values( sol, ev, grn, iz, nsta, verbose );
		
	/****************************/
	/*** free memory clean up ***/
	/****************************/
		if( verbose ) 
		{
		  fprintf( stdout, "%s: mtsim_subs.c: freeing memory inside invert\n\trows=%d\n\tcols=%d\n\n", 
				progname, rows, cols );
		  fflush(stdout);
		}

	/***
		fprintf( stdout, "%s: freeing memory for a_matrix\n", progname); free_matrix( a_matrix, 1, rows+1, 1, cols+1 );
		fprintf( stdout, "%s: freeing memory for u_matrix\n", progname); free_matrix( u_matrix, 1, rows+1, 1, cols+1 );
		fprintf( stdout, "%s: freeing memory for v_matrix\n", progname); free_matrix( v_matrix, 1, cols+1, 1, cols+1 );
		fprintf( stdout, "%s: freeing memory for cv_matrix\n", progname);free_matrix( cv_matrix, 1, cols+1, 1, cols+1 );
		fprintf( stdout, "%s: freeing memory for x_vector\n", progname); free_vector( x_vector, 1, cols+1 );
                fprintf( stdout, "%s: freeing memory for s_vector\n", progname); free_vector( s_vector, 1, rows+1 );
		fprintf( stdout, "%s: freeing memory for w_vector\n", progname); free_vector( w_vector, 1, cols+1 );
		fprintf( stdout, "%s: freeing memory for b_vector\n", progname); free_vector( b_vector, 1, rows+1 );
		fprintf( stdout, "%s: freeing memory for e_vector\n", progname); free_vector( e_vector, 1, cols+1 );
	***/

	}  /*** iz ***/

	if(verbose)
	{
	 fprintf( stdout, "%s: Exiting from invert()\n", progname );
	 fflush(stdout);
	}

} /*** END OF INVERT0() ***/


/***********************************************/
/*** INVERT_SIM()                            ***/
/***********************************************/

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
		int verbose ) 
{

/**************************************************/
/*** variable declarations ***/
/**************************************************/

	int mtdegfree = 6;
	FILE *fp;
        Solution *psol;
        int i, j, k;
	int rows, cols;
        float **a_matrix;   /*** A matrux with dimensions a[1..rows][1..cols] ***/
        float **u_matrix;   /*** temp space               u[1..rows][1..cols] ***/
        float **v_matrix;   /*** temp space               v[1..cols][1..cols] ***/
        float *w_vector;    /*** temp space               w[1..cols]          ***/
        float *b_vector;    /*** data                     b[1..rows]          ***/
        float *s_vector;    /*** synthetic                s[1..rows]          ***/
        float *x_vector;    /*** solution (mom ten)       x[1..cols]          ***/
        float *e_vector;    /*** error vector             e[1..cols]          ***/
        float **cv_matrix;  /*** covariance matrix       cv[1..cols][1..cols] ***/

/****************************/
/**** function prototypes ***/
/****************************/

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
        void compute_residual( float *b_vector, float *s_vector, int rows, float *res_vector, float fac );

/*** mtsim_subs.c ***/
        void fitness( Solution *sol, int iz, float *b_vector, float *s_vector, int rows, int mtdegfree, int verbose );

/*** mtsim_subs.c ***/
        void set_mt_sol( EventInfo *ev, Solution *sol, float *x_vector, int iz, int mtdegfree, int verbose );

/*** mtsim_subs.c ***/
        void write_mt_sim( FILE *fp, int nsol, Solution *sol, int verbose );

/*** mtsim_subs.c ***/
        void BS_replace( float *res, int rrows, float *best_b_vector, int brows, float *b_vector, int iseed, int verbose );

/*************************************************/
/*** get the total data length rows the matrix ***/
/*************************************************/

	cols = 6;
	rows = size_A_matrix( ev, grn, nsta, iz_best );

	if(verbose) 
	{
	  fprintf( stdout, "%s: invert_sim(): A matrix: rows=%d cols=%d\n", 
		progname, rows, cols);
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

/*** make the A-matrix and fill the data vector ***/

	make_amatrix( grn, ev, nsta, iz_best, a_matrix, b_vector, mtdegfree,
                        Distance_Normalize, DistNormR0, myfixisoz );

	for( k = 0; k < Nmodels; k++ )
	{

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

	} /*** loop over Nmodels bootstrap resamples ***/

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
	  fprintf( stdout, "%s: mtsim_subs.c: invert_sim_serial(): ", progname );
	  fprintf( stdout, " freeing memory inside invert\n\trows=%d\n\tcols=%d\n\n", rows, cols );
	  fflush(stdout);
	}
/***
	fprintf(stdout, "%s: freemem a_matrix\n", progname);free_matrix( a_matrix, 1, rows+1, 1, cols+1 );
	fprintf(stdout, "%s: freemem u_matrix\n", progname);free_matrix( u_matrix, 1, rows+1, 1, cols+1 );
        fprintf(stdout, "%s: freemem v_matrix\n", progname);free_matrix( v_matrix, 1, cols+1, 1, cols+1 );
        fprintf(stdout, "%s: freemem cv_matrix\n",progname);free_matrix( cv_matrix, 1, cols+1, 1, cols+1 );
        fprintf(stdout, "%s: freemem x_vector\n", progname);free_vector( x_vector, 1, cols+1 );
        fprintf(stdout, "%s: freemem s_vector\n", progname);free_vector( s_vector, 1, rows+1 );
        fprintf(stdout, "%s: freemem w_vector\n", progname);free_vector( w_vector, 1, cols+1 );
        fprintf(stdout, "%s: freemem b_vector\n", progname);free_vector( b_vector, 1, rows+1 );
        fprintf(stdout, "%s: freemem e_vector\n", progname);free_vector( e_vector, 1, cols+1 );
***/

	if(verbose)
	{
	  fprintf( stdout, "%s: Exiting from invert()\n", progname );
	  fflush(stdout);
	}

}  /*** END OF INVERT_SIM_SERIAL() ***/


/********************************************************************************************/
/*** subroutine write_mt_sim() ***/
/********************************************************************************************/

void write_mt_sim( FILE *fp, int nsol, Solution *sol, int verbose )
{
	int i;

	if( nsol == 1 )
	{

/*** k  eps   kiso  ffac   Mw   PDC   PCLVD PISO  VRED   S1    D1    R1   S2    D2     R2   lat lon ***/

	fprintf( fp, 
	  "%05d %6.3f %6.3f %8.2f %5.2f %3.0f %3.0f %3.0f %6.2f %3.0f %2.0f %4.0f %3.0f %2.0f %4.0f %.5f %.5f\n",
		nsol,
		sol->epsilon,
		sol->k,
		sol->f_factor,
		sol->mw,
		sol->PDC,
		sol->PCLVD,
		sol->PISO,
		sol->var_red,
		sol->stk0,
		sol->dip0,
		sol->rak0,
		sol->stk1,
		sol->dip1,
		sol->rak1,
		sol->lune_lat,
		sol->lune_lon );
	}
	else
	{

	  for( i = 0; i < nsol; i++ )
	  {
	    fprintf( fp,
          "%05d %6.3f %6.3f %8.2f %5.2f %3.0f %3.0f %3.0f %6.2f %3.0f %2.0f %4.0f %3.0f %2.0f %4.0f %.5f %.5f\n",
                i,
                sol[i].epsilon,
                sol[i].k,
                sol[i].f_factor,
                sol[i].mw,
                sol[i].PDC,
                sol[i].PCLVD,
                sol[i].PISO,
                sol[i].var_red,
                sol[i].stk0,
                sol[i].dip0,
                sol[i].rak0,
                sol[i].stk1,
                sol[i].dip1,
                sol[i].rak1,
                sol[i].lune_lat,
                sol[i].lune_lon );
	  }
	}

} /*** END OF WRITE_MT_SIM() ***/


/********************************************************************************************/
/*** subroutine vector_copy() ***/
/********************************************************************************************/

void vector_copy( int rows, float *out, float *in )
{
	int i;
	for( i = 0; i <= rows; i++ ) out[i] = in[i];

} /*** VECTOR_COPY() ***/


/********************************************************************************************/
/*** subroutine matrix_copy() ***/
/********************************************************************************************/

void matrix_copy( int rows, int cols, float **out, float **in )
{
	int i, j;
	for( j = 1; j <= cols; j++ )
	{
		for( i = 1; i <= rows; i++ )
		{
			out[i][j] = in[i][j];
		}
	}

} /*** MATRIX_COPY() ***/



/********************************************************************************************/
/*** subroutine inversion_init_mem() ***/
/********************************************************************************************/

void inversion_init_mem( int rows, int cols,
			float **a_matrix,
			float **u_matrix,
			float *w_vector,
			float *e_vector,
			float *x_vector,
			float **cv_matrix,
			float **v_matrix,
			float *b_vector,
			float *s_vector )
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
                e_vector[j] = 0;
                x_vector[j] = 0;
        }
        for( j=1; j<=cols; j++ )
        {
                for( i=1; i<=cols; i++ )
                {
                        cv_matrix[i][j] = 0;
                        v_matrix[i][j] = 0;
                }
        }
        for( i=1; i<=rows; i++ )
        {
                b_vector[i] = 0;
                s_vector[i] = 0;
        }

} /*** INVERSION_INIT_MEM() ***/



/********************************************************************************************/
/*** subroutine singular_value_decomposition() ***/
/********************************************************************************************/

void singular_value_decomposition( int rows, int cols,
				float **u_matrix,
				float *w_vector,
				float **v_matrix,
				float *b_vector,
				float *x_vector )
{
	int i, j;
	float wmax, wmin;
	const float tol = 1.0E-05;

	void svdcmp( float **, int, int, float *, float ** );
	void svbksb( float **, float *, float **, int, int, float *, float * );

	svdcmp( u_matrix, rows, cols, w_vector, v_matrix );

/********************************************************************************************/
/*** this is where we set the threshold for singular values allowed to be nonzero.        ***/
/*** The constant is typical but not universal. experiment with values other than 1.0E-06 ***/
/*** this will be the maximum singular value obtained                                     ***/
/********************************************************************************************/

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

} /*** END SINGULAR_VALUE_DECOMPOSITION() ***/



/********************************************************************************************/
/*** subroutine size_A_matrix() ***/
/********************************************************************************************/


int size_A_matrix( EventInfo *ev, Greens **grn, int nsta, int iz )
{
	int rows = 1;
	int ista;
	for( ista = 0; ista < nsta; ista++ )
	{
		if( ev[ista].iused == 1 )
		{
			rows += 3 * grn[ista][iz].nt;
		}
	}
	return rows;

}  /*** END SIZE_A_MATRIX() ***/


/********************************************************************************************/
/*** subroutine  set_mt_sol( ) ***/
/********************************************************************************************/

void set_mt_sol( EventInfo *ev, Solution *sol, float *x_vector, int iz, int mtdegfree, int verbose )
{
	MomentTensor Ma, Mn;

	void set_moment_tensor( MomentTensor *, MomentTensor *, float *, int, int );
	void normalize_moment_tensor( MomentTensor *, MomentTensor *, int );
	void mt2eig( MomentTensor, Solution *, int, int );
	void eig2iso( Solution *, int, int );
	void Eig2MajorDC( Solution *, int, int );
	void Eig2MinorDC( Solution *, int, int );
	void eig2lune_4mtinv( Solution *, int iz, int verbose );

	sol[iz].moment_tensor[1][1] = x_vector[1]; /* Mxx */
        sol[iz].moment_tensor[2][2] = x_vector[2]; /* Myy */
        sol[iz].moment_tensor[1][2] = x_vector[3]; /* Mxy */
        sol[iz].moment_tensor[1][3] = x_vector[4]; /* Mxz */
        sol[iz].moment_tensor[2][3] = x_vector[5]; /* Myz */
        sol[iz].moment_tensor[3][3] = x_vector[6]; /* Mzz */
        sol[iz].moment_tensor[2][1] = x_vector[3]; /* Myx */
        sol[iz].moment_tensor[3][1] = x_vector[4]; /* Mzx */
        sol[iz].moment_tensor[3][2] = x_vector[5]; /* Mzy */

	if( mtdegfree == 5 )
		sol[iz].moment_tensor[3][3] = -(x_vector[1] + x_vector[2]);

	if(verbose)
	  printf("Normalizing the seismic moment tensor\n" );

	set_moment_tensor( &Ma, &Mn, x_vector, mtdegfree, verbose );

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

/************************************************************************/
/*** eigenvalue analysis on moment tensor to get PDC and str/dip/rak  ***/
/*** for both nodal planes and T- P- and N-axis.  Clone moment tensor ***/
/*** so not to change the values by the next two routines.  Find      ***/
/*** eigenvalues and eigenvectors by decomposing Mij using eigenvalue ***/
/*** analysis                                                         ***/
/************************************************************************/
        mt2eig( Mn, sol, iz, verbose );
        eig2iso( sol, iz, verbose );
        Eig2MajorDC( sol, iz, verbose );
        Eig2MinorDC( sol, iz, verbose );
	eig2lune_4mtinv( sol, iz, verbose );

} /*** SET_MT_SOLUTION ***/



/********************************************************************************************/
/*** subroutine  fitness()   ***/
/********************************************************************************************/

void fitness( Solution *sol, int iz, float *b_vector, float *s_vector, int rows, 
	int mtdegfree, int verbose )
{
	float variance_reduction( float *, float *, int, int );
	float compute_l2norm_error( float *, float *, int );

	extern char progname[128];

	sol[iz].var_red = variance_reduction( b_vector, s_vector, 1, rows+1 );

	if(verbose) 
	{
	  printf( "%s: iz=%d %%VRED=%g\n",
              progname, iz, sol[iz].var_red );
	}

	sol[iz].l2norm_error = compute_l2norm_error( b_vector, s_vector, rows );

	if(verbose)
	{
	   printf( "%s: iz=%d %%L2NORM=%g\n",
             progname, iz,  sol[iz].l2norm_error );
	}

	sol[iz].total_fitness1 = 0;
	sol[iz].total_fitness2 = 0;

	if( mtdegfree == 1 || mtdegfree == 6 )
        {
                sol[iz].total_fitness1 = sol[iz].var_red;
                sol[iz].total_fitness2 = sol[iz].var_red;
        }
        else if( mtdegfree == 5 )
        {
           sol[iz].total_fitness1 =
                sol[iz].var_red / ( 101.0 - sol[iz].PDC );
                                                                                                                           
                                                                                                                   
           sol[iz].total_fitness2 =
                18 + sol[iz].var_red / ( 9 + (100-sol[iz].PDC) );
                                                                                                                           
                                                                                                                   
          if( sol[iz].total_fitness2 > 35 ) sol[iz].total_fitness2 = 35;
          if( sol[iz].total_fitness2 < 18 ) sol[iz].total_fitness2 = 18;
        }

} /*** END FITNESS() ***/



/********************************************************************************************/
/*** subroutine  compute_residual()   ***/
/********************************************************************************************/

void compute_residual( float *b_vector, float *s_vector, int rows, float *res_vector, float fac )
{
	int i;
	for( i = 1; i <= rows; i++ )
	{
		res_vector[i] = ( b_vector[i] - s_vector[i] );
	}

} /*** END COMPUTE_RESIDUAL() ***/


/********************************************************************************************/
/*** subroutine  dumpSAC()   ***/
/********************************************************************************************/

void dumpSAC( EventInfo *ev, Greens **grn, int nsta, int iz, int rows, float *s_vector, int verbose )
{
	int ista;
	char sac_file_name[128];
	
	void write_sac_file( char *, Sac_File *, int );
	void demultiplex( Greens **, EventInfo *, int, float *, int, int );
	void sac_minmax( float *, int, float *, float *, float * );
	
	for( ista=0; ista<nsta; ista++ )
	{
		grn[ista][iz].tra = calloc( grn[ista][iz].nt, sizeof(float) );
                grn[ista][iz].rad = calloc( grn[ista][iz].nt, sizeof(float) );
                grn[ista][iz].ver = calloc( grn[ista][iz].nt, sizeof(float) );
                ev[ista].syn_r.data = calloc( ev[ista].ns.s.npts, sizeof(float) );
                ev[ista].syn_t.data = calloc( ev[ista].ew.s.npts, sizeof(float) );
                ev[ista].syn_z.data = calloc( ev[ista].z.s.npts,  sizeof(float) );
	}

	if(verbose) 
		printf("%s: calling demultiplex\n", progname );

        demultiplex( grn, ev, iz, s_vector, rows, nsta );

	for( ista=0; ista<nsta; ista++ )
	{

       /*** SYN radial component ***/
          sprintf( sac_file_name, "%s.%s.%02d.%02d.syn.r.sac",
                ev[ista].stnm, ev[ista].net, iz, ista );
          ev[ista].syn_r.s = ev[ista].ns.s;
          sac_minmax( ev[ista].syn_r.data, ev[ista].syn_r.s.npts,
                &(ev[ista].syn_r.s.depmax), &(ev[ista].syn_r.s.depmin),
                &(ev[ista].syn_r.s.depmen) );
          write_sac_file( sac_file_name, &(ev[ista].syn_r), verbose );

        /*** SYN vertical component ***/
          sprintf( sac_file_name, "%s.%s.%02d.%02d.syn.z.sac",
                ev[ista].stnm, ev[ista].net, iz, ista );
          ev[ista].syn_z.s = ev[ista].z.s;
          sac_minmax( ev[ista].syn_z.data, ev[ista].syn_z.s.npts,
                &(ev[ista].syn_z.s.depmax), &(ev[ista].syn_z.s.depmin),
                &(ev[ista].syn_z.s.depmen) );
          write_sac_file( sac_file_name, &(ev[ista].syn_z), verbose );

        /*** SYN transverse component ***/
          sprintf( sac_file_name, "%s.%s.%02d.%02d.syn.t.sac",
                ev[ista].stnm, ev[ista].net, iz, ista );
          ev[ista].syn_t.s = ev[ista].ew.s;
          sac_minmax( ev[ista].syn_t.data, ev[ista].syn_t.s.npts,
                &(ev[ista].syn_t.s.depmax), &(ev[ista].syn_t.s.depmin),
                &(ev[ista].syn_t.s.depmen) );
          write_sac_file( sac_file_name, &(ev[ista].syn_t), verbose );

        /*** r,z,t data ***/
          sprintf( sac_file_name, "%s.%s.%02d.%02d.dat.r.sac",
                ev[ista].stnm, ev[ista].net, iz, ista );
          write_sac_file( sac_file_name, &(ev[ista].ns), verbose );
                                                                                                                                                                 
          sprintf( sac_file_name, "%s.%s.%02d.%02d.dat.z.sac",
                ev[ista].stnm, ev[ista].net, iz, ista );
          write_sac_file( sac_file_name, &(ev[ista].z), verbose );
                                                                                                                                                                 
          sprintf( sac_file_name, "%s.%s.%02d.%02d.dat.t.sac",
                ev[ista].stnm, ev[ista].net, iz, ista );
          write_sac_file( sac_file_name, &(ev[ista].ew), verbose );
	}

	for( ista=0; ista<nsta; ista++ )
	{
          free( ev[ista].syn_r.data );
          free( ev[ista].syn_r.data );
          free( ev[ista].syn_r.data );
          free( grn[ista][iz].rad );
          free( grn[ista][iz].tra );
          free( grn[ista][iz].ver );
	}

} /*** END DUMPSAC()  ***/

void Usage_Print()
{
	fprintf( stderr, "\nUSAGE: %s par= mtdegfree=(1,5,6)\n", progname );
        fprintf( stderr, "\t [no]verbose [no]dumpsac [no]fwd [no]gmtmap\n" );
        fprintf( stderr, "\t ts0=[0] fixz=[-99] [no]norm [no]shift ctol=[1] FixISOZ=[-99]\n" );

        fprintf( stderr, "\nREQUIRED PARAMETERS:\n" );
        fprintf( stderr, "\t par=(glib2inv.par) station parameter file\n" );
        fprintf( stderr, "\t mtdegfree=(1,5,6) 1=Isotropic MT, 5=Deviatoric MT, 6=Full MT\n" );

        fprintf( stderr, "\n OPTIONAL PARAMETERSL [DEFAULTS]\n" );
        fprintf( stderr, "\t ts0=[0]       Origin Time Shift Default is [0]\n" );
        fprintf( stderr, "\t fixz=[-99]    fix the depth Default is [-99] which turns off option\n" );
        fprintf( stderr, "\t [no]verbose   give verbose print to stdout. Default is off.\n" );

	fprintf( stderr, "\t [no]norm     distance normalization default is off\n" );
        fprintf( stderr, "\t R0=1.0       normalize Green functions to distance of R/R0 default is R=1 km required if norm is set\n" );

	fprintf( stderr, "\t [no]shift    shift the data automatically by cross correlation peak. default is off\n" );
        fprintf( stderr, "\t ctol=[0.6]  Correlation coefficient tolerance to shift the data when coef > ctol. defaut is off\n" );
        fprintf( stderr, "\t maxshift=   Maximum time in seconds a shift is allowed. default off\n" );

	fprintf( stderr, "\t FixISOZ=     fix the depth of the rex and zex Green's function.  Default is off\n" );
        fprintf( stderr, "\t [no]PltXcorLabel  Plot the time lag shift and cross correlation as a label in the PostScript plot [Default is on\n" );

        fprintf( stderr, "\t [no]use_snr    use peak-to-peak amplitude based Signal-Noise Ratio to make stations non-defining in inversion [default off]\n" );
        fprintf( stderr, "\t minsnr=3       minimum snr threshold.  all 3-components must be less than minsnr to set non-defining in inversion [default 3]\n" );
        fprintf( stderr, "\t               use_snr and minsnr only applies to stations that are defining and does not override users settings \n" );

        fprintf( stderr, "\n\n" );
}

/********************************************************************************************/
/*** subroutine  BS_replace()                      ***/
/********************************************************************************************/

void BS_replace( float *res, int rrows, float *best_b_vector, int brows, float *b_vector, 
	int iseed, int verbose )
{
	int i, j;
	float ratio;

	srand(iseed+1);
	ratio = 1/(float)RAND_MAX;
	rand();

	for( i = 0; i < brows; i++ )
	{
	  j = (int)rint( rrows * (ratio * (float) rand() ) );
	  b_vector[i] = best_b_vector[i] + res[j];
	}

} /*** END BS_REPLACE() ***/
