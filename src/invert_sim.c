
void invert_sim( EventInfo *ev, 
		Greens **grn, 
		int nsta, 
		int nz, 
		Solution *sol,
		int verbose,
		int mtdegfree,
		int Distance_Normalize,
		FixISOZ myfixisoz,
		char *output_pathname )
{

	int ista, iz;
	char sac_file_name[128];
	float wmax, wmin;
	MomentTensor Ma, Mn;
	void sac_absmax( float *, int, float * );

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
	float *res_vector;   /*** residual vector ***/
	
	float **matrix( int, int, int, int );
	float *vector( int, int );
	float variance_reduction( float *, float *, int, int );
	float compute_l2norm_error( float *, float *, int );
	void demultiplex( Greens **, EventInfo *, int, float *, int );
	void set_moment_tensor( MomentTensor *, MomentTensor *, float *, int, int );
	void normalize_moment_tensor( MomentTensor *, MomentTensor *, int );
	void calc_deviatoric_moment( MomentTensor *, int );
	void mt2eig( MomentTensor, Solution *, int, int );
	void eig2iso( Solution *, int, int );
	void Eig2MajorDC( Solution *, int, int );
	void Eig2MinorDC( Solution *, int, int );
	void make_amatrix( Greens **, EventInfo *, int, int, float **, float *, int, int, FixISOZ );
	void write_gmt_xy_values( char *, Solution *, EventInfo *, Greens **, int, int, int );
	void write_sac_file( char *, Sac_File *, int );
	void svbksb( float **, float *, float **, int, int, float *, float * );
	void svdcmp( float **, int, int, float *, float ** );
	void svdvar( float **, int, float *, float ** );
	void compute_residual( float *b_vector, float *s_vector, int rows, float *res_vector );
	void wrtnewsac( char *FO, float dt, int ns, float *ar, float b);
	void matmul( int, float **, int, float *, int, float * );
	void sac_minmax( float *, int, float *, float *, float * );

/**********************************************************************************
	MomentTensor Ma_err, Mn_err;	
	float max, tmp;
	float data_mean, data_variance, residual_variance;
	float mean( float *, int );
	float variance( float *, int, float );
	float root_mean_ssquare_variance( float *, float *, int );
	void diag( int, float **, float * );
*********************************************************************************/

/*************************************************************/
/*** set the iz to the iz_best passed to invert_sim via nz ***/
/*************************************************************/

	iz = nz;
	nz = 1;

/*************************************************/
/*** get the total data length rows the matrix ***/
/*************************************************/
	cols = 6;
	rows = 1;
	for(ista=0; ista<nsta; ista++ )
	{
		if( ev[ista].iused == 1 )
		{
			rows += 3*grn[ista][iz].nt;
		}
	}
	printf("%s: invert_sim(): A matrix: rows=%d cols=%d\n", progname, rows, cols);

/*************************/
/*** Allocating memory ***/
/*************************/
	printf("%s: allocating memory for iz=%d\n", progname, iz );
	fprintf(stdout, "%s: Allocating memory for a_matrix rows=%d cols=%d\n", progname, rows, cols );
	a_matrix  = matrix( 1, rows+1, 1, cols+1 );
	fprintf(stdout, "%s: Allocating memory for u_matrix\n", progname );
	u_matrix  = matrix( 1, rows+1, 1, cols+1 );
	fprintf(stdout, "%s: Allocating memory for v_matrix\n", progname );
	v_matrix  = matrix( 1, cols+1, 1, cols+1 );
	fprintf(stdout, "%s: Allocating memory for cv_matrix\n", progname );
	cv_matrix = matrix( 1, cols+1, 1, cols+1 );
	fprintf(stdout, "%s: Allocating memory for w_vector\n", progname );
	w_vector = vector( 1, cols+1 );
	fprintf(stdout, "%s: Allocating memory for e_vector\n", progname );
	e_vector = vector( 1, cols+1 );
	fprintf(stdout, "%s: Allocating memory for x_vector\n", progname );
	x_vector = vector( 1, cols+1 );
	fprintf(stdout, "%s: Allocating memory for b_vector\n", progname );
	b_vector = vector( 1, rows+1 );
	fprintf(stdout, "%s: Allocating memory for s_vector\n", progname );
	s_vector = vector( 1, rows+1 );
	res_vector = vector( 1, rows+1 );

/********************************/
/***** initialize memory ********/
/********************************/
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

/*********************************************************/
/*** set up the A matrix and data vector for inversion ***/
/*********************************************************/
	make_amatrix( grn, ev, nsta, iz, a_matrix, b_vector, mtdegfree, Distance_Normalize, myfixisoz );
	
/***************************************/
/*** copy A matrix into the U matrix ***/
/***************************************/
	for( j=1; j<=cols; j++ )
	{
		for( i=1; i<=rows; i++ )
		{
			u_matrix[i][j] = a_matrix[i][j];
		}
	}

/******************************/
/***** do the inversion *******/
/******************************/
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

/***************************************************************************************/
/*** multiply solution x vector with a matrix of greens function to get the s vector ***/
/*** or synthetics vector then compute the variance reduction and l2_norm error      ***/
/***************************************************************************************/

	matmul( 0, a_matrix, cols, x_vector, rows, s_vector );

	sol[iz].var_red = variance_reduction( b_vector, s_vector, 1, rows+1 );

	if(verbose) 
	{
		printf( "%s: iz=%d %%VRED=%g\n",
			progname, iz, sol[iz].var_red );
	}

	compute_residual( b_vector, s_vector, rows, res_vector );
	wrtnewsac( "res.sac", 1.0, rows, res_vector, 0.0);
	wrtnewsac( "dat.sac", 1.0, rows, b_vector, 0.0);
	wrtnewsac( "syn.sac", 1.0, rows, s_vector, 0.0);

/****************************************************/
/*** calculate the covariance matrix              ***/
/*** default is sigma = 1                         ***/
/*** get sigma from RMS preevent noise ? level    ***/
/*** error is 1.96 * sqrt( diag(CV_matrix) )      ***/
/****************************************************/
                                                                                                                                                                      
        svdvar( v_matrix, cols, w_vector, cv_matrix );
                                                                                                                                                                      
/******************************************************************/
/*** form the moment tensor from solution vector x              ***/
/*** col1-Mxx, col2-Myy, col3-Mxy, col4-Mxz, col5-Myz, col6-Mzz ***/
/******************************************************************/
 
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

	set_moment_tensor( &Ma, &Mn, x_vector, mtdegfree, verbose );
	normalize_moment_tensor( &Ma, &Mn, verbose );
	sol[iz].dmoment = Ma.moment;
        sol[iz].mw      = Ma.Mw;
        sol[iz].exponent = Ma.expon;
        sol[iz].abcassa = Ma.abcassa;

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

	mt2eig( Mn, sol, iz, verbose );
	eig2iso( sol, iz, verbose );
	Eig2MajorDC( sol, iz, verbose );
	Eig2MinorDC( sol, iz, verbose );

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

	sol[iz].evlo = ev[0].evlo;
        sol[iz].evla = ev[0].evla;
        sol[iz].evdp = ev[0].evdp;
        sol[iz].ot   = ev[0].ot_shift;

	write_gmt_xy_values( "./", sol, ev, grn, iz, nsta, verbose );

/****************************/
/*** free memory clean up ***/
/****************************/
	fprintf( stdout, "%s: freeing memory inside invert\n\trows=%d\n\tcols=%d\n\n",
		progname, rows, cols );
/*
	fprintf(stdout, "%s: freeing memory for a_matrix\n", progname);
	free_matrix( a_matrix, 1, rows+1, 1, cols+1 );
*/

	fprintf(stdout, "%s: freeing memory for u_matrix\n", progname);
	free_matrix( u_matrix, 1, rows+1, 1, cols+1 );
                                                                                                                                                                      
        fprintf(stdout, "%s: freeing memory for v_matrix\n", progname); free_matrix( v_matrix, 1, cols+1, 1, cols+1 );
        fprintf(stdout, "%s: freeing memory for cv_matrix\n", progname);free_matrix( cv_matrix, 1, cols+1, 1, cols+1 );
        fprintf(stdout, "%s: freeing memory for x_vector\n", progname); free_vector( x_vector, 1, cols+1 );
        fprintf(stdout, "%s: freeing memory for s_vector\n", progname); free_vector( s_vector, 1, rows+1 );
        fprintf(stdout, "%s: freeing memory for w_vector\n", progname); free_vector( w_vector, 1, cols+1 );
        fprintf(stdout, "%s: freeing memory for b_vector\n", progname); free_vector( b_vector, 1, rows+1 );
        fprintf(stdout, "%s: freeing memory for e_vector\n", progname); free_vector( e_vector, 1, cols+1 );
	fprintf( stdout, "%s: Exiting from invert()\n", progname );
	fflush( stdout );

}
