#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <math.h>

char progname[128];

typedef struct {

	int npts;

	float min_T, max_T, ave_T, adev_T, sdev_T;
        float min_k, max_k, ave_k, adev_k, sdev_k;
	float ave_lune_lat;
	float ave_lune_lon;
	float angular_stddev;

	float *eta;
	float *T;
	float *kiso;
	float *T1;
	float *k1;

	float *lune_lat;
	float *lune_lon;
	float *lune_lat_mu;
	float *lune_lon_mu;

} SrcType;

int main( int ac, char **av )
{
	SrcType st;
	int i, rows,cols=2;
	int c_rows, c_cols;
	float Nmodels = 1000;

	int ihudson = 0;
	int verbose = 0;
	int ilune   = 1;

	float major_axes, minor_axes, major_ax_Az, minor_ax_Az;

	float **A_matrix;
	float **C_matrix;
	float **A_transpose_matrix;

/**************************************************/
/*** subroutines: Function Interface Prototypes ***/
/**************************************************/

	void readInput( SrcType *st, char *inputfilename, int verbose );

	void transform_to_srctype_space( SrcType *st );

	void writeoutTEST( SrcType *st );

	void minmax2( float *data, int n, float *min, float *max );

	void moment2( float *data, int n, float *ave, float *adev, float *sdev );

	void matmul2(	float **A_transpose_matrix, int Atcols, int Atrows, 
			float **A_matrix, int Arows, int Acols, 
			float **C_matrix, int *c_rows, int *c_cols );

	void print_matrix( float **A_matrix, int rows, int cols, char *label );

	void scale_matrix( float **c, int rows, int cols, float scale );

	void transpose_matrix( int rows, int cols, float **in, float **out );

	void eigenvalues( float **C_matrix, int rows, int cols,
		float *major_axes, float *minor_axes,
		float *major_ax_Az, float *minor_ax_Az );

	void sphere_mean( int npts, float *lat, float *lon, float *avglat, float *avglon, 
		float *angular_standard_devation, int verbose );

	float **matrix( int start_row, int stop_row, int start_col, int stop_col );

	float *vector( int start, int stop );

	int setpar( int ac, char **av );
	int getpar(), mstpar();
	void endpar();

/**************************************************/
/*** start main program, get command line args  ***/
/**************************************************/
	strcpy( progname, av[0] );
	fprintf( stdout, "%s: ilune=%d\n", progname, ilune );

	setpar( ac,av);
	getpar( "hudson",  "b", &ihudson);
	getpar( "verbose", "b", &verbose);
	getpar( "lune",    "b", &ilune);
	endpar();

/**************************************************************/
/*** allocate memory for srctype struction arrarys          ***/
/*** read in file mtsim.out, translate epsilon to T = 2*eta ***/
/**************************************************************/

	st.eta  = calloc( 1001, sizeof(float) );
	st.T    = calloc( 1001, sizeof(float) );
	st.kiso = calloc( 1001, sizeof(float) );
	st.lune_lat = calloc( 1001, sizeof(float) );
	st.lune_lon = calloc( 1001, sizeof(float) );

	readInput( &st, "mtsim.out", verbose );

	if(verbose) writeoutTEST( &st );

	if(ihudson)
	{
		transform_to_srctype_space( &st );
	}

/******************************************************/
/*** calculate the mean, min and max for kiso and T ***/
/******************************************************/

	sphere_mean( st.npts, &(st.lune_lat[0]), &(st.lune_lon[0]), 
		&st.ave_lune_lat, &st.ave_lune_lon, &st.angular_stddev, verbose );

 	minmax2( &(st.kiso[0]), st.npts, &(st.min_k), &(st.max_k) );
	minmax2( &(st.T[0]),    st.npts, &(st.min_T), &(st.max_T) );

	moment2( &(st.kiso[0]), st.npts, &(st.ave_k), &(st.adev_k), &(st.sdev_k) );
	moment2( &(st.T[0]),    st.npts, &(st.ave_T), &(st.adev_T), &(st.sdev_T) );

	if(verbose) 
	{
	  fprintf( stdout, "npts=%d\n", st.npts );
	  fprintf( stdout, "ave  = %.2f %.2f \n", st.ave_T,  st.ave_k );
	  fprintf( stdout, "min  = %.2f %.2f \n", st.min_T,  st.min_k );
          fprintf( stdout, "max  = %.2f %.2f \n", st.max_T,  st.max_k );
	  fprintf( stdout, "adev = %.2f %.2f \n", st.adev_T, st.adev_k );
	  fprintf( stdout, "sdev = %.2f %.2f \n", st.sdev_T, st.sdev_k );
	  fprintf( stdout, "ave_lune_lat = %g avg_lune_lon = %g angular_standard_deviation = %g\n",
		st.ave_lune_lat, st.ave_lune_lon, st.angular_stddev );
	  
	  fflush(stdout);
	}

/***********************************/
/*** remove mean from kiso and T ***/
/***********************************/

	st.k1 = calloc( st.npts+1, sizeof(float) );
	st.T1 = calloc( st.npts+1, sizeof(float) );

	for( i = 1; i <= st.npts; i++ )
	{
		st.k1[i] = st.kiso[i] - st.ave_k;
		st.T1[i] = st.T[i]    - st.ave_T;
	}

/*******************************************/
/*** remove mean from lune_lat, lune_lon ***/
/*******************************************/

	st.lune_lat_mu = calloc( st.npts+1, sizeof(float) );
	st.lune_lon_mu = calloc( st.npts+1, sizeof(float) );

	for( i = 1; i <= st.npts; i++ )
	{
		st.lune_lat_mu[i] = st.lune_lat[i] - st.ave_lune_lat;
		st.lune_lon_mu[i] = st.lune_lon[i] - st.ave_lune_lon;

	/***
		fprintf( stdout, "i=%d lune_lat_mu=%g lune_lon_mu=%g ave_lune_lat=%g avg_lune_lon=%g %g %g\n",
			i,
			st.lune_lat_mu[i],
			st.lune_lon_mu[i],
			st.ave_lune_lat,
			st.ave_lune_lon,
			st.lune_lat[i],
			st.lune_lon[i] );
	***/

	}

/**************************************************/
/*** create A matrix A = [ T' kiso' ]           ***/
/**************************************************/

	rows = st.npts;
	A_matrix           = matrix( 0, rows+1, 0, cols+1 );
	A_transpose_matrix = matrix( 0, cols+1, 0, rows+1 );
	C_matrix           = matrix( 0, cols+1, 0, cols+1 );

	for( i = 1; i <= rows; i++ )
	{
		if( ilune )
		{
		  A_matrix[i][1] = st.lune_lon_mu[i];
		  A_matrix[i][2] = st.lune_lat_mu[i];
		}
		else
		{
		  A_matrix[i][1] = st.T1[i];
		  A_matrix[i][2] = st.k1[i];
		}
	}

	if(verbose) print_matrix( A_matrix, rows, cols, "A_matrix" );

	transpose_matrix( rows, cols, A_matrix, A_transpose_matrix );

	if(verbose) print_matrix( A_transpose_matrix, cols, rows, "A_transpose_matrix" );

/*****************************************************/
/*** create covariance matrix C = (A' * A)/Nmodels ***/
/*****************************************************/

	matmul2( A_transpose_matrix, cols, rows, A_matrix, rows, cols, C_matrix, &c_rows, &c_cols );

	Nmodels = (float)rows;

	scale_matrix( C_matrix, c_rows, c_cols, (float)Nmodels );

	if(verbose) print_matrix( C_matrix, c_rows, c_cols, "C_matrix" ); 

/**************************************************/
/*** calculate the eigenvalues and eigenvectors ***/
/**************************************************/

	eigenvalues( C_matrix, c_rows, c_cols,
		&major_axes, &minor_axes, &major_ax_Az, &minor_ax_Az );

/******************************/
/*** write out for GMT psxy ***/
/******************************/

	if( ilune )
	{
	  fprintf( stdout, "ellipse x0=%.4f y0=%.4f theta=%.2f a=%g b=%g\n",
                st.ave_lune_lon, 
		st.ave_lune_lat, 
		major_ax_Az,
                major_axes, 
		minor_axes );
	}
	else
	{
	  if(ihudson)
	  {
	    fprintf( stdout, "ellipse x0=%.2f y0=%.2f theta=%.0f a=%g b=%g\n",
		st.ave_T, 
		st.ave_k, 
		major_ax_Az, 
		major_axes, 
		minor_axes );
	  }
	  else
	  {
	    fprintf( stdout, "ellipse hudson x0=%.2f y0=%.2f theta=%.0f a=%g b=%g\n",
                st.ave_T, 
		st.ave_k, 
		major_ax_Az,
                major_axes, 
		minor_axes );
	  }
	}

/******************************/
/*** for GMT psxy -SE       ***/
/******************************/
	if(ilune)
	{
		fprintf( stdout, "%.4f %.4f %.1f %g %g\n",
                	st.ave_lune_lon, 
			st.ave_lune_lat, 
			major_ax_Az,
                	major_axes, 
			minor_axes );
	}
	else
	{
		fprintf( stdout, "%.2f %.2f %.0f %g %g\n",
			st.ave_T, 	
			st.ave_k, 
			major_ax_Az,
                	major_axes, 
			minor_axes );
	}

} /*** END MAIN() ***/

/***************************************************************/
/*** compute eigenvalues given 2X2 covariance matrix         ***/
/***************************************************************/

void eigenvalues( float **C_matrix, int rows, int cols,
	float *major_axes, float *minor_axes, 
	float *major_ax_Az, float *minor_ax_Az )
{
	float *eval1, *eval2;
        float **z1, **z2;
	int i, j, k;
	float R2D, majaxAz, minaxAz, minax, majax;
	int verbose = 0;

	float **matrix( int, int, int, int );
        float *vector( int, int );
 	void matrix_copy( int, int, float **, float ** );
        void vector_copy( int, float *, float * );
        void tred2( float **, int, float *, float * );
        void tqli( float *, float *, int, float ** );
	void eigsrt( float *, float **, int );
	float phase(float);

	if( rows != cols )
	{
		fprintf(stdout, "matrix must be square\n");
		exit(-1);
	}
	R2D = 180/M_PI;

	eval1 = vector( 0, rows+1 );
	eval2 = vector( 0, rows+1 );
	z1 = matrix( 0, rows+1, 0, rows+1 );
	z2 = matrix( 0, rows+1, 0, rows+1 );

	matrix_copy( rows, cols, z1, C_matrix );
	matrix_copy( rows, cols, z2, C_matrix );

	tred2( z2, rows, eval1, eval2 );
	tqli( eval1, eval2, rows, z2 );
	eigsrt( eval1, z2, rows );

	if(verbose)
	{
	  for( i=1; i<=rows; i++ )
	  {
		printf("eval1=%+7.3e z2=", eval1[i]);
		for( j=1; j<=cols; j++)
		{
			printf("%6.2f ", z2[i][j] );
		}
		printf("\n" );
	  }
	}

/***************************************************************/
/*** compute the angles for the major and minor axes         ***/
/***************************************************************/

	majaxAz = atan2( z2[1][1], z2[2][1] ) * R2D;
	minaxAz = atan2( z2[1][2], z2[2][2] ) * R2D;
	majaxAz = phase( majaxAz );
	minaxAz = phase( minaxAz );
	if(verbose) printf( "majaxAz=%.2f minaxAz=%.2f\n", majaxAz, minaxAz );
	*major_ax_Az = majaxAz;
	*minor_ax_Az = minaxAz;

/***************************************************************/
/** convert the eigenvalues to major and minor axes lengths  ***/
/** for 95% confidence ellipse                               ***/
/***************************************************************/

	majax = sqrt(eval1[1]) * 2.447; 
	minax = sqrt(eval1[2]) * 2.447;
	if(verbose) printf( "majax=%+7.3e minax=%+7.3e\n", majax, minax );
	*major_axes = majax;
	*minor_axes = minax;

} /*** END EIGENVALUES()  ***/

/***************************************************************/
/*** unwrap phase ***/
/***************************************************************/

float phase( float a )
{
	if( a < 0 )
	{
		while( a < 0 ) a += 360;
	}
	else if( a >= 360 )
	{
		while( a >= 360 ) a -= 360;
	}
	return a;

}  /*** END PHASE() ***/

/***************************************************************/
/*** print a matrix to the screen ***/
/***************************************************************/

void print_matrix( float **a, int rows, int cols, char *label )
{
	int i, j;

	printf( "print_matrix %s rows=%d cols=%d\n", label, rows, cols );
	fflush(stdout);

	for( i = 1; i <= rows; i++ )
	{
		printf( "%05d ", i );
		for( j = 1; j <= cols; j++ )
		{
			printf( "%+7.3e ", a[i][j] );
		}
		printf( "\n" );
	}
	fflush(stdout);

} /*** END print_matrix() ***/

/***************************************************************/
/*** multiply a scalar to a matrix ***/
/***************************************************************/

void scale_matrix( float **c, int rows, int cols, float scale )
{
	int i, j;
	float tmp;
	for( i = 1; i <= cols; i++ )
	{
		for( j = 1; j <= cols; j++ )
		{
			tmp = c[i][j];
			c[i][j] = tmp/scale;
		}
	}

} /*** END  scale_matrix() ***/

/***************************************************************/
/*** multiply two matrix together ***/
/***************************************************************/

void matmul2( float **m1, int m, int n, float **m2, 
	int p, int q, float **result, int *rows, int *cols )
{
	int i, j, k;
	if( n != p )
	{
	  fprintf( stderr, "error in matmul2 n=%d not equal to p=%d\n",
		n, p );
	  exit(-1);
	}

	*rows = m;
	*cols = q;

	for( i = 1; i <= m; i++ )
	{
	  for( j = 1; j <= q; j++ )
	  {
		result[i][j] = 0;
		for( k = 1; k <= n; k++ )
		{
		  result[i][j] += ( m1[i][k] * m2[k][j] );
		}
	  }
	}

} /*** matmul2() ***/

/***************************************************************/
/*** transpose matrix ***/
/***************************************************************/

void transpose_matrix( int rows, int cols, float **in, float **out )
{
	int i, j;
	for( j = 1; j <= cols; j++ )
	{
		for( i = 1; i <= rows; i++ )
		{
			out[j][i] = in[i][j];
		}
	}
}

/***************************************************************/
/*** read input from file   mtsim.out format from MTSIM.C   ***/
/***************************************************************/

void readInput( SrcType *st, char *inputfilename, int verbose )
{
	FILE *fp;
	char rec[256];
	int i = 1;

/*** dummy variables ***/

	int idum;
	float f_factor, mw, PDC, PCLVD, PISO, var_red, stk0, stk1, dip0, dip1, rak0, rak1;

	if( (fp = fopen( inputfilename, "r" )) == NULL )
	{
		fprintf( stderr, "ERROR cannot open file %s\n",
			inputfilename );
		exit(-1);
	}

	while( fgets( rec, 256, fp ) != NULL )
	{
		st->eta  = realloc( st->eta,  (i+1)*sizeof(float) );
		st->kiso = realloc( st->kiso, (i+1)*sizeof(float) );

		st->lune_lat = realloc( st->lune_lat, (i+1)*sizeof(float) );
		st->lune_lon = realloc( st->lune_lon, (i+1)*sizeof(float) );

		/***          1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 ***/
		sscanf( rec, "%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
			&idum,		/* 1 */
			&(st->eta[i]),	/* 2 */
			&(st->kiso[i]),	/* 3 */
			&f_factor,	/* 4 */
			&mw,		/* 5 */
			&PDC,		/* 6 */
			&PCLVD,		/* 7 */
			&PISO,		/* 8 */
			&var_red,	/* 9 */
			&stk0,		/* 10 */
			&dip0,		/* 11 */
			&rak0,		/* 12 */
			&stk1,		/* 13 */
			&dip1,		/* 14 */
			&rak1,		/* 15 */
			&(st->lune_lat[i]),	/* 16 */
			&(st->lune_lon[i]) );	/* 17 */

		 if(verbose)
		 {
		   fprintf( stdout, "%d %g %g %g %g\n",
			i, st->eta[i], st->kiso[i], st->lune_lat[i], st->lune_lon[i] );
		 }

		i++;
	}
	fclose(fp);
	st->npts = i - 1;

	if(verbose) fprintf( stdout, "npts read %d\n", st->npts );
	st->T = realloc( st->T, (st->npts + 1) * sizeof(float) );

/*** convert here from epsilon to T=2*eps ***/

	for( i = 1; i <= st->npts; i++ )
		st->T[i] = st->eta[i];
	/*
		st->T[i] = 2 * st->eta[i];
	*/

} /*** END readInput() ***/

/***************************************************************/
/*** a wrapper to transform t,k in linear space to src type space ***/
/***************************************************************/

void transform_to_srctype_space( SrcType *st )
{
	int i;
	float t, k, x = 0, y = 0;

	int hudson_transform(float,float,float*,float*);

	for( i = 1; i <= st->npts; i++ )
	{
		t = st->T[i];
		k = st->kiso[i];
		if( hudson_transform(t,k,&x,&y) != 0 )
                {
                  fprintf(stderr, 
  "ERROR, Hudson_transform row %d out of range, skipping row.\n", i );
                        continue;
                }
		st->T[i]    = x;
		st->kiso[i] = y;
	}

} /*** END transform_to_srctype_space() ***/

/***************************************************************/
/*** print readInput to the screen ***/
/***************************************************************/

void writeoutTEST( SrcType *st )
{
	int i;
	fprintf( stdout, "npts = %d\n", st->npts );
	for( i = 1; i <= st->npts; i++ )
	{
		fprintf( stdout, "%05d %.2f %.2f %.2f %.4f %.4f\n",
			i,
			st->eta[i],	
			st->kiso[i],	
			st->T[i],
			st->lune_lat[i],
			st->lune_lon[i] );
	}
	fflush(stdout);
}


/***************************************************************/
/*** simple minimum and maximum for a vector ***/
/***************************************************************/

void minmax2( float *data, int n, float *min, float *max )
{
	int i;
	float xmin,xmax;
	xmin = data[1];
	xmax = data[1];
	for( i = 1; i <= n; i++ )
	{
		if( data[i] > xmax ) xmax = data[i];
		if( data[i] < xmin ) xmin = data[i];
	}
	*min = xmin;
	*max = xmax;

} /*** END MINMAX2() ***/

/*************************************************************/
/*** this version of moment only calculates average,      ***/
/*** average deviation, standard deviation                ***/ 
/************************************************************/

void moment2( float *data, int n, float *ave, float *adev, float *sdev )
{
        int j;
        float s,p;
	float svar, skew, curt;

        if (n <= 1)
	{
	  printf("n must be at least 2 in MOMENT\n");
	  exit(-1);
	}
        s=0.0;
        for (j=1;j<=n;j++) s += data[j];
        *ave=s/n;
        *adev=(svar)=(skew)=(curt)=0.0;
        for (j=1;j<=n;j++) 
	{
                *adev += fabs(s=data[j]-(*ave));
                svar += (p=s*s);
                skew += (p *= s);
                curt += (p *= s);
        }
        *adev /= n;
        svar /= (n-1);
        *sdev=sqrt(svar);
        if (svar) 
	{
                skew /= (n*(svar)*(*sdev));
                curt=(curt)/(n*(svar)*(svar))-3.0;
        } 
	else 
	{
	  printf("No skew/kurtosis when variance = 0 (in MOMENT)\n");
	  exit(-1);
	}

}  /*** END MOMENT2() ***/

void vector_copy( int rows, float *out, float *in )
{
        int i;
        for( i = 0; i <= rows; i++ ) out[i] = in[i];
                                                                                                                                                                 
} /*** END VECTOR_COPY() ***/

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
                                                                                                                                                                 
} /*** END MATRIX_COPY() ***/



/********************************************************************************************/
/*** sphere_mean() - calculates the mean latitude and longitude on a sphere               ***/
/********************************************************************************************/

void sphere_mean( int npts, float *lat, float *lon, float *avglat, float *avglon, 
	float *angular_standard_deviation, int verbose )
{
	int i;
	float colat, phi, theta;
	float rs, xs, ys, zs;
	float xsum = 0, ysum = 0, zsum = 0;

	for ( i = 0; i < npts; i++ )
	{
		colat = 90 - lat[i];
		phi   = ( M_PI / 180 ) * colat;
		theta = ( M_PI / 180 ) * lon[i];
		xsum += cos( theta ) * sin( phi );
		ysum += sin( theta ) * sin( phi );
		zsum += cos( phi );
	}

	rs = sqrt( xsum*xsum + ysum*ysum + zsum*zsum );
	xs = xsum / rs;
	ys = ysum / rs;
	zs = zsum / rs;
	
	*avglon = atan2( ys, xs ) * ( 180 / M_PI );
	*avglat = 90 - ( acos( zs ) * ( 180 / M_PI ) );
	*angular_standard_deviation = acos( rs / npts ) * ( 180 / M_PI );

	if(verbose)
	{
	  fprintf( stdout, 
		"%s: npts=%d avglat=%g avglon=%g angular_standard_deviation=%g\n",
			progname, 
			npts,
			*avglat,
			*avglon,
			*angular_standard_deviation );
	}

} /*** END SPHERE_MEAN() ***/
