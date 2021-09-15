#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/nrutil.h"     /** numerical recipes **/
#include "../include/mt.h"         /** global datatype and structure declarations **/

char progname[128];

/*** this program inputs

     Total Moment
     DC     - S/D/R   %Mo
     CLVD   - Vertical, Horizontal, %Mo
     ISO    - %Mo

     output moment tensor
            moment tensor decomposition
            GMT plot deviatoric only
****/

int main( int ac, char **av )
{
/*** local variables and function prototypes ***/
	Tensor ISO_MT, CLVD_MT, DC_MT, MT;
	float MoTotal, str, dip, rak, pdc;
	char clvd_type[6];
	float pclvd, piso;

	void computeFullMT( Tensor isoM, Tensor dcM, Tensor clvdM, Tensor *M );
	void create_CLVD( float Mo, float pclvd, char *clvdType, Tensor *M );
	void create_DC( float Mo, float pdc, float str, float dip, float rak, Tensor *M );
	void create_ISO( float Mo, float piso, Tensor *M );
	void writeTensor( Tensor *M, char *label );

/*** variables and functions from mtinv toolkit libglib.a and mtinv_subs.o ***/
	MomentTensor Ma, Mn;
        float x_vector[7];
        int mtdegfree = 6;
        int verbose = 0;
        Solution *sol;
        int iz;
/* see include/mt.h Mo = math.pow( 10.0, 1.5*(Mw+10.73) ) = 1.2445146117713818e+16; Reference Mw = 0.0 */

	void set_moment_tensor( MomentTensor *Ma, float *x, int idf, int verbose );
	void normalize_moment_tensor( MomentTensor *Ma, MomentTensor *Mn, float Mtotal, int verbose );
	void mt2eig( MomentTensor Mn, Solution *sol, int iz, int verbose );
	void eig2iso( Solution *sol, int iz, int verbose );
	void Eig2MajorDC( Solution *sol, int iz, int verbose );
	void Eig2MinorDC( Solution *sol, int iz, int verbose );
	void eig2lune_4mtinv( Solution *sol, int iz, int verbose );

/*** misc ***/
	int setpar(int,char **), mstpar(), getpar();	
	void endpar();

/**********************/
/**** begin program ***/
/**********************/
	strcpy( progname, av[0] );

	setpar(ac,av);
	mstpar("Mo", "f", &MoTotal);

/*** getpar input Double Couple ***/
	mstpar("str", "f", &str );
	mstpar("dip", "f", &dip );
	mstpar("rak", "f", &rak );
	mstpar("pdc", "f", &pdc );

/*** getpar input Compensated Linear Vector Dipole ***/
/*** horizontal1, horizontal2 or vertical - negative ( -h1, -h2, -v ) ***/
/*** horizontal1, horizontal2 or vertical - positive ( +h1, +h2, +v ) ***/

	mstpar("clvd_type", "s", &clvd_type );
	mstpar("pclvd", "f", &pclvd );

/*** GetPar Input Isotropic ***/
	mstpar("piso", "f", &piso );

	getpar("verbose","b",&verbose);
	endpar();

	fprintf( stdout, "%s: Mo=%e str=%g dip=%g rak=%g pdc=%g clvd_type=%s pclvd=%g piso=%g\n",
		progname, MoTotal, str, dip, rak, pdc, clvd_type, pclvd, piso );

/*** create the moment tensor ***/
	
	create_ISO( MoTotal, piso, &ISO_MT );
	writeTensor( &ISO_MT, "isotropic" );

	create_DC( MoTotal, pdc, str, dip, rak, &DC_MT );
	writeTensor( &DC_MT, "double-couple" );

	create_CLVD( MoTotal, pclvd, clvd_type, &CLVD_MT );
	writeTensor( &CLVD_MT, "compensated linear vector dipole" );

	computeFullMT( ISO_MT, DC_MT, CLVD_MT, &MT );
	writeTensor( &MT, "full - moment tensor" );

/*** set and normalize the moment tensor ***/

	x_vector[1] = MT.xx;
	x_vector[2] = MT.yy;
	x_vector[3] = MT.xy;
	x_vector[4] = MT.xz;
	x_vector[5] = MT.yz;
	x_vector[6] = MT.zz;

	mtdegfree = 6;
	set_moment_tensor( &Ma, x_vector, mtdegfree, verbose );
	
	sol = (Solution *)malloc(2*sizeof(Solution));
	iz = 0;
	sol[iz].mt_type  = FULL_MT;

/*** decompose the moment tensor ***/

	mt2eig( Ma, sol, iz, verbose );
	eig2iso( sol, iz, verbose );

/*** normalize all moment tensor and eigenvalues ***/

	normalize_moment_tensor( &Ma, &Mn, sol[iz].Mtotal, verbose );
	sol[iz].dmoment  = Ma.moment;
        sol[iz].mw       = Ma.Mw;
        sol[iz].exponent = Ma.expon;
        sol[iz].abcassa  = Ma.abcassa;

	sol[iz].FullMT.eig[1].val /= pow(10.0, sol[iz].exponent );
        sol[iz].FullMT.eig[2].val /= pow(10.0, sol[iz].exponent );
        sol[iz].FullMT.eig[3].val /= pow(10.0, sol[iz].exponent );

        sol[iz].Dev.eig[1].val /= pow(10.0, sol[iz].exponent );
        sol[iz].Dev.eig[2].val /= pow(10.0, sol[iz].exponent );
        sol[iz].Dev.eig[3].val /= pow(10.0, sol[iz].exponent );

        sol[iz].FullMT.T.ev /= pow(10.0, sol[iz].exponent );
        sol[iz].FullMT.B.ev /= pow(10.0, sol[iz].exponent );
        sol[iz].FullMT.P.ev /= pow(10.0, sol[iz].exponent );

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

	Eig2MajorDC( sol, iz, verbose );
	Eig2MinorDC( sol, iz, verbose );
	eig2lune_4mtinv( sol, iz, verbose );
}

void computeFullMT( Tensor isoM, Tensor dcM, Tensor clvdM, Tensor *M )
{
	M->xx = isoM.xx + dcM.xx + clvdM.xx;
	M->yy = isoM.yy + dcM.yy + clvdM.yy;
	M->zz = isoM.zz + dcM.zz + clvdM.zz;
	M->xy = isoM.xy + dcM.xy + clvdM.xy;
	M->xz = isoM.xz + dcM.xz + clvdM.xz;
	M->yz = isoM.yz + dcM.yz + clvdM.yz;
	M->yx = M->xy;
	M->zx = M->xz;
	M->zy = M->yz;
}

/*** horizontal1, horizontal2 or vertical - negative ( -h1, -h2, -v ) ***/
/*** -2  0  0      +1  0  0      +1  0  0
      0 +1  0       0 -2  0       0 +1  0
      0  0 +1       0  0 +1       0  0 -2  ***/
/*** horizontal1, horizontal2 or vertical - positive ( +h1, +h2, +v ) ***/
/***
     +2  0  0      -1  0  0      -1  0  0
      0 -1  0       0 +2  0       0 -1  0
      0  0 -1       0  0 -1       0  0 +2  ***/

void create_CLVD( float Mo, float pclvd, char *clvdType, Tensor *M )
{
	float factor = 0.5;

	if( strcmp(clvdType,"+v") == 0 )
	{
		M->xx = -1.0 * Mo * factor * pclvd;
		M->yy = -1.0 * Mo * factor * pclvd;
		M->zz = +2.0 * Mo * factor * pclvd;
	}
	else if( strcmp(clvdType,"+h2") == 0 )	
	{
		M->xx = -1.0 * Mo * factor * pclvd;
                M->yy = +2.0 * Mo * factor * pclvd;
                M->zz = -1.0 * Mo * factor * pclvd;
	}
	else if( strcmp(clvdType,"+h1") == 0 )
	{
		M->xx = +2.0 * Mo * factor * pclvd;
                M->yy = -1.0 * Mo * factor * pclvd;
                M->zz = -1.0 * Mo * factor * pclvd;
	}
	else if( strcmp(clvdType,"-v") == 0 )
        {
		M->xx = +1.0 * Mo * factor * pclvd;
                M->yy = +1.0 * Mo * factor * pclvd;
                M->zz = -2.0 * Mo * factor * pclvd;
        }
        else if( strcmp(clvdType,"-h2") == 0 )  
        {
		M->xx = +1.0 * Mo * factor * pclvd;
                M->yy = -2.0 * Mo * factor * pclvd;
                M->zz = +1.0 * Mo * factor * pclvd;
        }
        else if( strcmp(clvdType,"-h1") == 0 )
        {
		M->xx = -2.0 * Mo * factor * pclvd;
                M->yy = +1.0 * Mo * factor * pclvd;
                M->zz = +1.0 * Mo * factor * pclvd;
        }
	M->xy = 0;
        M->xz = 0;
        M->yz = 0;
        M->yx = 0;
        M->zx = 0;
        M->zy = 0;
}

void create_DC( float Mo, float pdc, float str, float dip, float rak, Tensor *M )
{
	float d2r, strr, dipr, rakr, tol = 1.0E-07;
	void scaleTensor( Tensor *M, float scale );
	void writeTensor( Tensor *M, char *label );
	void floorTensor( Tensor *M, float tol );

	d2r   = M_PI / 180.0;
	strr  = str * d2r;
	dipr  = dip * d2r;
	rakr  = rak * d2r;

	M->xx = -(sin(dipr)*cos(rakr)*sin(2*strr)+sin(2*dipr)*sin(rakr)*sin(strr)*sin(strr));
        M->yy =  (sin(dipr)*cos(rakr)*sin(2*strr)-sin(2*dipr)*sin(rakr)*cos(strr)*cos(strr));
        M->zz =  (sin(2*dipr)*sin(rakr));
        M->xy =  (sin(dipr)*cos(rakr)*cos(2*strr)+0.5*sin(2*dipr)*sin(rakr)*sin(2*strr));
        M->xz = -(cos(dipr)*cos(rakr)*cos(strr)+cos(2*dipr)*sin(rakr)*sin(strr));
        M->yz = -(cos(dipr)*cos(rakr)*sin(strr)-cos(2*dipr)*sin(rakr)*cos(strr));
        M->yx = M->xy;
        M->zx = M->xz;
        M->zy = M->yz;

	/* writeTensor( M, "DC before scaling" ); */
	floorTensor( M, tol );
	scaleTensor( M, (Mo*pdc) );
}

void floorTensor( Tensor *M, float tol )
{
	if( fabs(M->xx) < tol ) M->xx = 0;
	if( fabs(M->xy) < tol ) M->xy = 0;
	if( fabs(M->xz) < tol ) M->xz = 0;
	if( fabs(M->yy) < tol ) M->yy = 0;
	if( fabs(M->yz) < tol ) M->yz = 0;
	if( fabs(M->zz) < tol ) M->zz = 0;
	M->yx = M->xy;
	M->zx = M->xz;
	M->zy = M->yz;
}

void scaleTensor( Tensor *M, float scale )
{
	M->xx *= scale;
	M->yy *= scale;
        M->zz *= scale;
        M->xy *= scale;
        M->xz *= scale;
        M->yz *= scale;
        M->yx = M->xy;
        M->zx = M->xz;
        M->zy = M->yz;
}

void create_ISO( float Mo, float piso, Tensor *M )
{
	/* float factor = 0.333333333333333333; */
	float factor = 1.0;

	M->xx = Mo * factor * piso;
	M->yy = Mo * factor * piso;
	M->zz = Mo * factor * piso;
	M->xy = 0;
	M->xz = 0;
	M->yz = 0;
	M->yx = 0;
	M->zx = 0;
	M->zy = 0;
}

void writeTensor( Tensor *M, char *label )
{
	fprintf( stdout, "%s\n", label );
	fprintf( stdout, "\t %6.2e %6.2e %6.2e\n", M->xx, M->xy, M->xz );
	fprintf( stdout, "\t %6.2e %6.2e %6.2e\n", M->yx, M->yy, M->yz );
	fprintf( stdout, "\t %6.2e %6.2e %6.2e\n", M->zx, M->zy, M->zz );
	fprintf( stdout, "\n" );
}
