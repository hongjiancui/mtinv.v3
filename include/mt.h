#include "../include/mt_version.h"
#include "../include/sac.h"
#include "../include/mytime.h"
#include "../include/peak2peak.h"

/* Mo = math.pow( 10.0, 1.5*(Mw+10.73) ) = 1.2445146117713818e+16; Reference Mw = 0.0 */
	static float base_moment = 1.2445146117713818e+16;

typedef struct { float x, y; } Vector;
typedef struct { float s, d, r; } Plane;
typedef struct { 
	float xx,xy,xz;
	float yx,yy,yz;
	float zx,zy,zz;
} Tensor;

typedef struct {
	Vector pt;
	int type;
	char ktype[3];
	float az,pl,ev;
	float rad,str,dip;
} Poles;

typedef struct {
	float val;
	float vec[4];
	float x,y,z;
} Vector3;

typedef struct {
	float xx,xy,xz;
	float yx,yy,yz;
	float zx,zy,zz;
	Tensor T;
	float moment, Mw, abcassa;
	int expon;
	float rr,tt,ff,rt,rf,tf;
	float mt[4][4];
} MomentTensor;

typedef struct {
	Vector3 eig[4];
	Poles T, P, B;
	Plane P1, P2;
	Tensor mt;
} MTdec;

typedef struct { float r,s,d; } Polar;

typedef struct {
	int np;
	Vector *p;
	Polar *v;
	int first;
	int inside;
	int Taxis;
	int Paxis;
	int Baxis;
} NodalPlane;

#define MAX_MODEL_LAYERS 1024
typedef struct {
        char modfile[256], modpath[256];
        int nlay;
	int maxlay;
        float thick[MAX_MODEL_LAYERS];
	float  ztop[MAX_MODEL_LAYERS];
	float    vp[MAX_MODEL_LAYERS];
	float    vs[MAX_MODEL_LAYERS];
	float    qa[MAX_MODEL_LAYERS];
	float    qb[MAX_MODEL_LAYERS];
	float   rho[MAX_MODEL_LAYERS];
	float sigma[MAX_MODEL_LAYERS];
} VelMod;

/* change 2048 to 4096 */
typedef struct {
        float rss[4096], rds[4096], rdd[4096], rep[4096];
	float zss[4096], zds[4096], zdd[4096], zep[4096];
	float tss[4096], tds[4096];
} Greens_Function;

typedef struct {
	char filename[256], stnm[8], net[8];
	float stla, stlo, stel, evla, evlo, evdp;
	float rdist, az, baz;
        float t0, dt, twin, fmax, damp, eps, smin, rigidity;
	float redv, ts0, tstart, tend;

/** newly added 2010/11/27 G. Ichinose see rayp_subs.c ***/
	float Ptakeoff, Prayparameter, Pttime, Praybottom;

	int kmax, nt;
        VelMod v;
        Greens_Function g;
	float *rad, *tra, *ver;
} Greens;

typedef struct {
	int id;
	char filename[256];
	char sacpzfile[256];

        MyTime ot,a,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9;
        MyTime ref, beg, end;

	Phase pha[2];
        float P2P_snr;

	Sac_Header s;
        float  *data;

} Sac_File;

typedef struct {
        char net[8], stnm[8], modfile[256];
        char data_filename[256], glib_filename[256], ginv_filename[256];
	char comment[256];
	char wavetype[32];

	long evid,orid;
	char dbuser[32],dbsid[32];
	float mb, ms, mbmle, msvmax, mwcoda;
	int ndef;
	double time;
	float lat, lon, depth;
	int grn;
	char grname[128];
	char grname_directory[128];

        float my_z;
        float lf,hf,dt,tr,tt;
        int npole,npass;
        int nt, iused, ienvelope;
        float trbndw,a;
	float ot_shift, ts0, tstart, tend, redv;
        float str,dip,rak,Mw,Mo;
        float stla,stlo,evla,evlo,evdp;
	float rdistkm, az, baz;

/*** every station will have a total response of all comp in GS ***/
	float vred,xcor,tlag;
	int ilag;

/*** added for variance reduction by each station and each components ***/
/*** see sol[iz].var_red for event values of variance reduction ***/
	float vred_sta, vred_zcmp, vred_rcmp, vred_tcmp;

/*** added for cross_correlation.c ***/
	float zxcor,rxcor,txcor;
	float ztlag, rtlag, ttlag;
	int izlag, irlag, itlag;

	int grd_mo_type; /*** 0=disp 1=vel ***/
	float mul_factor;  /*** on-the-fly scaling of all the data ***/
	float weight;
	float time_shift_all;

	/*** mini Greens ***/
	float rdist;

	MyTime ot, ot_orig;

	/* Sac_File r,t; */

	Sac_File ns, ew, z;
	Sac_File syn_r, syn_t, syn_z;

} EventInfo;

typedef struct {
        char id[8];
        float lat, lon, elev;
        char chan[256];
        char net[8];
        char description[512];
        MyTime beg;
        MyTime end;
        float dist, az, baz;
} Station;

typedef struct {
	float moment_tensor[4][4];
	float var_red, l2norm_error, total_fitness1, total_fitness2;
	float mw, mo;
	float abcassa;
	int exponent;
	/*** normalized moment tensor in cartesian cooridinates ***/
	float mxx, myy, mzz, mxy, mxz, myz;
	double dmoment;
	float stk0,dip0,rak0;
	float stk1,dip1,rak1;
	float evlo, evla, evdp, ot;
	int mt_type;
	float epsilon;
	float k;
	float fIso, fiso, fclvd;
	float f_factor;
	float lune_lat, lune_lon;

	/*** normalized moment tensor in spherical coordinates ***/
	float mrr, mtt,mff,mrt,mrf,mtf;
	float dmrr,dmtt,dmff,dmrt,dmrf,dmtf;
	float smxx,smyy,smzz,smxy,smxz,smyz;
	float Mdc,Mclvd,Miso,Mdev,Mtotal;
	float PDC,PCLVD,PISO,PDEV;
	int nsta;
	float maxgap;
	float dminkm;
	float taz,tdp,tev,paz,pdp,pev,baz,bdp,bev;
	float stkt,stkp,stkn,plnt,plnp,plnn;
	Tensor MT;
	MomentTensor M;
	MTdec FullMT;
	MTdec Dev;
	MTdec Maj;
	MTdec Min;
} Solution;

/*** fix depth of the isotropic greens functions ***/
typedef struct {
	int iswitch; /*** this is set 0 if z=0(OFF) else z>0 set to 1(on) ***/
	int indexz;  /*** this is z[indexz] in the array ***/
	float z;     /*** input the depth in km to check against Green's func lib ***/
} FixISOZ;

/****************************************************/
/*** Dziewonski et al., 1983                      ***/
/*** Moment tensor components R=rho T=theta F=phi ***/
/*** MXX MXY MXZ     MTT -MTF  MTR                ***/
/*** MYX MYY MYZ -> -MFT  MFF -MFR                ***/
/*** MZX MZY MZZ     MRT -MRF  MRR                ***/
/****************************************************/

#define DISPLACEMENT 0
#define VELOCITY 1
#define DEVIATORIC 0
#define FULL_MT 1
#define EXPLOSION 3
#define FORCE_EXPLOSION 1
#define DEVIATORIC_MOMENT_TENSOR 5
#define FULL_MOMENT_TENSOR 6

#define FULLMT  0
#define MAJORDC 1
#define MINORDC 2

#ifndef PORTRAIT
#define PORTRAIT 1
#endif

#ifndef LANDSCAPE
#define LANDSCAPE 0
#endif

#define SCREEN_PRINT 0 /*** gs -dEPSCrop , ps2pdf -dEPSCrop , ps2raster -A ***/
#define PAPER_PRINT  1  /*** ps2pdf or lpr ***/

/*** these for grid serach routines (GS) ***/
typedef struct {
	float s,d,r;
	int id, iz;
	float Mo,Mw,z;
	Tensor M;
	float vred,xcor;
	float Mdc,Mclvd,Miso,Mtotal;
	float pdc,pclvd,piso;
} SDRvector;

typedef struct {
	float str0, str1, dstr;
	float dip0, dip1, ddip;
	float rak0, rak1, drak;
	float wdc0, wdc1, dwdc;
	float Mw0,  Mw1,  dMw;
	int nstr,ndip,nrak,nwdc,nMw;
} SDRGrid;

/*** added 08/07/2008 to swap information about the Greens function depths ***/
typedef struct {
	float zmin, zinc, zmax;
	int nz;
} Depth_Info;
