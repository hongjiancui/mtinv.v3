#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include "sac.h"

#define NUMHEADBYTES 632 /* floats and longs are 440 rest are characters */

#define SWAP   0
#define NOSWAP 1

#ifndef BIG_ENDIAN 
#define BIG_ENDIAN    1
#endif

#ifndef LITTLE_ENDIAN
#define LITTLE_ENDIAN 0
#endif

int main( int ac, char **av )
{
	int i;
	Sac_Header sp;
	float *data;
	char filename[128];
	FILE *fp;
	int check_sac_byte_order( FILE * );
	float *swap_sac_byte_order( FILE *fp, Sac_Header *sp, float *data );

	if (ac == 1)
	{
		fprintf(stderr, "Usage: %s [sacfile(s)]\n", av[0]);
		fprintf(stderr, "\n\tsacswapbytes is intelligent and nondistructive\n" );
		fprintf(stderr, "\t written by Gene Ichinose pangea@ispwest.com Oct2006\n\n");
		exit(-1);
	}
	
	for( i=1; i<ac; i++ )
	{
		strcpy( filename, av[i] );
		if( (fp = fopen(filename, "rb")) == NULL )
		{
			fprintf(stderr, "%s: Warning Cannot Opening File %s\n", av[0], filename);
		}

		data =(float *)malloc(sizeof(float));

		if( check_sac_byte_order( fp ) == SWAP ) 
		{
			data = (float *)swap_sac_byte_order( fp, &sp, data );
		}
		else
		{
			fread(&sp,sizeof(Sac_Header),1,fp);
			data =(float *)realloc( data, sp.npts * sizeof(float) );
			fread(&data[0], sp.npts*sizeof(float), 1, fp );
		}
		fclose(fp);
		sscanf( sp.kstnm, "%s", sp.kstnm );
		sscanf( sp.kcmpnm, "%s", sp.kcmpnm );

		printf( "sta=%s cmp=%s nt=%d dt=%g b=%g min=%g max=%g\n",
			sp.kstnm, sp.kcmpnm, sp.npts, sp.delta, sp.b, sp.depmin, sp.depmax );

	/** write a new sac file **/
		fp = fopen( filename, "wb" );
		fwrite(&sp, sizeof(Sac_Header), 1, fp);
		fwrite(&data[0], sp.npts*sizeof(float), 1, fp);

		/* if( malloc_size( data ) > 0 ) free(data); */
		if( data != NULL ) free(data);
	}
}

float *swap_sac_byte_order( FILE *fp, Sac_Header *sp, float *data )
{
	int j;
	char cbuf[NUMHEADBYTES];

	float float_swap();
	int int_swap();
	long  long_swap();

/* set some sac header defaults */
	*sp = sac_null;
	
	rewind(fp);
/** read in header **/
	if( fread(cbuf, 440*sizeof(char), 1, fp) <= 0 )
	{
		fprintf(stderr, "A read error has occurred\n" );
		exit(-1);
	}

	sp->delta 	= float_swap(cbuf);
	sp->depmin	= float_swap(cbuf+4);
	sp->depmax	= float_swap(cbuf+8);
	sp->scale	= float_swap(cbuf+12);
	sp->odelta	= float_swap(cbuf+16);
	sp->b    	= float_swap(cbuf+20);
	sp->e   	= float_swap(cbuf+24);
	sp->o    	= float_swap(cbuf+28);
	sp->a    	= float_swap(cbuf+32);
	sp->internal1	= float_swap(cbuf+36);
	sp->t0    	= float_swap(cbuf+40);
	sp->t1    	= float_swap(cbuf+44);
	sp->t2    	= float_swap(cbuf+48);
	sp->t3    	= float_swap(cbuf+52);
	sp->t4    	= float_swap(cbuf+56);
	sp->t5    	= float_swap(cbuf+60);
	sp->t6    	= float_swap(cbuf+64);
	sp->t7    	= float_swap(cbuf+68);
	sp->t8    	= float_swap(cbuf+72);
	sp->t9    	= float_swap(cbuf+76);
	sp->f     	= float_swap(cbuf+80);
	sp->resp0	= float_swap(cbuf+84);
	sp->resp1	= float_swap(cbuf+88);
	sp->resp2	= float_swap(cbuf+92);
	sp->resp3	= float_swap(cbuf+96);
	sp->resp4	= float_swap(cbuf+100);
	sp->resp5	= float_swap(cbuf+104);
	sp->resp6	= float_swap(cbuf+108);
	sp->resp7	= float_swap(cbuf+112);
	sp->resp8	= float_swap(cbuf+116);
	sp->resp9	= float_swap(cbuf+120);
	sp->stla 	= float_swap(cbuf+124);
	sp->stlo 	= float_swap(cbuf+128);
	sp->stel 	= float_swap(cbuf+132);
	sp->stdp 	= float_swap(cbuf+136);
	sp->evla  	= float_swap(cbuf+140);
	sp->evlo  	= float_swap(cbuf+144);
	sp->evel 	= float_swap(cbuf+148);
	sp->evdp  	= float_swap(cbuf+152);
	sp->unused1	= float_swap(cbuf+156);
	sp->user0	= float_swap(cbuf+160);
	sp->user1	= float_swap(cbuf+164);
	sp->user2	= float_swap(cbuf+168);
	sp->user3	= float_swap(cbuf+172);
	sp->user4	= float_swap(cbuf+176);
	sp->user5	= float_swap(cbuf+180);
	sp->user6	= float_swap(cbuf+184);
	sp->user7	= float_swap(cbuf+188);
	sp->user8	= float_swap(cbuf+192);
	sp->user9	= float_swap(cbuf+196);
	sp->dist         = float_swap(cbuf+200);
	sp->az           = float_swap(cbuf+204);
	sp->baz          = float_swap(cbuf+208);
	sp->gcarc        = float_swap(cbuf+212);
	sp->internal2    = float_swap(cbuf+216);
	sp->internal3    = float_swap(cbuf+220);
	sp->depmen       = float_swap(cbuf+224);
	sp->cmpaz 	= float_swap(cbuf+228);
	sp->cmpinc	= float_swap(cbuf+232);
	sp->unused2	= float_swap(cbuf+236);
	sp->unused3	= float_swap(cbuf+240);
	sp->unused4	= float_swap(cbuf+244);
	sp->unused5	= float_swap(cbuf+248);
	sp->unused6	= float_swap(cbuf+252);
	sp->unused7	= float_swap(cbuf+256);
	sp->unused8	= float_swap(cbuf+260);
	sp->unused9	= float_swap(cbuf+264);
	sp->unused10	= float_swap(cbuf+268);
	sp->unused11	= float_swap(cbuf+272);
	sp->unused12	= float_swap(cbuf+276);

	sp->nzyear	= int_swap(cbuf+280);
	sp->nzjday	= int_swap(cbuf+284);
	sp->nzhour	= int_swap(cbuf+288);
	sp->nzmin 	= int_swap(cbuf+292);
	sp->nzsec 	= int_swap(cbuf+296);
	sp->nzmsec	= int_swap(cbuf+300);
	sp->internal4	= int_swap(cbuf+304);
	sp->internal5	= int_swap(cbuf+308);
	sp->internal6	= int_swap(cbuf+312);
	sp->npts         = int_swap(cbuf+316);
	sp->internal7    = int_swap(cbuf+320);
	sp->internal8    = int_swap(cbuf+324);
	sp->unused13     = int_swap(cbuf+328);
	sp->unused14     = int_swap(cbuf+332);
	sp->unused15     = int_swap(cbuf+336);
	sp->iftype       = int_swap(cbuf+340);
	sp->idep         = int_swap(cbuf+344);
	sp->iztype       = int_swap(cbuf+348);
	sp->unused16     = int_swap(cbuf+352);
	sp->iinst        = int_swap(cbuf+356);
	sp->istreg       = int_swap(cbuf+360);
	sp->ievreg       = int_swap(cbuf+364);
	sp->ievtyp       = int_swap(cbuf+368);
	sp->iqual        = int_swap(cbuf+372);
	sp->isynth       = int_swap(cbuf+376);
	sp->unused17     = int_swap(cbuf+380);
	sp->unused18     = int_swap(cbuf+384);
	sp->unused19     = int_swap(cbuf+388);
	sp->unused20     = int_swap(cbuf+392);
	sp->unused21     = int_swap(cbuf+396);
	sp->unused22     = int_swap(cbuf+400);
	sp->unused23     = int_swap(cbuf+404);
	sp->unused24     = int_swap(cbuf+408);
	sp->unused25     = int_swap(cbuf+412);
	sp->unused26     = int_swap(cbuf+416);
	sp->leven        = int_swap(cbuf+420);
	sp->lpspol       = int_swap(cbuf+424);
	sp->lovrok       = int_swap(cbuf+428);
	sp->lcalda       = int_swap(cbuf+432);
	sp->unused27     = int_swap(cbuf+436);

	fread(cbuf, 8*sizeof(char), 1, fp);
	bcopy( cbuf, sp->kstnm, 8 );

	fread(cbuf, 16*sizeof(char), 1, fp);
	bcopy( cbuf, sp->kevnm, 16 );

	fread(cbuf, 8*sizeof(char), 1, fp);
	bcopy( cbuf, sp->khole, 8 );

	fread(cbuf, 8*sizeof(char), 1, fp);
	bcopy( cbuf, sp->ko, 8 );

	fread(cbuf, 8*sizeof(char), 1, fp);
	bcopy( cbuf, sp->ka, 8 );

	fread(cbuf, 8*sizeof(char), 1, fp);
	bcopy( cbuf, sp->kt0, 8 );

	fread(cbuf, 8*sizeof(char), 1, fp);
	bcopy( cbuf, sp->kt1, 8 );

	fread(cbuf, 8*sizeof(char), 1, fp);
	bcopy( cbuf, sp->kt2, 8 );

	fread(cbuf, 8*sizeof(char), 1, fp);
	bcopy( cbuf, sp->kt3, 8 );

	fread(cbuf, 8*sizeof(char), 1, fp);
	bcopy( cbuf, sp->kt4, 8 );

	fread(cbuf, 8*sizeof(char), 1, fp);
	bcopy( cbuf, sp->kt5, 8 );
                                                                                                                             
	fread(cbuf, 8*sizeof(char), 1, fp);
	bcopy( cbuf, sp->kt6, 8 );
                                                                                                                             
	fread(cbuf, 8*sizeof(char), 1, fp);
	bcopy( cbuf, sp->kt7, 8 );
                                                                                                                             
	fread(cbuf, 8*sizeof(char), 1, fp);
	bcopy( cbuf, sp->kt8, 8 );
                                                                                                                             
	fread(cbuf, 8*sizeof(char), 1, fp);
	bcopy( cbuf, sp->kt9, 8 );
                                                                                                                             
	fread(cbuf, 8*sizeof(char), 1, fp);
	bcopy( cbuf, sp->kf, 8 );
                                                                                                                             
	fread(cbuf, 8*sizeof(char), 1, fp);
	bcopy( cbuf, sp->kuser0, 8 );
                                                                                                                             
	fread(cbuf, 8*sizeof(char), 1, fp);
	bcopy( cbuf, sp->kuser1, 8 );
                                                                                                                             
	fread(cbuf, 8*sizeof(char), 1, fp);
	bcopy( cbuf, sp->kuser2, 8 );
                                                                                                                             
	fread(cbuf, 8*sizeof(char), 1, fp);
	bcopy( cbuf, sp->kcmpnm, 8 );

	fread(cbuf, 8*sizeof(char), 1, fp);
	bcopy( cbuf, sp->knetwk, 8 );
                                                                                                                             
	fread(cbuf, 8*sizeof(char), 1, fp);
	bcopy( cbuf, sp->kdatrd, 8 );
                                                                                                                             
	fread(cbuf, 8*sizeof(char), 1, fp);
	bcopy( cbuf, sp->kinst, 8 );

/** read the data **/
	data = realloc( data, sp->npts * sizeof(float) );
	for ( j=0; j<sp->npts; j++)
	{
		fread(cbuf, sizeof(char), 4, fp);
		data[j] = float_swap(cbuf);
	}
	return((float *)data);
}

int check_sac_byte_order( fp )
FILE *fp;
{
	char cbuf[4];
	float test0, test1;
	float float_swap(), float_noswap();
	int My_Endianness;
	int TestByteOrder(void);

	if( TestByteOrder() == BIG_ENDIAN )
	{
		My_Endianness = BIG_ENDIAN;
		printf("I am Big Endian... ");
	}
	else if( TestByteOrder() == LITTLE_ENDIAN )
	{
		My_Endianness = LITTLE_ENDIAN;
		printf("I am Little Endian... ");
	}

	fread(cbuf, 4*sizeof(char), 1, fp);
	test0  = float_noswap(cbuf);
	test1  = float_swap(cbuf);

	if( test0 > 1.0E-08 && test0 < 1.0E+08 ) 
	{
		if( My_Endianness == BIG_ENDIAN )
			printf("SAC file is big endian and no byte swapping needed.\n");
		else
			printf("SAC file is little endian and no byte swapping needed.\n");

		rewind(fp);
		return NOSWAP;
	}
	else if( test1 > 1.0E-08 && test1 < 1.0E+08 )
	{
		if( My_Endianness == BIG_ENDIAN )
			printf("SAC file is little endian need to swap bytes\n");
		else
			printf("SAC file is big endian need to swap bytes\n");
		rewind(fp);
		return SWAP;
	}
	else
	{	
		printf("Sorry... cant tell what byte order?\n");
		exit(-1);
	}
}

float float_noswap( cbuf )
char cbuf[];
{
        union {
                char cval[4];
                float fval;
        } f_union;
        f_union.cval[0] = cbuf[0];
        f_union.cval[1] = cbuf[1];
        f_union.cval[2] = cbuf[2];
        f_union.cval[3] = cbuf[3];
        return(f_union.fval);
}

int int_swap( cbuf )
char cbuf[];
{
        union {
                char cval[4];
                int ival;
        } int_union;
        int_union.cval[3] = cbuf[0];
        int_union.cval[2] = cbuf[1];
        int_union.cval[1] = cbuf[2];
        int_union.cval[0] = cbuf[3];
        return(int_union.ival);
}
                                                                                                                          
                                                                                                                          
long long_swap( cbuf )
char cbuf[];
{
        union {
                char cval[4];
                long lval;
        } l_union;
                                                                                                                                     
        l_union.cval[3] = cbuf[0];
        l_union.cval[2] = cbuf[1];
        l_union.cval[1] = cbuf[2];
        l_union.cval[0] = cbuf[3];
        return(l_union.lval);
}
                                                                                                                                     
float float_swap( cbuf )
char cbuf[];
{
        union {
                char cval[4];
                float fval;
        } f_union;
                                                                                                                                     
        f_union.cval[3] = cbuf[0];
        f_union.cval[2] = cbuf[1];
        f_union.cval[1] = cbuf[2];
        f_union.cval[0] = cbuf[3];
        return(f_union.fval);
}

int TestByteOrder(void)
{
        short int word = 0x0001;
        char *byte = (char *) &word;
        return(byte[0] ? LITTLE_ENDIAN : BIG_ENDIAN);
}
