#include <stdio.h>
#include <stdlib.h>
#include "sac.h"

float *readsac( Sac_Header *s, char *filename )
{
        FILE *fp;
        float *data;
        if( (fp = fopen( filename, "rb" )) == NULL )
        {
                fprintf(stderr, "error opening file %s\n", filename );
                exit(-1);
        }
        fread( s, sizeof(Sac_Header), 1, fp );
        data = (float *)calloc(s->npts, sizeof(float) );
        fread( &data[0], s->npts * sizeof(float), 1, fp );
        fclose(fp);
        return (float *)data;
}

void wrtoldsac( char *FO, Sac_Header *s, float *data )
{
	FILE           *iopo;
	if ((iopo = fopen(FO, "w")) == NULL)
	{
		printf("cannot open output file %s\n", FO);
		exit(1);
	}
	fwrite(s, sizeof(Sac_Header), 1, iopo);
	fwrite(&data[0], s->npts*sizeof(float), 1, iopo);
	fclose(iopo);
	return;
}

void wrtnewsac( char *FO, float dt, int ns, float *ar, float b)
{
    Sac_Header      sp;
    float           zero = 0.0;
    FILE           *iopo;
 
    if ((iopo = fopen(FO, "w")) == NULL)
    {
        printf("cannot open output file %s\n", FO);
        exit(1);
    }

/* convert parameters to sac header parameters */
    sp = sac_null;
 
/* floats */
    sp.b = b;
    sp.delta = dt;  /* sampling rate in seconds per sample */
    sp.e = b + (ns * sp.delta);
 
/* integers */
    sp.nvhdr  = 6;           /* LLL sets it */
    sp.norid  = 0;           /* LLL sets it */
    sp.nevid  = 0;           /* LLL sets it */
    sp.iftype = ITIME;  /* data type: IRLIM spec re im | IAMPH amp ph | IXY general x,y  */
    sp.idep   = IUNKN;  /* not disp vel acc or volts */
    sp.iztype = IB;  /* types: IUNKN,IB,IDAY,IO,IA,ITn  */
    sp.ievtyp = IUNKN;  /* type of event IQUAKE - earthquake */
    sp.npts   = ns;
/* logicals  */
    sp.leven = TRUE;        /* is data evenly sampled  */
    sp.lpspol = FALSE;          /* LLL sets it */
    sp.lcalda = TRUE;           /* should az,baz,gcarc,dist be calculated?  */
    sp.lhdr5  = FALSE;        /* LLL sets it */
 
/* write binary sac header */
    if (fwrite(&sp, sizeof(Sac_Header), 1, iopo) != 1)
    {
        perror("fwrite");
        printf("FATAL ERROR: bad  write %s, header\n", FO);
        exit(4);
    }
/* write floating point data */
    if (fwrite(&ar[0], ns * sizeof(float), 1, iopo) != 1)
    {
        printf("ns=%d, %lu, nbytes=%lu\n", ns, (unsigned long)ns, (unsigned long)ns * sizeof(float));
        printf("FATAL ERROR: bad write %s\n", FO);
        exit(5);
    }
    fclose(iopo);
    return;
}

void set_sac_minmax( Sac_Header *s, float *data )
{
	int nt, i;
	float min,max,sum=0;
	nt = s->npts;
	min = data[0];
	max = data[0];
	sum = data[0];
	for( i=1; i<nt; i++ )
	{
		if( data[i] < min ) min = data[i];
		if( data[i] > max ) max = data[i];
		sum += data[i];
	}
	s->depmax = max;
	s->depmin = min;
	s->depmen = sum/nt;
	return;
}
