#include <stdio.h>
#include <stdlib.h>

int main(int ac, char **av)
{
	int iz, nz;
	float *zvec;
	char zstr[8];

	int setpar(int ac, char **av);
	int mstpar();
	void endpar();

	setpar( ac, av );

	mstpar( "nz", "d", &nz );

	zvec = (float *)calloc( nz, sizeof(float) ); 

	sprintf( zstr, "vf[%d]", nz );

	fprintf( stdout, "zstr=%s nz=%d\n", zstr, nz );

	mstpar( "z", zstr, zvec );

	endpar();

	for( iz = 0; iz < nz; iz++ )
	{
		fprintf( stdout,  "iz=%d z=%g\n", iz, zvec[iz] );
	}
}
