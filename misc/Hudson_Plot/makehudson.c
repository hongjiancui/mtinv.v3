#include <stdio.h>
#include <stdlib.h>

int main(int ac,char **av)
{
	char region[64];
	char proj[64];
	char yshift[32];
	char xshift[32];
	char lineweight[32];
	char psfile[128];	
	char orientation[8];
	int init = 0;
	FILE *fp;

	int setpar(int,char **);
	int mstpar();
	int getpar();
	void endpar();

	setpar(ac,av);
	mstpar("R", "s", region);
	mstpar("J", "s", proj);
	mstpar("Y", "s", yshift);
	mstpar("X", "s", xshift);
	mstpar("W", "s", lineweight);
	mstpar("ps", "s", psfile );
	getpar("init", "b", &init );
	mstpar("orientation", "s", orientation );
	endpar();

	fp = fopen("tmp.csh","w");
	fprintf( fp, "#!/bin/csh\n" );

	if( init )
	{
	  fprintf( fp, "psxy %s %s -M %s %s -K %s %s >! %s << EOF\n",
	    region, proj, lineweight, orientation, xshift, yshift, psfile );
	}
	else
	{
	  fprintf( fp, "psxy %s %s -M %s %s -O -K %s %s >> %s << EOF\n",
	    region, proj, lineweight, orientation, xshift, yshift, psfile );
	}

	fprintf( fp, ">\n" );
	fprintf( fp, "  0.000000   1.000000\n" );
	fprintf( fp, "  1.333333   0.333333\n" );
	fprintf( fp, "  0.000000  -1.000000\n" );
	fprintf( fp, " -1.333333  -0.333333\n" );
	fprintf( fp, "  0.000000   1.000000\n" );
	fprintf( fp, ">\n" );
	fprintf( fp, "  0.000000   1.000000\n" );
	fprintf( fp, "  0.000000  -1.000000\n" );
	fprintf( fp, ">\n" );
	fprintf( fp, " -1.000000   0.000000\n" );
	fprintf( fp, "  1.000000   0.000000\n" );
	fprintf( fp, ">\n" );
	fprintf( fp, " -0.500000   0.500000\n" );
	fprintf( fp, "  0.000000   0.500000\n" );
	fprintf( fp, "  0.666667   0.666667\n" );
	fprintf( fp, ">\n" );
	fprintf( fp, " -0.666667  -0.666667\n" );
	fprintf( fp, "  0.000000  -0.500000\n" );
	fprintf( fp, "  0.500000  -0.500000\n" );
	fprintf( fp, ">\n" );
	fprintf( fp, "  0.000000   1.000000\n" );
	fprintf( fp, "  0.571429   0.142857\n" );
	fprintf( fp, "  0.000000  -1.000000\n" );
	fprintf( fp, ">\n" );
	fprintf( fp, "  0.000000   1.000000\n" );
	fprintf( fp, " -0.571429  -0.142857\n" );
	fprintf( fp, "  0.000000  -1.000000\n" );
	fprintf( fp, "EOF\n" );

	fprintf( fp, "pstext %s %s -N -O -K >> %s << EOF\n", region, proj, psfile );
	fprintf( fp, " -0.030000   1.030000 12 0 31 BR EXPLOSION\n" );
	fprintf( fp, " -0.474444   0.585556 12 0 31 BC +TC       \n" );
	fprintf( fp, " -0.696667   0.373333 12 0 31 BC +VD       \n" );
	fprintf( fp, " -1.000000  -0.030000 12 0 31 BR +CLVD     \n" );
	fprintf( fp, "  0.030000   0.030000 12 0 31 BL DC       \n" );
	fprintf( fp, "  1.030000   0.000000 12 0 31 TL -CLVD    \n" );
	fprintf( fp, "  0.696667  -0.333333 12 0 31 TL -VD      \n" );
	fprintf( fp, "  0.474444  -0.555556 12 0 31 TL -TC      \n" );
	fprintf( fp, "  0.030000  -1.000000 12 0 31 TL IMPLOSION\n" );
	fprintf( fp, "EOF\n" );

	fprintf( fp, "psxy %s %s -Sc0.05i -W1p/0/0/0 -G0/0/0 -O -P -K >> %s << EOF\n",
		region, proj, psfile );
	fprintf( fp, "  0.000000   1.000000 12 0 31 BR EXPLOSION\n" );
        fprintf( fp, " -0.444444   0.555556 12 0 31 BC TC       \n" );
        fprintf( fp, " -0.666667   0.333333 12 0 31 BC VD       \n" );
        fprintf( fp, " -1.000000  -0.000000 12 0 31 BR CLVD     \n" );
        fprintf( fp, "  0.000000   0.000000 12 0 31 BL DC       \n" );
        fprintf( fp, "  1.000000  -0.000000 12 0 31 TL -CLVD    \n" );
        fprintf( fp, "  0.666667  -0.333333 12 0 31 TL -VD      \n" );
        fprintf( fp, "  0.444444  -0.555556 12 0 31 TL -TC      \n" );
        fprintf( fp, "  0.000000  -1.000000 12 0 31 TL IMPLOSION\n" );
        fprintf( fp, "EOF\n" );

	fclose(fp);
	system("csh ./tmp.csh");
}
