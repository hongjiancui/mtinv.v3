#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char progname[128];

int main( int ac, char **av )
{
	int i, j;
	int num, modulus, nloop = 126;

	char color[64][42] = { "black", "red", "blue", "green", "brown", "purple", 
				"cyan", "yellow", "pink", "gray", "orange", "navy",
				"lavender", "royalblue", "skyblue", "turquoise", "seagreen",
				"khaki", "goldenrod", "salmon", "maroon", "black", "red",
				"blue",	"green", "brown", "purple", "cyan", "yellow",
				"pink", "gray", "orange", "khaki", "goldenrod", "salmon", "maroon" };

	strcpy( progname, av[0] );
	if( ac == 3 )
	{
		modulus = atoi( av[1] );
		num     = atoi( av[2] );
		
		for( j = 0, i = 0; i < nloop; i++ )
		{
			if( ( i % modulus ) == 0 ) j++;

			if( i == num - 1 )
			{
				fprintf( stdout, "%s", color[j] );
			}
		}
	}
	else	
	{
		fprintf( stderr, "%s: error needs 2 args: modulus num\n", progname );
		exit(-1);
	}
	exit(0);
}
